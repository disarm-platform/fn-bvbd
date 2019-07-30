import json
import urllib.request
import geopandas
import pandas as pd
import numpy as np
import disarm_gears
from datetime import timedelta
from scipy.spatial import distance_matrix


def run_function(params: dict):
    #
    # 1. Handle input
    #

    # Set random seed
    np.random.seed(1000)

    end_date = params.get('end_date')
    observed_periods = params.get('observed_periods')
    point_data = params.get('point_data')
    layer_names = params.get('layer_names')

    #
    # 2. Process
    #

    # Load stored data
    link_1 = 'https://www.dropbox.com/s/ab7i0ynbn2ei77z/tha_attrib_covs.json?dl=1'
    with urllib.request.urlopen(link_1) as url:
        tha_attrib = pd.read_json(json.loads(url.read().decode()))
    tha_attrib.dropna(axis=0, inplace=True)  # Some covariates are NA
    tha_attrib.reset_index(inplace=True, drop=True)

    # Define prediction grid
    pred_grid = tha_attrib.copy()
    pred_grid['total_cases'] = np.nan
    pred_grid['incidence'] = np.nan
    pred_grid['prob_cases'] = np.nan
    pred_grid['risk_class'] = np.nan

    # Add covariates to the attributes table
    # At the moment the covariates are already loaded through tha_attrib_covs.csv

    # Define timeline and spatiotemporal grid
    timeline = disarm_gears.frames.Timeframe(start=None, end=end_date, length=observed_periods, by='day', step=28)
    input_data = geopandas.GeoDataFrame.from_features(point_data['features'])
    input_data['knot'] = timeline.which_knot(np.array(input_data.date))
    input_data = input_data.loc[input_data.knot > -1, :]
    input_data['village_id'] = -1
    input_data['lng'] = input_data.geometry.x
    input_data['lat'] = input_data.geometry.y
    input_xy = np.vstack([input_data.geometry.x, input_data.geometry.y]).T
    input_data.drop(labels=['geometry'], axis=1, inplace=True)

    # Define polygons of observed cases
    G = disarm_gears.frames.SparseGrid(x_lim=(97, 106), y_lim=(5, 21), n_cols=9, n_rows=16, tag_prefix='G')
    G.add_polygon_from_xy(input_xy)
    H = G.get_simplified()
    h_frame = disarm_gears.frames.TilePattern(geometries=H.geometry)
    input_data['cluster_id'] = h_frame.locate(np.array(input_data[['lng', 'lat']]))
    ix_cluster = h_frame.locate(np.array(tha_attrib[['lng', 'lat']]))
    tha_attrib['cluster_id'] = ix_cluster
    tha_attrib = tha_attrib.loc[tha_attrib.cluster_id > -1, :]

    # Filter data per cluster and fit model
    for cluster_i in range(H.shape[0]):

        input_data_i = input_data.loc[input_data.cluster_id == cluster_i].copy()
        polyg_data_i = tha_attrib.loc[tha_attrib.cluster_id == cluster_i].copy()
        polyg_series = pd.concat([polyg_data_i] * observed_periods)
        polyg_series['knot'] = np.hstack([np.repeat(ki, polyg_data_i.shape[0]) for ki in np.arange(observed_periods)])
        polyg_series['total_cases'] = 0

        dmat = distance_matrix(input_data_i[['lng', 'lat']], polyg_data_i[['lng', 'lat']])
        input_data_i.loc[input_data_i.index, 'village_id'] = polyg_data_i.index[
            np.array([di.argmin() for di in dmat])].copy()

        # Add input data to spatiotemporal grid
        input_data_i = input_data_i.groupby(by=['village_id', 'knot'],
                                            as_index=False).sum()[['village_id', 'knot', 'total_cases']]

        for i, rowi in input_data_i.iterrows():
            polyg_series.loc[np.logical_and(polyg_series.index == rowi.village_id,
                                            polyg_series.knot == rowi.knot), 'total_cases'] += rowi.total_cases
        pred_grid_i = polyg_data_i.copy()
        pred_grid_i['knot'] = observed_periods

        # Candidate GAM models
        gam_flist = []
        gam_sp = f"total_cases ~ offset(log(population)) + te(lng, lat, bs='gp', m=list(c(2, -1, 2)), d=2, k=25)"
        gam_flist.append(gam_sp)
        if layer_names is not None:
            gam_spc = [gam_sp] + [f'{i}' for i in layer_names]
            gam_spc = '+'.join(gam_spc)
            gam_flist.append(gam_spc)
        rho1 = 1.5
        gam_st = f"total_cases ~ offset(log(population)) + te(lng, lat, knot, bs='gp', m=list(c(2, 2, 2), c(2, {rho1}, 2)), d=c(2, 1), k=c(25, 3))"
        gam_flist.append(gam_st)
        if layer_names is not None:
            gam_stc = [gam_st] + [f'{i}' for i in layer_names]
            gam_stc = '+'.join(gam_stc)
            gam_flist.append(gam_stc)


        # Fit candidate models and keep the best
        for j, gam_formula in enumerate(gam_flist):
            gam_j = disarm_gears.r_plugins.r_methods.mgcv_fit(data=polyg_series, formula=gam_formula,
                                                              family='poisson', weights=None, method='REML',
                                                              bam=False)
            y_j = disarm_gears.r_plugins.r_methods.mgcv_predict(gam_j, data=pred_grid_i, response_type='response')
            phi_j = 10000 * y_j / pred_grid_i.population
            ix_ubre = disarm_gears.r_plugins.r_methods.get_names(gam_j).index('gcv.ubre')
            gcv = gam_j[ix_ubre][0]
            if j == 0:
                gam_star = gam_j
                gcv_star = gcv
                y_star = y_j.copy()
                phi_star = phi_j.copy()
            elif gcv_star < gcv:
                gcv_star = gcv
                gam_star = gam_j
                y_star = y_j.copy()
                phi_star = phi_j.copy()

        n_samples = 500
        y_sims = disarm_gears.r_plugins.r_methods.mgcv_posterior_samples(gam_star, data=pred_grid_i,
                                                                         n_samples=n_samples, response_type='link')
        prob_cases = np.sum(y_sims + np.log(pred_grid_i.population[None, :]) >= 0, axis=0) / n_samples
        risk_class = np.zeros_like(prob_cases)
        for i, pi in enumerate([.1, .25, .5, .75, .975]):
            risk_class[prob_cases >= pi] = i + 1

        pred_grid.loc[pred_grid_i.index, 'total_cases'] = y_star
        pred_grid.loc[pred_grid_i.index, 'incidence'] = phi_star
        pred_grid.loc[pred_grid_i.index, 'prob_cases'] = prob_cases
        pred_grid.loc[pred_grid_i.index, 'risk_class'] = risk_class


    #
    # 3. Package output
    #

    pred_date = timeline.end + timedelta(28)
    pred_grid['date'] = pred_date.strftime('%Y-%m-%d')
    output_gdf = geopandas.GeoDataFrame(pred_grid, geometry=geopandas.points_from_xy(pred_grid.lng, pred_grid.lat))
    slimmer_gdf = output_gdf.drop(['lat', 'lng', 'population'], axis=1)

    # return response.get('point_data')
    return json.loads(slimmer_gdf.to_json())
