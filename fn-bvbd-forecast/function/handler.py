import json
import geopandas
import pandas as pd
import numpy as np
import disarm_gears
import uuid
import requests
from datetime import datetime, timedelta


def run_function(params: dict):
    #
    # 1. Handle input
    #

    # Set random seed
    np.random.seed(1000)

    # redirecting STDOUT to avoid over-chatty packages/libraries
    original = sys.stdout
    sys.stdout = open('dummy-stdout-file', 'w')

    end_date = params.get('end_date')
    observed_periods = params.get('observed_periods')
    point_data = params.get('point_data')
    layer_names = params.get('layer_names')

    #
    # 2. Process
    #

    # Load stored data
    tha_shp = geopandas.GeoDataFrame.from_file('help_data/tha_villages/tha_villages.shp')
    #tha_attrib = pd.read_csv('help_data/tha_pred_points.csv')
    tha_attrib = pd.read_csv('help_data/tha_attrib_covs.csv')

    # Add covariates to the attributes table
    #if layer_names is not None:
    if False:
        # Call fn-covariate-extractor
        open_faas_link = 'http://faas.srv.disarm.io/function/fn-covariate-extractor'
        req_options = {
            'points': "https://www.dropbox.com/s/1x57v4ilyk4jbe6/tha_centroids.json?dl=1",
            'layer_names': layer_names
        }
        covs_response = requests.post(open_faas_link, json=req_options)
        covs_response_json = covs_response.json()
        if covs_response_json['type'] == 'error':
            msg = "Problem with remote function call: " + covs_response_json['result']
            raise Exception(msg)
        covs_result = covs_response.json()['result']
        covs_gdf = geopandas.GeoDataFrame.from_features(covs_result['features'])
        covs_data = pd.DataFrame(covs_gdf[[col for col in covs_gdf.columns if col != covs_gdf._geometry_column_name]])

        # Merge output into input_data
        tha_attrib = pd.concat([tha_attrib, covs_data], axis=1)

    # Define timeline and spatiotemporal grid
    tha_frame = disarm_gears.frames.TilePattern(geometries=tha_shp.geometry, attributes=tha_attrib)
    timeline = disarm_gears.frames.Timeframe(start=None, end=end_date, length=observed_periods, by='day', step=28)
    spacetime_grid = tha_frame.make_attributes_series(knots=np.arange(observed_periods), var_name='knot')
    pred_grid = tha_frame.make_attributes_series(knots=np.array([6]), var_name='knot')

    # Make a GeoPandas DataFrame from input data
    input_data = geopandas.GeoDataFrame.from_features(point_data['features'])
    input_xy = np.vstack([input_data.geometry.x, input_data.geometry.y]).T
    input_data.drop(labels=['geometry'], axis=1, inplace=True)
    input_data['knot'] = timeline.which_knot(np.array(input_data.date))
    input_data['polygon_id'] = tha_frame.locate(input_xy)
    input_data = input_data.loc[input_data.knot > -1, :]
    input_data = input_data.loc[input_data.polygon_id > -1, :]

    # Add input data to spatiotemporal grid
    input_data = input_data.groupby(by=['polygon_id', 'knot'], as_index=False).sum()[['polygon_id', 'knot', 'total_cases']]
    spacetime_grid['total_cases'] = 0
    for i, rowi in input_data.iterrows():
        spacetime_grid.loc[np.logical_and(spacetime_grid.index == rowi.polygon_id,
                                          spacetime_grid.knot == rowi.knot), 'total_cases'] += rowi.total_cases

    # Define and fit mgcv model
    rho0 = 7.
    rho1 = 1.5
    gam_formula = f"total_cases ~ offset(log(population)) + te(lng, lat, knot, bs='gp', m=list(c(2, {rho0}, 2), c(2, {rho1}, 2)), d=c(2, 1), k=c(9, 3))"
    if layer_names is not None:
        gam_formula = [gam_formula] + [f'{i}' for i in layer_names]
        gam_formula = '+'.join(gam_formula)

    gam = disarm_gears.r_plugins.r_methods.mgcv_fit(data=spacetime_grid, formula=gam_formula, family='poisson',
                                                    weights=None, method='fREML', bam=True, chunk_size=20000)

    # Make predictions and simulations
    y_pred = disarm_gears.r_plugins.r_methods.mgcv_predict(gam, data=pred_grid, response_type='response')
    y_sims = disarm_gears.r_plugins.r_methods.mgcv_posterior_samples(gam, data=pred_grid, n_samples=n_samples,
                                                                     response_type='link')

    #
    # 3. Package output
    #

    phi = 10000 * y_pred / pred_grid.population
    pred_grid['total_cases'] = y_pred
    pred_grid['incidence'] = phi

    n_samples = 500
    prob_cases = np.sum(y_sims + np.log(pred_grid.population[None, :]) >= 0, axis=0) / n_samples
    risk_class = np.zeros_like(prob_cases)
    for i, pi in enumerate([.1, .25, .5, .75, .975]):
        risk_class[prob_cases >= pi] = i + 1
    pred_grid['prob_cases'] = prob_cases
    pred_grid['risk_class'] = risk_class
    pred_date = timeline.end + timedelta(28)
    pred_grid['date'] = pred_date.strftime('%Y-%m-%d')

    output_gdf = geopandas.GeoDataFrame(pred_grid, geometry=geopandas.points_from_xy(pred_grid.lng, pred_grid.lat))
    slimmer_gdf = output_gdf.drop(['lat', 'lng', 'knot', 'population'], axis=1)

    # Restore STDOUT
    sys.stdout = original

    # return response.get('point_data')
    return json.loads(slimmer_gdf.to_json())
