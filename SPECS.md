# BVBD Incidence predictor

Give us a bunch of GeoJSON points with number of positive cases within the
as well as a GeoJSON of prediction points, and we'll predict the probability of occurrence at each prediction point.

This function depends on a pre-defined framework of village locations and their population

## Parameters

A nested JSON object containing:
- `point_data` - {GeoJSON FeatureCollection} Required. Features with following properties:
  - `date` - {string yyyy-mm-dd} Required. Date of reference of the cases observed
  - `total_cases` - {integer} Required. Number of individuals infected at each location/date (only non-zero cases). By default places/dates with no cases reported are assumed to have zero cases.
  - `id` - {string} Optional id for each point. Must be unique. If not provided, 1:n (where n is the number of Features in the FeatureCollection) will be used.

- `end_date` - {string yyyy-mm-dd} Required. For the point_data provided, the last day of information included in it. This function will make a forecast of the number of cases 28 days after this date.
- `observed_periods` - {integer} Required. Number of periods of 4 weeks included in the point_data, until end_date.
- `layer_names` - {array of strings} Optional. Default is to run with only latitude and longitude. Names relating to the covariate to use to model and predict. See [here](https://github.com/disarm-platform/fn-covariate-extractor/blob/master/SPECS.md) for options.


## Constraints

- maximum number of points/features
- maximum number of layers is XX
- can only include points within a single country

## Response

`point_data` {GeoJSON FeatureCollection} with the following fields: 
- `id` - as defined by user or 1:n (where n is the number of Features in the FeatureCollection)
- `date` - Last day of the 4-weeks period covered by the predictions
- `total_cases` - predicted number of cases
- `incidence` - predicted incidence of malaria
- `prob_cases` - probability of having at least one case
- `risk_class` - categories of risk based on prob_cases as follows
  - 0: prob_cases < .1
  - 1: prob_cases .25
  - 2: prob_cases < .5
  - 3: prob_cases < .75
  - 4: prob_cases < .975
  - 5: prob_cases >= .975
