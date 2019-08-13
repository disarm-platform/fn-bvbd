# BVBD Incidence predictor

Give us a set of GeoJSON points with the number of positive cases and we'll
predict the disease incidence and the probability of having one or more new
cases per village.

The function requires an input file with positive cases only. Villages with no
cases reported in the input file will be considered as having zero cases.

When processing the data, the number of cases are aggregated by village in
periods of four weeks. That is how they should be reported in the input file.

This function depends on a pre-defined framework of village locations, their
population and bioclimate data.

Different spatial and spatiotemporal models are fitted to the data. The results
reported back correspond to the model selected as most adequate according to a
GCV criterion. If the model chosen is spatiotemporal, the results reported are
the predictions corresponding to the last period. If the model chosen is
spatial only, then all periods (if more than one) are treated as different
observations of the same spatial process.

## Parameters

A nested JSON object containing:
- `point_data` - {GeoJSON FeatureCollection} Required. Features with following
properties:
  - `date` - {string yyyy-mm-dd} Required. Reference date of the cases observed
  - `total_cases` - {integer} Required. Number of individuals infected at each
location/date (only non-zero cases). By default places/dates with no cases
reported are assumed to have zero cases.

- `end_date` - {string yyyy-mm-dd} Required. For the point_data provided, the
last day of information included in it.
- `observed_periods` - {integer} Required. Number of periods of 4 weeks included
in the point_data, until end_date.
- `layer_names` - {array of strings} Optional. Default is to run with only
latitude and longitude. Names relating to the covariate to use to model and
predict. See [here]
(https://github.com/disarm-platform/fn-covariate-extractor/blob/master/SPECS.md)
for options.


## Constraints

- observed_periods > 1
- maximum number of periods < 10

## Response

`point_data` {GeoJSON FeatureCollection} with the following fields: 
- `id` - as defined by user or 1:n (where n is the number of Features in the
FeatureCollection)
- `date` - Last day of the 4-weeks period covered by the predictions
- `total_cases` - predicted number of cases
- `incidence` - predicted incidence of malaria for each 10000 inhabitants
- `prob_cases` - probability of having at least one case
- `risk_class` - categories of risk based on prob_cases as follows
  - 0: prob_cases < .1
  - 1: prob_cases .25
  - 2: prob_cases < .5
  - 3: prob_cases < .75
  - 4: prob_cases < .975
  - 5: prob_cases >= .975
