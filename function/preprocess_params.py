def required_exists(key, params):
    if key not in params:
        raise ValueError(f'Required param \'{key}\' not received.')

def preprocess(params: dict):
    required_exists('point_data', params)
    required_exists('observed_periods', params)
    required_exists('end_date', params)
