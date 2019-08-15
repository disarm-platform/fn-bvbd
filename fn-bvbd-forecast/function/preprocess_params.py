from preprocess_helpers import write_temp_from_url_or_base64, required_exists, is_type


def preprocess(params: dict):
    required_exists('point_data', params)
    required_exists('observed_periods', params)
    required_exists('end_date', params)
