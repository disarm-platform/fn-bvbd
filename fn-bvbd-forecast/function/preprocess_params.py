from preprocess_helpers import write_temp_from_url_or_base64, required_exists, is_type


def preprocess(params: dict):
    required_exists('point_data', params)

    # Optional, but if exists, must be string
    is_type('uncertainty_type', params, str)