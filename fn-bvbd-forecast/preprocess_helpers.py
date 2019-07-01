import base64
import hashlib
import json
import os
import uuid
from urllib.parse import urlparse
from urllib.request import urlretrieve

# Mutates the params, replacing with local temporary file if required
import config


# Required param must exist
def required_exists(key, params):
    if key not in params:
        raise ValueError(f'Required param \'{key}\' not received.')


# If exists, param must be of given type
def is_type(key, params: dict, param_type):
    param = params.get(key)

    if param is None:
        return

    if not isinstance(params[key], param_type):
        raise ValueError(f'Params \'{key}\' is not of type {param_type}')

def write_temp_from_url_or_base64(key, params):
    value = params[key]

    if not isinstance(value, str):
        # Is an object, guessing is JSON
        write_to_file(key, params)
    elif is_url(value):
        download_to_file(key, params)
    else:
        decode_base64_to_file(key, params)


# Extracts the value from params and writes string to a temp file
# Mutates params
def write_to_file(key, params):
    value = params[key]

    filename = temp_filename()

    # Write STRING to file
    open(filename, 'w').write(json.dumps(value))

    params[key] = filename


# Extracts the URL from params, downloads from URL and writes bytes to a temp file
# Mutates params
def download_to_file(key: str, params: dict):
    url = params[key]

    # Hash the URL to create a temp filename
    hashed_url = hash_url(url)

    filename = os.path.join(config.TEMP, hashed_url)

    # Download from URL to temporary file
    urlretrieve(url, filename)

    # Replace params[key] with temporary file name
    params[key] = filename


# Extracts the value from params, decodes from base64 string and writes to a temp file
# Mutates params
def decode_base64_to_file(key: str, params: dict):
    value = params[key]

    # Decode base64 string
    decoded = base64.b64decode(bytes(value, "utf-8"))

    filename = temp_filename()

    # Write BYTES to file
    open(filename, 'wb').write(decoded)

    # Write to temporary file
    # Replace params[key] with temporary file name
    params[key] = filename


# https://stackoverflow.com/a/52455972/678255
def is_url(url):
    try:
        result = urlparse(url)
        return all([result.scheme, result.netloc])
    except ValueError:
        return False


# Return MD5 hash of URL
def hash_url(url: str) -> str:
    return hashlib.md5(url.encode('utf-8')).hexdigest()


# Create a random filename in the temporary directory
# Creating this using tempfile.NamedTemporaryFile() or similar results in the file
# being closed as soon as the context manager closes: doesn't seem to be another way
# to keep the temporary file around long enough to be used in the runner.
def temp_filename():
    random_string = uuid.uuid4().hex
    return os.path.join(config.TEMP, random_string)
