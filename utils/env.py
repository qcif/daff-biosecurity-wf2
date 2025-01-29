"""Environment utils."""

import os


def getenv(name):
    """Get environment variable by name."""
    value = os.getenv(name)
    if value is None:
        raise EnvironmentError(f"Environment variable {name} not set")
    return value
