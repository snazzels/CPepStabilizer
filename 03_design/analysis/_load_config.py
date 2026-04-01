"""Load config.yaml from the repository root."""

import os
import yaml


def load_config():
    """Walk up from this file's location to find and load config.yaml."""
    here = os.path.dirname(os.path.abspath(__file__))
    path = here
    for _ in range(6):
        candidate = os.path.join(path, "config.yaml")
        if os.path.isfile(candidate):
            with open(candidate) as f:
                return yaml.safe_load(f)
        path = os.path.dirname(path)
    raise FileNotFoundError(f"config.yaml not found above {here}")
