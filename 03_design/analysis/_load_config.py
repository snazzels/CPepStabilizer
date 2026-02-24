"""
Minimal config.yaml reader — no pyyaml dependency.
Handles the two-level structure used in this project's config.yaml.
"""
import re
from pathlib import Path


def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config.yaml"
    result = {}
    section = None
    with open(config_path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.lstrip().startswith('#'):
                continue
            # Top-level section key (no leading spaces, ends with colon)
            if re.match(r'^[a-z_]+:$', line):
                section = line[:-1]
                result[section] = {}
            # Nested key: value (2-space indent)
            elif line.startswith('  ') and ':' in line and section is not None:
                key, _, raw = line.strip().partition(': ')
                if not raw:
                    continue
                raw = raw.strip().strip('"\'')
                if raw == 'true':
                    val = True
                elif raw == 'false':
                    val = False
                else:
                    try:
                        val = int(raw)
                    except ValueError:
                        try:
                            val = float(raw)
                        except ValueError:
                            val = raw
                result[section][key] = val
    return result
