#!/usr/bin/env bash
set -euo pipefail

_find_config() {
    local dir="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
    for _ in 1 2 3 4 5; do
        [ -f "$dir/config.yaml" ] && echo "$dir/config.yaml" && return
        dir="$(dirname "$dir")"
    done
    echo "ERROR: config.yaml not found" >&2; exit 1
}

CONFIG="$(_find_config)"
REPO_ROOT="$(dirname "$CONFIG")"

_cfg() {
    python3 -c "
import re, sys
section, key, path = sys.argv[1], sys.argv[2], sys.argv[3]
in_sec = False
for line in open(path):
    if re.match(rf'^{section}:', line): in_sec = True
    elif in_sec and re.match(r'^\S', line): in_sec = False
    elif in_sec and re.match(rf'\s+{key}:', line):
        print(line.split(':', 1)[1].strip().strip('\"\''))
        sys.exit()
" "$1" "$2" "$CONFIG"
}

PYTHON_GPU="${PYTHON_GPU:-$(_cfg environments python_gpu)}"
export PATH="$(dirname "$PYTHON_GPU"):$PATH"

cd "$REPO_ROOT/02_pocket_analysis"
"$PYTHON_GPU" af2bind.py
