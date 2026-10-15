"""Regenerate reference data: `python -m rsv_reference [function ...]` (see cli.py)."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
