"""Folder layout shared by the scripts (paths are relative to the repository root)."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"

for _d in (PROCESSED, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)
