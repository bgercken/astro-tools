#!/usr/bin/env python3
"""
scan_dso_targets.py

Scan a drive/folder for astrophotography files and identify likely DSO targets.

Usage:
    python scan_dso_targets.py "E:\\AstroData"

Optional:
    python scan_dso_targets.py "E:\\AstroData" --csv report.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional

# Optional FITS support
try:
    from astropy.io import fits
    ASTROPY_AVAILABLE = True
except Exception:
    ASTROPY_AVAILABLE = False

SUPPORTED_EXTENSIONS = {
    ".fit", ".fits", ".fts",
    ".xisf",
    ".cr2", ".cr3", ".nef", ".arw", ".dng", ".tif", ".tiff"
}

# Words that are common in astro filenames but are not target names
STOP_WORDS = {
    "light", "dark", "flat", "bias", "darkflat", "dark-flat",
    "bin1", "bin2", "bin3", "bin4",
    "gain0", "gain100", "gain101", "gain120", "gain200",
    "ha", "h", "oiii", "o3", "sii", "s2", "l", "lum", "r", "g", "b",
    "mono", "osc",
    "rc51", "rc71", "rc91", "redcat51", "redcat71", "redcat91",
    "asi2600mm", "asi2600mc", "2600mm", "2600mc",
    "deg", "c", "f",
    "autorun"
}

TARGET_PATTERNS = [
    r"\b(M\s?\d{1,3})\b",
    r"\b(NGC\s?\d{1,4})\b",
    r"\b(IC\s?\d{1,4})\b",
    r"\b(SH2[-\s]?\d{1,3})\b",
    r"\b(B\d{1,3})\b",
    r"\b(Caldwell\s?\d{1,3})\b",
    r"\b(C\s?\d{1,3})\b",
]

def normalize_spaces(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()

def clean_target(text: str) -> str:
    text = text.replace("_", " ").replace("-", " ")
    text = normalize_spaces(text)
    return text.strip(" .")

def read_fits_object_header(path: Path) -> Optional[str]:
    if not ASTROPY_AVAILABLE:
        return None

    if path.suffix.lower() not in {".fit", ".fits", ".fts"}:
        return None

    try:
        with fits.open(path, ignore_missing_simple=True) as hdul:
            header = hdul[0].header
            for key in ("OBJECT", "OBJNAME", "TARGET"):
                value = header.get(key)
                if value:
                    value = str(value).strip()
                    if value:
                        return clean_target(value)
    except Exception:
        return None

    return None

def extract_target_from_filename(filename: str) -> Optional[str]:
    base = Path(filename).stem
    base = clean_target(base)

    # First try well-known catalog patterns
    for pattern in TARGET_PATTERNS:
        match = re.search(pattern, base, flags=re.IGNORECASE)
        if match:
            return normalize_catalog_name(match.group(1))

    # Then try a heuristic:
    # split filename into chunks and keep likely descriptive words
    parts = re.split(r"[ _]+", base)
    filtered = []

    for part in parts:
        p = part.strip().lower()
        if not p:
            continue
        if p in STOP_WORDS:
            continue
        if re.fullmatch(r"\d+(\.\d+)?s", p):   # exposure like 300s
            continue
        if re.fullmatch(r"\d{8}(-\d{6})?", p): # date stamps
            continue
        if re.fullmatch(r"\d{1,4}", p):        # plain numbers
            continue
        if re.fullmatch(r"\d+deg", p):
            continue
        if re.fullmatch(r"-?\d+(\.\d+)?c", p):
            continue
        filtered.append(part)

    if not filtered:
        return None

    # Keep first few descriptive tokens
    candidate = " ".join(filtered[:4])
    candidate = clean_target(candidate)

    if len(candidate) < 2:
        return None

    return candidate

def normalize_catalog_name(name: str) -> str:
    name = name.upper().replace("  ", " ")
    name = re.sub(r"\s+", " ", name).strip()

    # Normalize forms like "M42" -> "M 42"
    for prefix in ("M", "NGC", "IC", "SH2", "B", "C"):
        m = re.fullmatch(rf"{prefix}\s?(\d{{1,4}})", name)
        if m:
            return f"{prefix} {m.group(1)}"
    return name

def identify_target(path: Path) -> str:
    # 1. FITS header
    header_target = read_fits_object_header(path)
    if header_target:
        return header_target

    # 2. Filename heuristic
    filename_target = extract_target_from_filename(path.name)
    if filename_target:
        return filename_target

    return "Unknown"

def scan_folder(root: Path, include_unknown=False):
    grouped = defaultdict(list)

    for dirpath, _, filenames in os.walk(root):
        for filename in filenames:
            path = Path(dirpath) / filename
            if path.suffix.lower() not in SUPPORTED_EXTENSIONS:
                continue

            target = identify_target(path)

            if target == "Unknown" and not include_unknown:
                continue

            grouped[target].append(path)

    return grouped

def print_report(grouped):
    if not grouped:
        print("No supported image files found.")
        return

    print("\nTargets found:\n")
    for target in sorted(grouped, key=lambda k: (-len(grouped[k]), k.lower())):
        files = grouped[target]
        folders = sorted({str(p.parent) for p in files})
        print(f"{target}: {len(files)} file(s)")
        for folder in folders[:5]:
            print(f"  {folder}")
        if len(folders) > 5:
            print(f"  ... plus {len(folders) - 5} more folder(s)")
        print()

def write_csv(grouped, csv_path: Path):
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Target", "FileCount", "ExampleFolder"])
        for target in sorted(grouped, key=lambda k: (-len(grouped[k]), k.lower())):
            files = grouped[target]
            example_folder = str(files[0].parent) if files else ""
            writer.writerow([target, len(files), example_folder])

def main():
    parser = argparse.ArgumentParser(description="Scan astro files and identify DSO targets.")
    parser.add_argument("root", help="Root folder/drive to scan, e.g. E:\\AstroData")
    parser.add_argument("--csv", help="Optional CSV report output path")
    parser.add_argument("--include_unknown", action="store_true", help="Include files where target could not be identified" )
    args = parser.parse_args()

    root = Path(args.root)
    if not root.exists():
        print(f"Path does not exist: {root}")
        return 1

    print(f"Scanning: {root}")
    if not ASTROPY_AVAILABLE:
        print("Note: astropy is not installed, so FITS headers will not be read.")
        print("      Install with: pip install astropy")

    grouped = scan_folder(root, include_unknown=args.include_unknown)
    print_report(grouped)

    if args.csv:
        csv_path = Path(args.csv)
        write_csv(grouped, csv_path)
        print(f"CSV written to: {csv_path}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())