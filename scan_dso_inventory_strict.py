#!/usr/bin/env python3
"""
scan_dso_inventory_strict.py

Scan an astrophotography archive and summarize:
- Target names (strict catalog-style matching only)
- Filters used
- Exposure totals
- Imaging dates
- File counts by filter
- Per-night breakdowns

Target detection is intentionally strict:
accepted examples:
    M31, M 31
    NGC2244, NGC 2244
    IC443, IC 443
    Sh2-240, SH2 240
    C49, Caldwell 49
    B33

This avoids false positives from:
    458.3ms
    working folders
    final edits
    processed data names
    random descriptive words

Examples:
    python scan_dso_inventory_strict.py "E:\\Astro"
    python scan_dso_inventory_strict.py "E:\\Astro" --csv "C:\\Temp\\astro_inventory.csv"
    python scan_dso_inventory_strict.py "E:\\Astro" --include-unknown
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from datetime import datetime
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

FILTER_ALIASES = {
    "L": {"l", "lum", "luminance"},
    "R": {"r", "red"},
    "G": {"g", "green"},
    "B": {"b", "blue"},
    "Ha": {"ha", "hα", "halpha", "h-alpha", "h_alpha"},
    "OIII": {"oiii", "o3", "o-iii", "oxygeniii"},
    "SII": {"sii", "s2", "s-ii", "sulfurii", "sulphurii"},
    "Hb": {"hb", "hbeta", "h-beta", "h_beta"},
    "DualBand": {"lextreme", "l-extreme", "lenhance", "l-enhance", "duoband", "dualband"},
    "RGB": {"rgb"},
    "OSC": {"osc", "color","colour"},
}

CALIBRATION_WORDS = {
    "dark", "flat", "bias", "darkflat", "dark-flat", "dark_flat",
    "masterdark", "masterflat", "masterbias",
    "flatdark", "flat-dark", "flat_dark"
}

HEADER_DATE_KEYS = ["DATE-OBS", "DATEOBS", "DATE_OBS"]
HEADER_EXPOSURE_KEYS = ["EXPTIME", "EXPOSURE"]
HEADER_FILTER_KEYS = ["FILTER", "FILT"]
HEADER_TARGET_KEYS = ["OBJECT", "OBJNAME", "TARGET"]
HEADER_BAYER_KEYS = ["BAYERPAT", "BAYERPATTERN", "BAYER_PATTERN"]
HEADER_CAMERA_KEYS = ["INSTRUME", "INSTRUMENT", "CAMERA", "XPIXSZ"]

DATE_PATTERNS = [
    r"\b(20\d{2})(\d{2})(\d{2})\b",   # YYYYMMDD
    r"\b(20\d{2})-(\d{2})-(\d{2})\b", # YYYY-MM-DD
    r"\b(20\d{2})_(\d{2})_(\d{2})\b", # YYYY_MM_DD
]

EXPOSURE_PATTERNS = [
    r"\b(\d+(?:\.\d+)?)s\b",               # 300s
    r"\bexp(?:osure)?[_\- ]?(\d+(?:\.\d+)?)\b",
]

# Strict catalog-style target patterns.
# These are used both on FITS header values and on path/filename strings.
STRICT_TARGET_PATTERNS = [
    ("M",        re.compile(r"(?<![A-Z0-9])M\s*?(\d{1,3})(?!\d)", re.IGNORECASE)),
    ("NGC",      re.compile(r"(?<![A-Z0-9])NGC\s*?(\d{1,4})(?!\d)", re.IGNORECASE)),
    ("IC",       re.compile(r"(?<![A-Z0-9])IC\s*?(\d{1,4})(?!\d)", re.IGNORECASE)),
    ("SH2",      re.compile(r"(?<![A-Z0-9])SH\s*2[-\s]*?(\d{1,3})(?!\d)", re.IGNORECASE)),
    ("SH2",      re.compile(r"(?<![A-Z0-9])SH2[-\s]*?(\d{1,3})(?!\d)", re.IGNORECASE)),
    ("C",        re.compile(r"(?<![A-Z0-9])C(?:ALDWELL)?\s*?(\d{1,3})(?!\d)", re.IGNORECASE)),
    ("B",        re.compile(r"(?<![A-Z0-9])B\s*?(\d{1,3})(?!\d)", re.IGNORECASE)),
]

COMMON_NAME_MAP = {
    # Galaxies
    "M 31": "Andromeda Galaxy",
    "M 32": "Andromeda Galaxy Companion",
    "M 33": "Triangulum Galaxy",
    "M 49": "Virgo Galaxy",
    "M 51": "Whirlpool Galaxy",
    "M 58": "Virgo Galaxy",
    "M 59": "Virgo Galaxy",
    "M 60": "Virgo Galaxy",
    "M 61": "Virgo Galaxy",
    "M 63": "Sunflower Galaxy",
    "M 64": "Black Eye Galaxy",
    "M 65": "Leo Triplet Galaxy",
    "M 66": "Leo Triplet Galaxy",
    "M 74": "Phantom Galaxy",
    "M 77": "Cetus A Galaxy",
    "M 81": "Bode's Galaxy",
    "M 82": "Cigar Galaxy",
    "M 83": "Southern Pinwheel Galaxy",
    "M 84": "Virgo Galaxy",
    "M 85": "Virgo Galaxy",
    "M 86": "Virgo Galaxy",
    "M 87": "Virgo A Galaxy",
    "M 88": "Virgo Galaxy",
    "M 89": "Virgo Galaxy",
    "M 90": "Virgo Galaxy",
    "M 91": "Virgo Galaxy",
    "M 94": "Croc's Eye Galaxy",
    "M 95": "Barred Spiral Galaxy",
    "M 96": "Leo I Galaxy",
    "M 98": "Virgo Galaxy",
    "M 99": "Coma Pinwheel Galaxy",
    "M 100": "Virgo Galaxy",
    "M 101": "Pinwheel Galaxy",
    "M 102": "Spindle Galaxy",
    "M 104": "Sombrero Galaxy",
    "M 106": "Intermediate Spiral Galaxy",
    "M 108": "Surfboard Galaxy",
    "M 109": "Barred Spiral Galaxy",
    "M 110": "Andromeda Galaxy Companion",

    # Nebulae
    "M 1": "Crab Nebula",
    "M 8": "Lagoon Nebula",
    "M 16": "Eagle Nebula",
    "M 17": "Omega Nebula",
    "M 20": "Trifid Nebula",
    "M 27": "Dumbbell Nebula",
    "M 42": "Orion Nebula",
    "M 43": "De Mairan's Nebula",
    "M 57": "Ring Nebula",
    "M 76": "Little Dumbbell Nebula",
    "M 78": "Reflection Nebula",
    "M 97": "Owl Nebula",

    # Star clusters (open & globular)
    "M 2": "Globular Cluster",
    "M 3": "Globular Cluster",
    "M 4": "Globular Cluster",
    "M 5": "Globular Cluster",
    "M 6": "Butterfly Cluster",
    "M 7": "Ptolemy Cluster",
    "M 9": "Globular Cluster",
    "M 10": "Globular Cluster",
    "M 11": "Wild Duck Cluster",
    "M 12": "Globular Cluster",
    "M 13": "Hercules Globular Cluster",
    "M 14": "Globular Cluster",
    "M 15": "Globular Cluster",
    "M 18": "Open Cluster",
    "M 19": "Globular Cluster",
    "M 21": "Open Cluster",
    "M 22": "Sagittarius Globular Cluster",
    "M 23": "Open Cluster",
    "M 24": "Sagittarius Star Cloud",
    "M 25": "Open Cluster",
    "M 26": "Open Cluster",
    "M 28": "Globular Cluster",
    "M 29": "Cooling Tower Cluster",
    "M 30": "Globular Cluster",
    "M 34": "Open Cluster",
    "M 35": "Open Cluster",
    "M 36": "Open Cluster",
    "M 37": "Open Cluster",
    "M 38": "Starfish Cluster",
    "M 39": "Open Cluster",
    "M 40": "Winnecke 4 (Double Star)",
    "M 41": "Open Cluster",
    "M 44": "Beehive Cluster",
    "M 45": "Pleiades",
    "M 46": "Open Cluster",
    "M 47": "Open Cluster",
    "M 48": "Open Cluster",
    "M 50": "Open Cluster",
    "M 52": "Open Cluster",
    "M 53": "Globular Cluster",
    "M 54": "Globular Cluster",
    "M 55": "Globular Cluster",
    "M 56": "Globular Cluster",
    "M 67": "Open Cluster",
    "M 68": "Globular Cluster",
    "M 69": "Globular Cluster",
    "M 70": "Globular Cluster",
    "M 71": "Globular Cluster",
    "M 72": "Globular Cluster",
    "M 73": "Asterism",
    "M 75": "Globular Cluster",
    "M 79": "Globular Cluster",
    "M 80": "Globular Cluster",
    "M 92": "Globular Cluster",
    "M 93": "Open Cluster",
    "M 103": "Open Cluster",
    "M 105": "Elliptical Galaxy",
    "M 107": "Globular Cluster",

    "NGC 2244": "Rosette Nebula",
    "NGC 2237": "Rosette Nebula",
    "IC 443": "Jellyfish Nebula",
    "IC 1805": "Heart Nebula",
    "IC 1848": "Soul Nebula",
    "NGC 1499": "California Nebula",
    "NGC 7000": "North America Nebula",
    "NGC 6960": "Veil Nebula (Western)",
    "NGC 6992": "Veil Nebula (Eastern)",
    "NGC 6888": "Crescent Nebula",
    "NGC 281": "Pacman Nebula",
    "NGC 7635": "Bubble Nebula",

    "B 33": "Horsehead Nebula",
    "SH2-240": "Spaghetti Nebula",
    "SH2-155": "Cave Nebula",

    "C 49": "Rosette Nebula",  # Caldwell mapping
}

# Optional alias normalization. Add to this as needed.
# Left side must already be a normalized catalog target.
TARGET_ALIAS_MAP = {
    # "NGC 2237": "Rosette / NGC 2237",
    # "IC 443": "Jellyfish / IC 443",
    # "IC 1805": "Heart / IC 1805",
    # "IC 1848": "Soul / IC 1848",
}

@dataclass
class FileInfo:
    path: Path
    target: str
    filter_name: str
    exposure_s: Optional[float]
    image_date: Optional[str]
    file_type: str

@dataclass
class TargetSummary:
    files: list[FileInfo] = field(default_factory=list)

    def total_exposure_s(self) -> float:
        return sum(f.exposure_s or 0.0 for f in self.files)

    def file_count(self) -> int:
        return len(self.files)

def clean_text(text: str) -> str:
    text = text.replace("_", " ").replace("-", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text.strip(" .")

def format_seconds(seconds: float) -> str:
    seconds = int(round(seconds))
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    secs = seconds % 60

    if hours > 0:
        return f"{hours}h {minutes}m"
    if minutes > 0:
        return f"{minutes}m {secs}s"
    return f"{secs}s"

def normalize_target(prefix: str, number: str) -> str:
    number = str(int(number))  # strips leading zeroes safely
    prefix = prefix.upper()

    if prefix == "SH2":
        return f"Sh2-{number}"
    if prefix == "C":
        return f"C {number}"
    return f"{prefix} {number}"

def apply_target_alias(target: str) -> str:
    return TARGET_ALIAS_MAP.get(target, target)

def extract_strict_target(text: str) -> Optional[str]:
    if not text:
        return None

    # Normalize separators for easier matching, but keep original enough
    # that compact forms like M31 / NGC2244 still match.
    candidate = str(text).strip()

    for prefix, pattern in STRICT_TARGET_PATTERNS:
        match = pattern.search(candidate)
        if match:
            return apply_target_alias(normalize_target(prefix, match.group(1)))

    return None

def looks_like_calibration(path: Path) -> bool:
    tokens = re.split(r"[ _\-.\\/]+", str(path).lower())
    return any(tok in CALIBRATION_WORDS for tok in tokens if tok)

def read_fits_header_value(path: Path, keys: list[str]) -> Optional[str]:
    if not ASTROPY_AVAILABLE:
        return None
    if path.suffix.lower() not in {".fit", ".fits", ".fts"}:
        return None

    try:
        with fits.open(path, ignore_missing_simple=True) as hdul:
            header = hdul[0].header
            for key in keys:
                value = header.get(key)
                if value is not None:
                    value = str(value).strip()
                    if value:
                        return value
    except Exception:
        return None

    return None

def read_fits_target(path: Path) -> Optional[str]:
    value = read_fits_header_value(path, HEADER_TARGET_KEYS)
    if not value:
        return None
    return extract_strict_target(value)

def extract_target_from_path(path: Path) -> Optional[str]:
    # Search progressively broader scopes:
    # 1. filename stem
    # 2. parent folder names
    # 3. full path
    #
    # This allows things like:
    #   M31AndromedaFinal.fit   -> M 31
    #   E:\Astro\M31AndromedaFinal\Lights\... -> M 31
    # while still ignoring unrelated junk that has no valid catalog token.

    candidates = [
        path.stem,
        path.name,
        path.parent.name,
        str(path),
    ]

    for item in candidates:
        target = extract_strict_target(item)
        if target:
            return target

    return None

def identify_target(path: Path) -> str:
    target = read_fits_target(path)
    if target:
        return target

    target = extract_target_from_path(path)
    if target:
        return target

    return "Unknown"

def normalize_filter_name(value: str) -> str:
    raw = clean_text(value).lower().replace(" ", "")
    for canonical, aliases in FILTER_ALIASES.items():
        if raw in aliases:
            return canonical

    if raw in {"ha7nm", "ha5nm", "ha3nm"}:
        return "Ha"
    if raw in {"oiii7nm", "oiii5nm", "oiii3nm"}:
        return "OIII"
    if raw in {"sii7nm", "sii5nm", "sii3nm"}:
        return "SII"

    return clean_text(value)

def read_fits_filter(path: Path) -> Optional[str]:
    value = read_fits_header_value(path, HEADER_FILTER_KEYS)
    if not value:
        return None
    return normalize_filter_name(value)

def extract_filter_from_filename(filename: str) -> Optional[str]:
    base = Path(filename).stem.lower()
    tokens = re.split(r"[ _\-.]+", base)

    for token in tokens:
        token = token.strip().lower()
        if not token:
            continue
        for canonical, aliases in FILTER_ALIASES.items():
            if token in aliases:
                return canonical

    patterns = [
        (r"l[-_ ]?extreme", "DualBand"),
        (r"l[-_ ]?enhance", "DualBand"),
        (r"dual[-_ ]?band", "DualBand"),
        (r"duoband", "DualBand"),    
        (r"(?:^|[_\-. ])ha(?:$|[_\-. ])", "Ha"),
        (r"(?:^|[_\-. ])oiii(?:$|[_\-. ])", "OIII"),
        (r"(?:^|[_\-. ])o3(?:$|[_\-. ])", "OIII"),
        (r"(?:^|[_\-. ])sii(?:$|[_\-. ])", "SII"),
        (r"(?:^|[_\-. ])s2(?:$|[_\-. ])", "SII"),
        (r"(?:^|[_\-. ])hb(?:$|[_\-. ])", "Hb"),
        (r"(?:^|[_\-. ])lum(?:$|[_\-. ])", "L"),
        (r"(?:^|[_\-. ])l(?:$|[_\-. ])", "L"),
        (r"(?:^|[_\-. ])r(?:$|[_\-. ])", "R"),
        (r"(?:^|[_\-. ])g(?:$|[_\-. ])", "G"),
        (r"(?:^|[_\-. ])b(?:$|[_\-. ])", "B"),
    ]

    for pattern, result in patterns:
        if re.search(pattern, base, flags=re.IGNORECASE):
            return result

    return None

def identify_filter(path: Path) -> str:
    # 1. Explixir FITS/header filter
    filt = read_fits_filter(path)
    if filt:
        return filt

    filt = extract_filter_from_filename(path.name)
    if filt:
        # If the named filter is a dual-band style filter, normalize it
        raw = filt.lower().replace(" ", "")
        if raw in {"lextreme", "l-extreme", "lenhance", "l-enhance", "duoband", "dualband"}:
            return "DualBand"
        return filt
        
    # 2. Filename-based filter
    filt = extract_filter_from_filename(path.name)
    if filt:
        return filt

    # 3. If it looks like OSC, classify as OSC
    if is_osc_file(path):
        return "OSC"

    return "Unknown"

def read_fits_exposure(path: Path) -> Optional[float]:
    value = read_fits_header_value(path, HEADER_EXPOSURE_KEYS)
    if not value:
        return None

    try:
        return float(value)
    except Exception:
        return None

def extract_exposure_from_filename(filename: str) -> Optional[float]:
    base = Path(filename).stem.lower()

    # Reject millisecond-like naming such as 458.3ms; we do not want that
    # interpreted as a target, and we also do not want it treated as light exposure
    # here unless your naming convention explicitly uses ms.
    if re.search(r"\b\d+(?:\.\d+)?ms\b", base):
        return None

    for pattern in EXPOSURE_PATTERNS:
        match = re.search(pattern, base, flags=re.IGNORECASE)
        if match:
            try:
                return float(match.group(1))
            except Exception:
                pass

    return None

def identify_exposure(path: Path) -> Optional[float]:
    exp = read_fits_exposure(path)
    if exp is not None:
        return exp

    exp = extract_exposure_from_filename(path.name)
    if exp is not None:
        return exp

    return None

def parse_date_to_yyyy_mm_dd(value: str) -> Optional[str]:
    value = value.strip()

    # Try ISO-ish forms first
    for splitter in ("T", " "):
        if splitter in value:
            value = value.split(splitter, 1)[0]
            break

    for fmt in ("%Y-%m-%d", "%Y_%m_%d", "%Y%m%d"):
        try:
            dt = datetime.strptime(value, fmt)
            return dt.strftime("%Y-%m-%d")
        except Exception:
            pass

    for pattern in DATE_PATTERNS:
        match = re.search(pattern, value)
        if match:
            y, m, d = match.groups()
            try:
                dt = datetime(int(y), int(m), int(d))
                return dt.strftime("%Y-%m-%d")
            except Exception:
                return None

    return None

def read_fits_date(path: Path) -> Optional[str]:
    value = read_fits_header_value(path, HEADER_DATE_KEYS)
    if not value:
        return None
    return parse_date_to_yyyy_mm_dd(value)

def extract_date_from_filename_or_path(path: Path) -> Optional[str]:
    search_text = str(path)

    for pattern in DATE_PATTERNS:
        match = re.search(pattern, search_text)
        if match:
            y, m, d = match.groups()
            try:
                dt = datetime(int(y), int(m), int(d))
                return dt.strftime("%Y-%m-%d")
            except Exception:
                pass

    return None

def identify_date(path: Path) -> Optional[str]:
    d = read_fits_date(path)
    if d:
        return d

    d = extract_date_from_filename_or_path(path)
    if d:
        return d

    return None

def classify_file_type(path: Path) -> str:
    if looks_like_calibration(path):
        return "Calibration"
    return "Light"

def scan_folder(root: Path, include_unknown: bool = False, lights_only: bool = True) -> dict[str, TargetSummary]:
    grouped: dict[str, TargetSummary] = defaultdict(TargetSummary)

    for dirpath, _, filenames in os.walk(root):
        for filename in filenames:
            path = Path(dirpath) / filename

            if path.suffix.lower() not in SUPPORTED_EXTENSIONS:
                continue

            file_type = classify_file_type(path)
            if lights_only and file_type != "Light":
                continue

            target = identify_target(path)
            if target == "Unknown" and not include_unknown:
                continue

            info = FileInfo(
                path=path,
                target=target,
                filter_name=identify_filter(path),
                exposure_s=identify_exposure(path),
                image_date=identify_date(path),
                file_type=file_type,
            )

            grouped[target].files.append(info)

    return grouped

def format_target_display(target: str) -> str:
    """Helper function to map a target name to a common name.
    
    Args:
        target (str): The standard notation for a target for eample "M 51".
    Returns:
        str: The common name for the target or the standard name.
    """
    common = COMMON_NAME_MAP.get(target)
    if common:
        return f"{target} ({common})"
    return target

def print_report(grouped: dict[str, TargetSummary]) -> None:
    if not grouped:
        print("No matching image files found.")
        return

    print("\nTargets found:\n")

    sorted_targets = sorted(
        grouped.items(),
        key=lambda kv: (-kv[1].file_count(), kv[0].lower())
    )

    for target, summary in sorted_targets:
        files = summary.files
        total_exposure = summary.total_exposure_s()

        # print(f"{target}")
        display_target = format_target_display(target)
        print(f"{display_target}")
        print(f"  Files: {len(files)}")
        print(f"  Total integration: {format_seconds(total_exposure)}")

        by_filter = defaultdict(list)
        for f in files:
            by_filter[f.filter_name].append(f)

        print("  Filters:")
        for filt, filt_files in sorted(by_filter.items(), key=lambda kv: (kv[0] != "Unknown", kv[0])):
            filt_count = len(filt_files)
            filt_exp = sum(x.exposure_s or 0.0 for x in filt_files)
            print(f"    {filt}: {filt_count} file(s), {format_seconds(filt_exp)}")

        by_date = defaultdict(list)
        for f in files:
            by_date[f.image_date or "UnknownDate"].append(f)

        print("  Nights:")
        for date_key, date_files in sorted(by_date.items(), key=lambda kv: kv[0]):
            night_exp = sum(x.exposure_s or 0.0 for x in date_files)
            filter_counts = Counter(x.filter_name for x in date_files)
            filter_bits = ", ".join(
                f"{flt}={count}"
                for flt, count in sorted(filter_counts.items(), key=lambda kv: kv[0])
            )
            print(f"    {date_key}: {len(date_files)} file(s), {format_seconds(night_exp)} [{filter_bits}]")

        folders = sorted({str(f.path.parent) for f in files})
        print("  Example folders:")
        for folder in folders[:5]:
            print(f"    {folder}")
        if len(folders) > 5:
            print(f"    ... plus {len(folders) - 5} more")
        print()

def write_csv(grouped: dict[str, TargetSummary], csv_path: Path) -> None:
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Target",
            "TotalFiles",
            "TotalExposureSeconds",
            "TotalExposureFormatted",
            "Filter",
            "FilterFiles",
            "FilterExposureSeconds",
            "FilterExposureFormatted",
            "Date",
            "DateFiles",
            "DateExposureSeconds",
            "DateExposureFormatted",
        ])

        sorted_targets = sorted(
            grouped.items(),
            key=lambda kv: (-kv[1].file_count(), kv[0].lower())
        )

        for target, summary in sorted_targets:
            files = summary.files
            total_exp = summary.total_exposure_s()

            by_filter = defaultdict(list)
            by_date = defaultdict(list)

            for fi in files:
                by_filter[fi.filter_name].append(fi)
                by_date[fi.image_date or "UnknownDate"].append(fi)

            all_filters = sorted(by_filter.keys())
            all_dates = sorted(by_date.keys())
            max_rows = max(len(all_filters), len(all_dates), 1)

            filter_rows = []
            for filt in all_filters:
                filt_files = by_filter[filt]
                filt_exp = sum(x.exposure_s or 0.0 for x in filt_files)
                filter_rows.append((
                    filt,
                    len(filt_files),
                    filt_exp,
                    format_seconds(filt_exp),
                ))

            date_rows = []
            for date_key in all_dates:
                date_files = by_date[date_key]
                date_exp = sum(x.exposure_s or 0.0 for x in date_files)
                date_rows.append((
                    date_key,
                    len(date_files),
                    date_exp,
                    format_seconds(date_exp),
                ))

            for i in range(max_rows):
                fr = filter_rows[i] if i < len(filter_rows) else ("", "", "", "")
                dr = date_rows[i] if i < len(date_rows) else ("", "", "", "")

                writer.writerow([
                    target if i == 0 else "",
                    len(files) if i == 0 else "",
                    round(total_exp, 2) if i == 0 else "",
                    format_seconds(total_exp) if i == 0 else "",
                    fr[0], fr[1], round(fr[2], 2) if fr[2] != "" else "", fr[3],
                    dr[0], dr[1], round(dr[2], 2) if dr[2] != "" else "", dr[3],
                ])

def normalize_user_target_input(value: str) -> Optional[str]:
    """
    Convert user input like:
        M31, m 31, ngc2244, IC 443
    into normalized form used internally:
        M 31, NGC 2244, IC 443
    """
    if not value:
        return None

    return extract_strict_target(value)
    
def normalize_user_filter_input(value: str) -> Optional[str]:
    """
    Convert user input like:
        ha, Ha, h-alpha, oiii, O3, lum, red
    into normalized internal form:
        Ha, OIII, SII, L, R, G, B, Hb, DualBand, RGB
    """
    if not value:
        return None

    raw = clean_text(value).lower().replace(" ", "").replace("-", "")

    # 1. Alias-based matching (primary)
    for canonical, aliases in FILTER_ALIASES.items():
        if raw in aliases:
            return canonical

    # 2. Direct mapping fallback (CLI-friendly)
    canonical_map = {
        "ha": "Ha",
        "halpha": "Ha",
        "oiii": "OIII",
        "o3": "OIII",
        "sii": "SII",
        "s2": "SII",
        "hb": "Hb",
        "hbeta": "Hb",
        "l": "L",
        "lum": "L",
        "luminance": "L",
        "r": "R",
        "red": "R",
        "g": "G",
        "green": "G",
        "b": "B",
        "blue": "B",
        "rgb": "RGB",
        "dualband": "DualBand",
        "duoband": "DualBand",
        "osc": "OSC",
    }
    
    return canonical_map.get(raw)

def filter_grouped_by_filter(
    grouped: dict[str, TargetSummary],
    only_filter: str
) -> dict[str, TargetSummary]:
    filtered_grouped: dict[str, TargetSummary] = {}

    for target, summary in grouped.items():
        matching_files = [f for f in summary.files if f.filter_name == only_filter]
        if matching_files:
            filtered_grouped[target] = TargetSummary(files=matching_files)

    return filtered_grouped


def is_osc_file(path: Path) -> bool:
    # Strongest signal: FITS Bayer pattern
    bayer = read_fits_header_value(path, HEADER_BAYER_KEYS)
    if bayer:
        return True

    # Filename/path hints
    search_text = str(path).lower()

    osc_hints = [
        "osc",
        "one shot color",
        "oneshotcolor",
        "color",
        "2600mc",
        "2600mcpro",
        "2600mcair",
        "2600mc duo",
        "2600mcduo",
        "533mc",
        "533mcpro",
        "294mc",
        "071mc",
        "585mc",
        "662mc",
    ]

    return any(hint in search_text for hint in osc_hints)

def safe_folder_name(name: str) -> str:
    """
    Make a Windows-safe folder name.
    """
    name = name.replace(":", "_")
    name = name.replace("\\", "_")
    name = name.replace("/", "_")
    name = name.replace("*", "_")
    name = name.replace("?", "_")
    name = name.replace('"', "_")
    name = name.replace("<", "_")
    name = name.replace(">", "_")
    name = name.replace("|", "_")
    return name.strip()


def copy_grouped_files_to_cache(
    grouped: dict[str, TargetSummary],
    cache_root: Path,
    dry_run: bool = False,
) -> None:
    """
    Copy matched files into a local cache organized as:

        cache_root/
          Target/
            Date/
              Filter/
                filename.fit

    Example:
        D:/AstroCache/
          NGC 2239/
            2026-02-12/
              H/
              O/
              S/
    """

    total_files = sum(len(summary.files) for summary in grouped.values())
    copied = 0
    skipped = 0

    print()
    print(f"Preparing to copy {total_files} file(s) to: {cache_root}")
    print()

    cache_root.mkdir(parents=True, exist_ok=True)

    for target, summary in grouped.items():
        target_folder = safe_folder_name(target)

        for info in summary.files:
            date_folder = safe_folder_name(info.image_date or "UnknownDate")
            filter_folder = safe_folder_name(info.filter_name or "UnknownFilter")

            dest_dir = cache_root / target_folder / date_folder / filter_folder
            dest_path = dest_dir / info.path.name

            if dry_run:
                print(f"[DRY RUN] {info.path} -> {dest_path}")
                continue

            dest_dir.mkdir(parents=True, exist_ok=True)

            # Skip if same size already exists
            if dest_path.exists():
                try:
                    if dest_path.stat().st_size == info.path.stat().st_size:
                        skipped += 1
                        continue
                except Exception:
                    pass

            try:
                shutil.copy2(info.path, dest_path)
                copied += 1
            except Exception as exc:
                print(f"ERROR copying: {info.path}")
                print(f"  {exc}")

    if dry_run:
        print()
        print("Dry run complete. No files were copied.")
    else:
        print()
        print(f"Copy complete. Copied: {copied}, skipped existing: {skipped}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Scan an astro archive and summarize targets, filters, exposure, and dates."
    )
    parser.add_argument(
        "root",
        help=r"Root folder/drive to scan, e.g. E:\Astro"
    )
    parser.add_argument(
        "--csv",
        help="Optional CSV output path"
    )
    parser.add_argument(
        "--include-unknown",
        action="store_true",
        help="Include files where target could not be identified"
    )
    parser.add_argument(
        "--include-calibration",
        action="store_true",
        help="Include darks/flats/bias frames too"
    )
    parser.add_argument(
        "--only-target",
        help="Show only a specific target (e.g. M31, NGC2244, IC443)"
    )
    parser.add_argument(
        "--only-filter",
        help="Show only a specific filter (e.g. Ha, OIII, SII, L, R, G, B)"
    )
    parser.add_argument(
    "--copy-to",
    help=r"Copy matched files to a local cache folder, e.g. D:\AstroCache"
    )
    parser.add_argument(
    "--dry-run-copy",
    action="store_true",
    help="Show what would be copied without actually copying files"
)

    args = parser.parse_args()

    root = Path(args.root)
    if not root.exists():
        print(f"Path does not exist: {root}")
        return 1

    print(f"Scanning: {root}")

    if not ASTROPY_AVAILABLE:
        print("Note: astropy is not installed, so FITS headers will not be read.")
        print("Install with: py -m pip install astropy")

    grouped = scan_folder(
        root=root,
        include_unknown=args.include_unknown,
        lights_only=not args.include_calibration,
    )

    # Apply single-target filter if requested
    if args.only_target:
        normalized_target = normalize_user_target_input(args.only_target)

        if not normalized_target:
            print(f"Could not recognize target format: {args.only_target}")
            return 1

        if normalized_target not in grouped:
            print(f"Target '{normalized_target}' not found in dataset.")
            return 0

        grouped = {normalized_target: grouped[normalized_target]}

    # Apply single-filter filter if requested
    if args.only_filter:
        normalized_filter = normalize_user_filter_input(args.only_filter)

        if not normalized_filter:
            print(f"Could not recognize filter: {args.only_filter}")
            return 1

        grouped = filter_grouped_by_filter(grouped, normalized_filter)

        if not grouped:
            if args.only_target:
                print(f"No files found for target '{normalized_target}' with filter '{normalized_filter}'.")
            else:
                print(f"No files found with filter '{normalized_filter}'.")
            return 0

    print_report(grouped)

    if args.csv:
        csv_path = Path(args.csv)
        write_csv(grouped, csv_path)
        print(f"CSV written to: {csv_path}")

    if args.copy_to:
        cache_root = Path(args.copy_to)
        copy_grouped_files_to_cache(
            grouped=grouped,
            cache_root=cache_root,
            dry_run=args.dry_run_copy,
        )

    return 0

if __name__ == "__main__":
    raise SystemExit(main())