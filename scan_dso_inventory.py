#!/usr/bin/env python3
"""
scan_dso_inventory.py

Scan an astrophotography archive and summarize:
- Target names
- Filters used
- Exposure totals
- Imaging dates
- File counts by filter
- Per-night breakdowns

Works best with FITS files, but can also fall back to filename parsing.

Examples:
    python scan_dso_inventory.py "E:\\Astro"
    python scan_dso_inventory.py "E:\\Astro" --csv "C:\\Temp\\astro_inventory.csv"
    python scan_dso_inventory.py "E:\\Astro" --include-unknown
"""

from __future__ import annotations

import argparse
import csv
import os
import re
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

STOP_WORDS = {
    "light", "dark", "flat", "bias", "darkflat", "dark-flat", "dark_flat",
    "bin1", "bin2", "bin3", "bin4",
    "gain0", "gain100", "gain101", "gain120", "gain200",
    "mono", "osc", "autorun", "image", "images",
    "temp", "cool", "cooled",
    "session", "seq", "sequence",
    "deg", "adus", "adu",
    "raw", "stack", "processed",
    "master", "integration", "registered", "calibrated",
    "redcat", "ultracat", "wifd",
    "asi2600mm", "asi2600mc", "2600mm", "2600mc",
    "asi533mc", "533mc", "asi533mm", "533mm",
    "asi6200mm", "asi6200mc", "6200mm", "6200mc",
    "rc51", "rc71", "rc91", "uc56", "uc76",
}

FILTER_ALIASES = {
    "L": {"l", "lum", "luminance"},
    "R": {"r", "red"},
    "G": {"g", "green"},
    "B": {"b", "blue"},
    "Ha": {"ha", "hα", "halpha", "h-alpha", "h_alpha"},
    "OIII": {"oiii", "oiii3", "o3", "o-iii", "oxygen", "oxygeniii"},
    "SII": {"sii", "s2", "s-ii", "sulfurii", "sulphurii"},
    "Hb": {"hb", "hbeta", "h-beta", "h_beta"},
    "DualBand": {"lextreme", "l-extreme", "lenhance", "l-enhance", "duoband", "dualband"},
    "RGB": {"rgb"},
}

TARGET_PATTERNS = [
    r"\b(M\s?\d{1,3})\b",
    r"\b(NGC\s?\d{1,4})\b",
    r"\b(IC\s?\d{1,4})\b",
    r"\b(SH2[-\s]?\d{1,3})\b",
    r"\b(B\d{1,3})\b",
    r"\b(CALDWELL\s?\d{1,3})\b",
    r"\b(C\s?\d{1,3})\b",
]

DATE_PATTERNS = [
    r"\b(20\d{2})(\d{2})(\d{2})\b",          # YYYYMMDD
    r"\b(20\d{2})-(\d{2})-(\d{2})\b",        # YYYY-MM-DD
    r"\b(20\d{2})_(\d{2})_(\d{2})\b",        # YYYY_MM_DD
]

EXPOSURE_PATTERNS = [
    r"\b(\d+(?:\.\d+)?)s\b",                 # 300s
    r"\bexp(?:osure)?[_\- ]?(\d+(?:\.\d+)?)\b",
]

HEADER_DATE_KEYS = [
    "DATE-OBS",
    "DATEOBS",
    "DATE_OBS",
]

HEADER_EXPOSURE_KEYS = [
    "EXPTIME",
    "EXPOSURE",
    "EXPOSURE",
]

HEADER_FILTER_KEYS = [
    "FILTER",
    "FILT",
]

HEADER_TARGET_KEYS = [
    "OBJECT",
    "OBJNAME",
    "TARGET",
]

CALIBRATION_WORDS = {
    "dark", "flat", "bias", "darkflat", "dark-flat", "dark_flat",
    "masterdark", "masterflat", "masterbias"
}


@dataclass
class FileInfo:
    path: Path
    target: str
    filter_name: str
    exposure_s: Optional[float]
    image_date: Optional[str]  # YYYY-MM-DD
    file_type: str


@dataclass
class TargetSummary:
    files: list[FileInfo] = field(default_factory=list)

    def total_exposure_s(self) -> float:
        return sum(f.exposure_s or 0.0 for f in self.files)

    def file_count(self) -> int:
        return len(self.files)


def normalize_spaces(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def clean_text(text: str) -> str:
    text = text.replace("_", " ").replace("-", " ")
    text = normalize_spaces(text)
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


def normalize_catalog_name(name: str) -> str:
    name = name.upper().replace("  ", " ")
    name = re.sub(r"\s+", " ", name).strip()

    if name.startswith("CALDWELL"):
        m = re.fullmatch(r"CALDWELL\s?(\d{1,3})", name)
        if m:
            return f"Caldwell {m.group(1)}"
        return name.title()

    for prefix in ("M", "NGC", "IC", "SH2", "B", "C"):
        m = re.fullmatch(rf"{prefix}\s?(\d{{1,4}})", name)
        if m:
            return f"{prefix} {m.group(1)}"

    return name.title()


def looks_like_calibration(path: Path) -> bool:
    parts = re.split(r"[ _\-.]+", str(path).lower())
    return any(part in CALIBRATION_WORDS for part in parts)


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

    value = clean_text(value)

    for pattern in TARGET_PATTERNS:
        match = re.search(pattern, value, flags=re.IGNORECASE)
        if match:
            return normalize_catalog_name(match.group(1))

    return value if value else None


def extract_target_from_filename(filename: str) -> Optional[str]:
    base = clean_text(Path(filename).stem)

    for pattern in TARGET_PATTERNS:
        match = re.search(pattern, base, flags=re.IGNORECASE)
        if match:
            return normalize_catalog_name(match.group(1))

    parts = re.split(r"[ _]+", base)
    filtered = []

    for part in parts:
        p = part.strip().lower()
        if not p:
            continue
        if p in STOP_WORDS:
            continue
        if p in CALIBRATION_WORDS:
            continue
        if re.fullmatch(r"\d+(\.\d+)?s", p):
            continue
        if re.fullmatch(r"\d{8}(-\d{6})?", p):
            continue
        if re.fullmatch(r"\d{1,6}", p):
            continue
        if re.fullmatch(r"\d+deg", p):
            continue
        if re.fullmatch(r"-?\d+(\.\d+)?c", p):
            continue
        filtered.append(part)

    if not filtered:
        return None

    candidate = clean_text(" ".join(filtered[:4]))
    if len(candidate) < 2:
        return None

    return candidate


def identify_target(path: Path) -> str:
    target = read_fits_target(path)
    if target:
        return target

    target = extract_target_from_filename(path.name)
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

    # safer patterns so single-letter filters don't match inside words
    patterns = [
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
    filt = read_fits_filter(path)
    if filt:
        return filt

    filt = extract_filter_from_filename(path.name)
    if filt:
        return filt

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

    # FITS style like 2026-04-19T02:13:44.123
    for fmt in (
        "%Y-%m-%dT%H:%M:%S.%f",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d",
    ):
        try:
            dt = datetime.strptime(value[:len(fmt.replace("%f", "000000"))], fmt)
            return dt.strftime("%Y-%m-%d")
        except Exception:
            pass

    # Regex fallback
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

        print(f"{target}")
        print(f"  Files: {len(files)}")
        print(f"  Total integration: {format_seconds(total_exposure)}")

        # Filters
        by_filter = defaultdict(list)
        for f in files:
            by_filter[f.filter_name].append(f)

        print("  Filters:")
        for filt, filt_files in sorted(by_filter.items(), key=lambda kv: (kv[0] != "Unknown", kv[0])):
            filt_count = len(filt_files)
            filt_exp = sum(x.exposure_s or 0.0 for x in filt_files)
            print(f"    {filt}: {filt_count} file(s), {format_seconds(filt_exp)}")

        # Nights
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

        # Example folders
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

    print_report(grouped)

    if args.csv:
        csv_path = Path(args.csv)
        write_csv(grouped, csv_path)
        print(f"CSV written to: {csv_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())