from pathlib import Path
import argparse

def rename_files(root: Path, filter_tag: str, dry_run: bool = True):
    count = 0

    # Normalize filter tag (remove spaces, uppercase)
    filter_tag = filter_tag.strip()

    for path in root.rglob("*.fit"):
        name = path.name

        # Only modify files missing filter after camera name
        # if "2600MM_" in name and f"2600MM_{filter_tag}_" not in name:
        if "2600MM_" in name and not any(f"2600MM_{f}_" in name for f in ["H", "O", "S", "Ha", "OIII", "SII", "L", "R", "G", "B"]):
            # Avoid double-inserting if already has another filter
            if "2600MM_" in name:
                new_name = name.replace("2600MM_", f"2600MM_{filter_tag}_")

                new_path = path.with_name(new_name)

                if dry_run:
                    print(f"[DRY RUN] {path} -> {new_path}")
                else:
                    path.rename(new_path)
                    print(f"Renamed: {path.name} -> {new_name}")

                count += 1

    print(f"\nProcessed {count} file(s)")


def main():
    parser = argparse.ArgumentParser(description="Insert missing filter tag into FITS filenames")
    parser.add_argument("root", help="Root folder to scan")
    parser.add_argument("--filter", required=True, help="Filter tag to insert (e.g. H, O, S, Ha, OIII)")
    parser.add_argument("--apply", action="store_true", help="Actually rename files (default is dry run)")
    args = parser.parse_args()

    root = Path(args.root)

    if not root.exists():
        print(f"Path does not exist: {root}")
        return

    rename_files(
        root=root,
        filter_tag=args.filter,
        dry_run=not args.apply
    )


if __name__ == "__main__":
    main()