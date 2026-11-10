    

#!/usr/bin/env python3
"""Download ALICE Grid AO2Ds into D0/data/NNN, resuming safely.

An existing ROOT file is considered *plausible*, not completely verified, if it
has a ROOT header and is >= 1 MiB. Use ROOT to inspect any suspect files.
"""
import argparse
import shutil
import subprocess
import time
from pathlib import Path

GRID = '/alice/sim/2026/LHC26g5/0/571211'
MIN_SIZE = 1024 * 1024


def plausible_root(path):
    try:
        if path.stat().st_size < MIN_SIZE:
            return False
        with path.open('rb') as f:
            return f.read(4) == b'root'
    except (OSError, ValueError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parent / 'data')
    parser.add_argument('--start', type=int, default=23)
    parser.add_argument('--target', type=int, default=200)
    parser.add_argument('--max-id', type=int, default=999)
    parser.add_argument('--retries', type=int, default=3)
    parser.add_argument('--min-free-gb', type=float, default=5)
    args = parser.parse_args()
    data = args.data_dir.expanduser().resolve()
    data.mkdir(parents=True, exist_ok=True)
    if not shutil.which('alien.py'):
        parser.error('alien.py is not on PATH. Activate the ALICE environment.')
    if args.start < 0 or args.target < 1 or args.max_id < args.start or args.retries < 1:
        parser.error('Invalid start, target, max-id, or retries')

    # Only count three-digit directories with plausible AO2D.root files.
    have = {p.parent.name for p in data.glob('[0-9][0-9][0-9]/AO2D.root')
            if plausible_root(p)}
    print(f'Found {len(have)} plausible local AO2Ds. Target: {args.target}.', flush=True)
    print('Note: ROOT header/size check does not prove full tree integrity.', flush=True)
    if len(have) >= args.target:
        return 0
    failed = []
    for number in range(args.start, args.max_id + 1):
        if len(have) >= args.target:
            break
        idx = f'{number:03d}'
        if idx in ('001', '013') or idx in have:
            continue
        folder = data / idx
        final = folder / 'AO2D.root'
        if final.exists() and not plausible_root(final):
            print(f'{idx}: existing file is suspect; leaving untouched and skipping. Inspect it manually.', flush=True)
            failed.append(idx)
            continue
        folder.mkdir(parents=True, exist_ok=True)
        temp = folder / 'AO2D.root.part'
        remote = f'alien://{GRID}/{idx}/AO2D.root'
        success = False
        for attempt in range(1, args.retries + 1):
            available = shutil.disk_usage(data).free / (1024**3)
            if available < args.min_free_gb:
                print(f'STOP: only {available:.1f} GiB free (threshold {args.min_free_gb:.1f} GiB).', flush=True)
                return 2
            temp.unlink(missing_ok=True)
            print(f'{idx}: downloading (attempt {attempt}/{args.retries}); {len(have)}/{args.target} locally', flush=True)
            try:
                cp = subprocess.run(['alien.py', 'cp', remote, 'file://' + str(temp)],
                                    timeout=1800)
                if cp.returncode == 0 and plausible_root(temp):
                    temp.replace(final)
                    have.add(idx)
                    print(f'{idx}: saved {final} ({final.stat().st_size / 1e6:.1f} MB)', flush=True)
                    success = True
                    break
                print(f'{idx}: copy failed or returned invalid/incomplete file (exit {cp.returncode})', flush=True)
            except subprocess.TimeoutExpired:
                print(f'{idx}: copy timed out', flush=True)
            except OSError as exc:
                print(f'{idx}: copy error: {exc}', flush=True)
            finally:
                temp.unlink(missing_ok=True)
            if attempt < args.retries:
                time.sleep(30 * attempt)
        if not success:
            failed.append(idx)
    print(f'Finished: {len(have)}/{args.target} plausible local AO2Ds.', flush=True)
    if failed:
        print('Unsuccessful/suspect IDs (can retry later): ' + ', '.join(failed), flush=True)
    return 0 if len(have) >= args.target else 2


if __name__ == '__main__':
    raise SystemExit(main())

