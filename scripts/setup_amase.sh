#!/bin/bash
# Install Astro AMASE for con-line-search, inside .local-tools/ (git-ignored):
#
#   .local-tools/astro_amase/    the package source (github.com/zfried/astro_amase)
#   .local-tools/amase-venv/     its own Python 3.11 environment
#   .local-tools/amase-data/     its database (CDMS/JPL catalogues, molsim objects)
#
# Astro AMASE needs Python 3.11: its astrochem_embedding dependency requires
# rdkit-pypi, which has no Python 3.12 builds. uv provides Python 3.11 if the
# computer has none. Re-running skips finished steps. Needs internet; the
# downloads are large (PyTorch, RDKit and the AMASE database).
set -euo pipefail
export GIT_PAGER=cat PAGER=cat   # never stop in a pager (git log would open less)
ROOT="$(cd -- "$(dirname -- "$0")/.." && pwd)"
TOOLS="$ROOT/.local-tools"
SRC="$TOOLS/astro_amase"
VENV="$TOOLS/amase-venv"
DATA="$TOOLS/amase-data"
DB_URL="${AMASE_DB_URL:-https://www.dropbox.com/scl/fo/r5trqr0ua8y4yxbcwc47k/AEXQaGWZ88rjBhZpx0wadP0?rlkey=11uzyacot910qttpqhv3g6byr&dl=1}"
mkdir -p "$TOOLS"
step() { printf '\n== %s\n' "$*"; }

step "1/5 uv (Python environment manager)"
if command -v uv >/dev/null 2>&1; then UV="$(command -v uv)"
elif [[ -x "$TOOLS/uv-env/bin/uv" ]]; then UV="$TOOLS/uv-env/bin/uv"
else
  base=""
  for p in "$ROOT/.venv/bin/python" python3.12 python3.13 python3.11 python3; do
    if command -v "$p" >/dev/null 2>&1; then base="$(command -v "$p")"; break; fi
  done
  [[ -n "$base" ]] || { echo "No Python found to install uv."; exit 1; }
  "$base" -m venv "$TOOLS/uv-env"
  "$TOOLS/uv-env/bin/python" -m pip install -q --upgrade pip uv
  UV="$TOOLS/uv-env/bin/uv"
fi
echo "Using $UV"

step "2/5 Astro AMASE source"
if [[ -d "$SRC/.git" ]]; then git -C "$SRC" pull -q --ff-only || echo "(could not update; using the existing copy)"
else git clone -q --depth 1 https://github.com/zfried/astro_amase "$SRC"
fi
git --no-pager -C "$SRC" log -1 --format='commit %h (%ad)' --date=short

step "3/5 Python 3.11 environment and packages (several minutes the first time)"
# uv's own Python builds include Tk; group_selfies (an AMASE dependency) imports
# tkinter, which some Python installs (e.g. Homebrew's) lack.
if [[ -x "$VENV/bin/python" ]] && ! "$VENV/bin/python" -c 'import tkinter' 2>/dev/null; then
  echo "The existing environment has no Tk support; rebuilding it."
  rm -rf "$VENV"
fi
[[ -x "$VENV/bin/python" ]] || UV_PYTHON_PREFERENCE=only-managed "$UV" venv -q --python 3.11 "$VENV"
"$UV" pip install --python "$VENV/bin/python" -e "$SRC"

step "4/5 Astro AMASE database"
find_db() { find "$DATA" -name all_cdms_final_official.csv -print -quit 2>/dev/null | xargs -I{} dirname {}; }
ZIP="$DATA/amase-database.zip"
if [[ -z "$(find_db)" && ! -f "$ZIP" ]]; then
  mkdir -p "$DATA"
  echo "Downloading the database (about 5 GB) from the link in the Astro AMASE README…"
  curl -L --fail --progress-bar -o "$ZIP.part" "$DB_URL"
  mv "$ZIP.part" "$ZIP"
fi
if [[ -f "$ZIP" ]]; then
  # Unpack with Python: Dropbox folder zips contain a "/" entry that makes macOS
  # unzip stop with an error. Files already unpacked with the right size are kept;
  # the zip is deleted once every file is in place.
  echo "Unpacking (checking every file)…"
  "$VENV/bin/python" - "$ZIP" "$DATA" <<'PYCODE'
import os, sys, zipfile
archive, dest = sys.argv[1], os.path.realpath(sys.argv[2])
with zipfile.ZipFile(archive) as z:
    for info in z.infolist():
        name = info.filename.lstrip('/')
        target = os.path.realpath(os.path.join(dest, name))
        if not name or not target.startswith(dest + os.sep):
            continue                      # the "/" entry, or anything outside the folder
        if info.is_dir():
            os.makedirs(target, exist_ok=True); continue
        if os.path.exists(target) and os.path.getsize(target) == info.file_size:
            continue
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with z.open(info) as src, open(target, 'wb') as out:
            while chunk := src.read(1 << 20):
                out.write(chunk)
        if os.path.getsize(target) != info.file_size:
            sys.exit(f'Incomplete file after unpacking: {name}')
for z in [n for n in os.listdir(dest) if n.endswith('.zip') and n != os.path.basename(archive)]:
    zipfile.ZipFile(os.path.join(dest, z)).extractall(dest)   # nested zips, if any
print('All files unpacked.')
PYCODE
  rm -f "$ZIP"
fi
DB="$(find_db)"
[[ -n "$DB" ]] || { echo "The database download does not contain all_cdms_final_official.csv."; exit 1; }
echo "Database: $DB"

step "5/5 Check that Astro AMASE loads (downloads its embedding model once)"
# Run from a neutral folder: inside .local-tools the source folder would shadow the package.
(cd / && "$VENV/bin/python" -c "import astro_amase; from astro_amase import assign_observations; print('Astro AMASE', astro_amase.__version__, 'is ready.')")

printf '{\n "python": "%s",\n "data_dir": "%s"\n}\n' "$VENV/bin/python" "$DB" > "$TOOLS/amase_paths.json"
echo
echo "Done. The app will use these paths automatically:"
echo "  Python:   $VENV/bin/python"
echo "  Database: $DB"
