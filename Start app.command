#!/bin/bash
set -e
cd -- "$(dirname -- "$0")"

if [[ ! -x .venv/bin/python ]]; then
  echo "Create the project environment first; see Setup in README.md."
  read -r -p "Press Return to close. "
  exit 1
fi

if [[ -z "${SOFIA_EXECUTABLE:-}" && -x .local-tools/sofia ]]; then
  export SOFIA_EXECUTABLE="$PWD/.local-tools/sofia"
fi

echo "Open http://127.0.0.1:8501 in your browser."
echo "Keep this terminal open while using the app. Press Control-C to stop it."
exec .venv/bin/python -m streamlit run app.py --server.port 8501
