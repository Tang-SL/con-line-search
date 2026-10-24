#!/bin/bash
# Start, open, stop or check the con-line-search app.
#
#   scripts/launch.sh start    set up on first use, start in the background, open the browser
#   scripts/launch.sh open     open the browser on the running app
#   scripts/launch.sh stop     stop the app
#   scripts/launch.sh status   exit 0 if the app is running
#
# The server runs in the background (log: .local-tools/app.log), so no terminal
# window has to stay open. CONLINES_PORT changes the port (default 8501).
set -u
ROOT="$(cd -- "$(dirname -- "$0")/.." && pwd)"
cd -- "$ROOT" || exit 1
PORT="${CONLINES_PORT:-8501}"
URL="http://127.0.0.1:${PORT}"
STATE="$ROOT/.local-tools"
LOG="$STATE/app.log"
PIDFILE="$STATE/app.pid"
PY="$ROOT/.venv/bin/python"
STAMP="$STATE/app.code-stamp"
mkdir -p "$STATE"

say() { printf '%s\n' "$*"; }

healthy() {
  curl -fs --max-time 2 "$URL/_stcore/health" 2>/dev/null | grep -q ok
}

open_url() {
  if command -v open >/dev/null 2>&1; then open "$URL"
  elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$URL" >/dev/null 2>&1
  else say "Open $URL in your browser."
  fi
}

remember_location() {
  # Lets a copy of the Mac app outside this folder find the project.
  local dir="$HOME/Library/Application Support/con-line-search"
  [[ "$(uname)" == Darwin ]] || dir="${XDG_CONFIG_HOME:-$HOME/.config}/con-line-search"
  mkdir -p "$dir" 2>/dev/null && printf '%s\n' "$ROOT" > "$dir/project-folder" 2>/dev/null
}

find_python() {
  local candidate
  for candidate in python3.12 python3.13 python3 \
      /Library/Frameworks/Python.framework/Versions/3.12/bin/python3.12 \
      /opt/homebrew/bin/python3.12 /usr/local/bin/python3.12; do
    if command -v "$candidate" >/dev/null 2>&1 && \
       "$candidate" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 12) else 1)' 2>/dev/null; then
      command -v "$candidate"; return 0
    fi
  done
  return 1
}

setup() {
  if [[ ! -x "$PY" ]]; then
    local base
    if ! base="$(find_python)"; then
      say "Python 3.12 or newer is needed to set up the app (https://www.python.org/downloads/)."
      return 1
    fi
    say "First start: creating the project environment with $base"
    "$base" -m venv "$ROOT/.venv" || return 1
  fi
  if ! "$PY" -c 'import streamlit, plotly, adplib' >/dev/null 2>&1; then
    say "Installing the app's Python packages (a few minutes, needs internet)…"
    "$PY" -m pip install --upgrade pip >/dev/null || return 1
    "$PY" -m pip install -r "$ROOT/requirements-ui.txt" || return 1
  fi
  return 0
}

sofia_check() {
  if [[ -z "${SOFIA_EXECUTABLE:-}" && -x "$STATE/sofia" ]]; then
    export SOFIA_EXECUTABLE="$STATE/sofia"
  fi
  if [[ -z "${SOFIA_EXECUTABLE:-}" ]] && ! command -v sofia >/dev/null 2>&1; then
    say "Note: SoFiA 2 was not found. The app opens, but detection needs it:"
    say "      put the compiled 'sofia' in .local-tools/ (see README, Setup)."
  fi
}

code_stamp() {
  # Changes whenever the app's Python code changes (pulls, branch switches, edits).
  "$PY" - "$ROOT" <<'PYCODE' 2>/dev/null
import hashlib, pathlib, sys
root = pathlib.Path(sys.argv[1])
files = [root/'app.py', *sorted((root/'conlines').rglob('*.py'))]
h = hashlib.sha256()
for f in files:
    if f.exists():
        h.update(str(f.relative_to(root)).encode()); h.update(f.read_bytes())
print(h.hexdigest())
PYCODE
}

running_pid() {
  local pid
  [[ -f "$PIDFILE" ]] || return 1
  pid="$(cat "$PIDFILE" 2>/dev/null)"
  [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null && printf '%s' "$pid"
}

start() {
  remember_location
  if healthy; then
    if [[ -x "$PY" && -f "$STAMP" && "$(cat "$STAMP")" == "$(code_stamp)" ]]; then
      say "The app is already running at $URL"
      open_url; return 0
    fi
    # Streamlit keeps imported modules in memory, so a server started before the
    # code changed would mix old and new code. Restart it.
    say "The running app was started with older code; restarting it."
    stop || return 1
  fi
  setup || return 1
  sofia_check
  say "Starting the app…"
  # New session (setsid): closing the Terminal window that started it neither
  # stops the server nor asks about running processes.
  nohup "$PY" -c 'import os, sys; os.setsid(); os.execv(sys.argv[1], sys.argv[1:])' \
        "$PY" -m streamlit run "$ROOT/app.py" --server.port "$PORT" --server.headless true \
        > "$LOG" 2>&1 < /dev/null &
  local pid=$!
  printf '%s\n' "$pid" > "$PIDFILE"
  code_stamp > "$STAMP"
  disown "$pid" 2>/dev/null || true
  local i
  for i in $(seq 1 120); do
    if healthy; then
      say "The app is running at $URL"
      open_url; return 0
    fi
    if ! kill -0 "$pid" 2>/dev/null; then
      say "The app stopped while starting. Last lines of $LOG:"
      tail -n 15 "$LOG"; return 1
    fi
    sleep 1
  done
  say "The app did not answer within two minutes. See $LOG"
  return 1
}

stop() {
  local pid
  if pid="$(running_pid)"; then
    kill "$pid" 2>/dev/null
    for _ in $(seq 1 10); do kill -0 "$pid" 2>/dev/null || break; sleep 1; done
    kill -0 "$pid" 2>/dev/null && kill -9 "$pid" 2>/dev/null
  elif healthy && command -v lsof >/dev/null 2>&1; then
    pid="$(lsof -ti "tcp:$PORT" -sTCP:LISTEN 2>/dev/null | head -n 1)"
    [[ -n "$pid" ]] && kill "$pid" 2>/dev/null
  fi
  rm -f "$PIDFILE" "$STAMP"
  sleep 1
  if healthy; then say "Could not stop the app on port $PORT."; return 1; fi
  say "The app is stopped."
}

case "${1:-start}" in
  start) start ;;
  open) if healthy; then open_url; else start; fi ;;
  stop) stop ;;
  status) healthy ;;
  *) say "Usage: $0 [start|open|stop|status]"; exit 2 ;;
esac
