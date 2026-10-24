#!/bin/bash
# Double-click to start con-line-search. The first start sets up the Python
# environment. The app keeps running in the background after this window closes;
# use "Stop app.command" (or the Mac app) to stop it.
cd -- "$(dirname -- "$0")" || exit 1
if scripts/launch.sh start; then
  echo
  echo "The app keeps running after you close this window."
  echo "To stop it, double-click 'Stop app.command' or the CON Line Search app."
  sleep 3
else
  echo
  read -r -p "The app could not start (see above). Press Return to close. "
  exit 1
fi
