#!/bin/bash
# Double-click to install Astro AMASE for the "Astro AMASE" identification method.
cd -- "$(dirname -- "$0")" || exit 1
if scripts/setup_amase.sh; then
  echo
  echo "Back in the app, open Step 4 › Astro AMASE settings: the paths are filled in."
  read -r -p "Press Return to close. "
else
  echo
  read -r -p "Setup stopped (see above). Fix the problem and run it again; finished steps are skipped. Press Return to close. "
  exit 1
fi
