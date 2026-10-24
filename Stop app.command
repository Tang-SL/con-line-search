#!/bin/bash
# Double-click to stop con-line-search.
cd -- "$(dirname -- "$0")" || exit 1
scripts/launch.sh stop
sleep 2
