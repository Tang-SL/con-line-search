#!/bin/bash
# Double-click (or use the CON Line Agent app) to start the line-identification agent:
# Claude Code, opened in this folder, where the line-id-agent skill and scripts/agent live.
cd -- "$(dirname -- "$0")" || exit 1
folder="$(pwd)"
mkdir -p "$HOME/Library/Application Support/con-line-search" && printf '%s' "$folder" > "$HOME/Library/Application Support/con-line-search/project-folder"

printf '\n  \033[1mCON Line Agent\033[0m: line identification in ALMA cubes\n  %s\n\n' "$folder"
if [[ ! -x .venv/bin/python ]]; then
  echo "  Note: the Python environment (.venv) is not set up yet. Start the CON Line Search app once, or let"
  echo "  the agent set it up (it runs scripts/launch.sh)."
fi
if [[ -z "${ADS_API_TOKEN:-}${ADS_DEV_KEY:-}" && ! -s "$HOME/.ads/dev_key" && ! -s .local-tools/ads_token ]] && \
   ! grep -qs 'ADS_API_TOKEN' "$HOME/.zshrc" "$HOME/.zprofile" "$HOME/.bash_profile"; then
  echo "  Note: no NASA ADS token found (~/.ads/dev_key): the agent cannot search the literature."
fi
cat <<'TXT'
  Ask in plain words, for example:
    What lines are in IRAS 17208-0014? The data are in ~/Documents/CONs/data/IRAS17208-0014
    Rerun the identification strategy for runs/IC860_B5_20261105

TXT

# Start Claude Code through the login shell, so that PATH and ADS_API_TOKEN from ~/.zshrc apply.
shell="${SHELL:-/bin/zsh}"
exec "$shell" -lic '
  claude_bin="$(command -v claude 2>/dev/null)"
  for c in "$HOME/.local/bin/claude" "$HOME/.claude/local/claude" /opt/homebrew/bin/claude /usr/local/bin/claude "$HOME/.npm-global/bin/claude"; do
    [[ -z "$claude_bin" && -x "$c" ]] && claude_bin="$c"
  done
  if [[ -z "$claude_bin" ]]; then
    echo "  Claude Code is not installed (the claude command was not found)."
    echo "  Opening the installation guide; afterwards start CON Line Agent again."
    open "https://docs.claude.com/en/docs/claude-code"
    read "?  Press Return to close. " 2>/dev/null || read -r -p "  Press Return to close. "
    exit 1
  fi
  exec "$claude_bin"
'
