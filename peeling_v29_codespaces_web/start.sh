#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "--background" ]]; then
  if pgrep -f "streamlit run app.py" >/dev/null 2>&1; then
    exit 0
  fi
  nohup streamlit run app.py >/tmp/peeling-streamlit.log 2>&1 &
  exit 0
fi

exec streamlit run app.py
