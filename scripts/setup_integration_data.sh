#!/usr/bin/env bash
set -euo pipefail

# Download a SeqRepo snapshot used by integration tests and print export command.

SEQREPO_ROOT="${SEQREPO_ROOT:-/usr/local/share/seqrepo}"
SEQREPO_INSTANCE="${SEQREPO_INSTANCE:-2024-12-20}"
NO_SUDO="${NO_SUDO:-0}"

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage: scripts/setup_integration_data.sh

Environment variables:
  SEQREPO_ROOT      Installation root for SeqRepo snapshots.
                    Default: /usr/local/share/seqrepo
  SEQREPO_INSTANCE  SeqRepo instance to download.
                    Default: 2024-12-20
  NO_SUDO           Set to 1 to skip sudo for mkdir/chmod.
                    Default: 0

After completion, export GA4GH_VRS_DATAPROXY_URI as printed by this script.
EOF
  exit 0
fi

if [[ "${NO_SUDO}" == "1" ]]; then
  mkdir -p "${SEQREPO_ROOT}"
  chmod a+w "${SEQREPO_ROOT}"
else
  sudo mkdir -p "${SEQREPO_ROOT}"
  sudo chmod a+w "${SEQREPO_ROOT}"
fi

uv run seqrepo --root-directory "${SEQREPO_ROOT}" pull --instance-name "${SEQREPO_INSTANCE}"

echo
echo "Set this before running integration tests:"
echo "export GA4GH_VRS_DATAPROXY_URI=seqrepo+file://${SEQREPO_ROOT}/${SEQREPO_INSTANCE}"

