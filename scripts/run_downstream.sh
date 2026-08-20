#!/bin/bash
# Per-group downstream pipeline: filter/dedup -> stats -> peaks+FRiP -> QC -> report.
# Run run_upstream.sh first (once, for the whole project).
#
# Usage: ./run_downstream.sh <group>      e.g. ./run_downstream.sh cutrun
#                                               ./run_downstream.sh atac
set -eu
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
GROUP="${1:?Usage: $0 <group>  (e.g. cutrun or atac)}"

"${SCRIPT_DIR}/run_filter_dedup.sh" "$GROUP"
"${SCRIPT_DIR}/build_stats.sh" "$GROUP"
"${SCRIPT_DIR}/run_macs2_frip.sh" "$GROUP"
"${SCRIPT_DIR}/run_qc_step2.sh" "$GROUP"
python3 "${SCRIPT_DIR}/make_report.py" "$GROUP"
echo "=== pipeline finished for group=$GROUP -> qc/${GROUP}/report.html ==="
