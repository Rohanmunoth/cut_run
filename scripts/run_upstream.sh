#!/bin/bash
# Shared upstream steps, run ONCE for the whole project (both groups' fastqs
# together): organize by assay -> trim -> align. Downstream steps are
# per-group -- see run_downstream.sh.
set -eu
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"${SCRIPT_DIR}/organize_by_assay.sh"
"${SCRIPT_DIR}/run_trim_galore.sh"
"${SCRIPT_DIR}/run_bowtie2.sh"
echo "=== upstream finished ==="
