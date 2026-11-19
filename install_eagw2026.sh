#!/usr/bin/env bash
# Create/update the Conda environment used by this repository.
#
# Usage:
#   bash install_eagw2026.sh
#
# The active package list lives in installtions.txt.  Entries commented out
# there (currently lalsuite and teobresums) remain intentionally optional.

set -euo pipefail

ENV_NAME="eagw2026"
PYTHON_VERSION="3.11"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REQUIREMENTS_FILE="${SCRIPT_DIR}/installtions.txt"

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: Conda was not found on PATH." >&2
    echo "Install Miniconda or Anaconda, then run this script again." >&2
    exit 1
fi

if [[ ! -f "${REQUIREMENTS_FILE}" ]]; then
    echo "Error: requirements file not found: ${REQUIREMENTS_FILE}" >&2
    exit 1
fi

cd "${SCRIPT_DIR}"

if conda run --no-capture-output -n "${ENV_NAME}" python --version >/dev/null 2>&1; then
    echo "Updating existing Conda environment '${ENV_NAME}'..."
else
    echo "Creating Conda environment '${ENV_NAME}' with Python ${PYTHON_VERSION}..."
    conda create --yes --name "${ENV_NAME}" "python=${PYTHON_VERSION}" pip
fi

echo "Upgrading pip build tooling..."
conda run --no-capture-output -n "${ENV_NAME}" \
    python -m pip install --upgrade pip setuptools wheel

echo "Installing packages listed in ${REQUIREMENTS_FILE}..."
conda run --no-capture-output -n "${ENV_NAME}" \
    python -m pip install --upgrade -r "${REQUIREMENTS_FILE}"

echo
echo "Installation complete. Activate the environment with:"
echo "  conda activate ${ENV_NAME}"
