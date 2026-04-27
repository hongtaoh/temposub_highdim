#!/bin/bash
set -e  # Exit immediately on error

echo "run_mlhc.sh started at $(date)"
echo "Running in directory: $(pwd)"
echo "Running with arguments: $@"

# Prevent user-level site packages from interfering
export PYTHONNOUSERSITE=1

# ==============================================================================
# 📂 Prepare directories
# ==============================================================================
mkdir -p logs
chmod 755 logs
echo "Created logs directory at $(pwd)/logs"

mkdir -p algo_results

# ==============================================================================
# 🐍 Conda Env Extraction
# ==============================================================================
ENV_TARBALL="/staging/hhao9/env.tar.gz"
ENV_DIR=".conda_env"
PYTHON_EXEC=""

rm -rf "$ENV_DIR"

if [[ -f "$ENV_TARBALL" ]]; then
    echo "Extracting environment from $ENV_TARBALL..."
    mkdir -p "$ENV_DIR"
    tar -xzf "$ENV_TARBALL" -C "$ENV_DIR"
    PYTHON_EXEC="$ENV_DIR/bin/python"
    echo "Using extracted environment at $PYTHON_EXEC"
else
    echo "❌ $ENV_TARBALL not found — aborting"
    exit 1
fi


# ==============================================================================
# 🧪 Final sanity check
# ==============================================================================
echo "=== ENVIRONMENT VALIDATION ==="
echo "Python path: $PYTHON_EXEC"
echo "Python version: $($PYTHON_EXEC --version)"

if ! "$PYTHON_EXEC" -c "from kde_ebm import mixture_model; from pySuStaIn.MixtureSustain import MixtureSustain; import bebms;" &>/dev/null; then
    echo "❌ Final environment validation failed — aborting"
    exit 1
fi

# ==============================================================================
# 📦 Extract data
# ==============================================================================
DATA_TARBALL="/staging/hhao9/temposub_highdim_data.tar.gz"

if [[ -f "$DATA_TARBALL" ]]; then
    echo "📦 Extracting $DATA_TARBALL..."
    tar -xzf "$DATA_TARBALL"
    # If extraction creates "temposub_highdim_data", rename to "data"
    if [[ -d "temposub_highdim_data" ]]; then
        rm -rf data   # remove old data folder if it exists
        mv temposub_highdim_data data
        echo "Renamed temposub_highdim_data -> data"
    fi
else
    echo "❌ $DATA_TARBALL not found — aborting"
    exit 1
fi


# ==============================================================================
# See files present
# ==============================================================================
# echo "Files present:"
# ls -l

# ==============================================================================
# ▶️ Run Python Script
# ==============================================================================
echo "=== STARTING MAIN SCRIPT ==="
TQDM_DISABLE=1 "$PYTHON_EXEC" ./run_mlhc.py "$@"

# ==============================================================================
# 🧹 Cleanup pickle files for this run
# ==============================================================================
echo "Cleaning up pickle files for $1"

for d in algo_results/sustain_gmm/pickle_files algo_results/sustain_kde/pickle_files
do
    if [[ -d "$d" ]]; then
        # Use quotes around the path pattern in case of spaces
        rm -f "$d"/*"$1"* 2>/dev/null
        echo "Deleted files matching *$1* in $d"
    fi
done


echo "✅ Script completed at $(date)"
