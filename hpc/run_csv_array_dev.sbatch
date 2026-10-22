#!/bin/bash
#SBATCH --job-name=mixmustr
#SBATCH --partition=cpuq
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=1-6400
#SBATCH --output=logs_dev_v1/%x_%A_%a.out
#SBATCH --error=logs_dev_v1/%x_%A_%a.err

set -euo pipefail

CSV="${1:?Usage: sbatch run_csv_array_dev.sbatch ./data/jobs_dev.csv}"
RSCRIPT="${RSCRIPT:-./worker_dev.R}"

mkdir -p logs_dev_v1

# count data rows (excluding header) to avoid running beyond EOF
NROWS=$(awk 'NR>1{n++} END{print n+0}' "$CSV")

TASK_ID="${SLURM_ARRAY_TASK_ID}"
if (( TASK_ID < 1 || TASK_ID > NROWS )); then
  echo "Task ${TASK_ID} out of range (1..${NROWS}); exiting."
  exit 0
fi

# Extract the (TASK_ID)-th data row (skip header). Robust for simple CSV (no embedded commas).
LINE=$(awk -v tid="$TASK_ID" 'NR==tid+1{print; exit}' "$CSV")

# Split CSV into fields
IFS=',' read -r tracers source_names degree sim_id model_id truth_stream G M input_file prop_seed tracer_seed stan_seed noises <<< "$LINE"

echo "==== SLURM ENVIRONMENT VARIABLES ===="
env | grep '^SLURM_' | sort
echo "======================================"

echo "CSV: $CSV"
echo "Task: $TASK_ID / $NROWS"

echo "tracers=$tracers"
echo "source_names=$source_names"
echo "degree=$degree"
echo "sim_id=$sim_id"
echo "model_id=$model_id"
echo "truth_stream=$truth_stream" 
echo "G=$G"
echo "M=$M"
echo "input_file=$input_file"
echo "prop_seed=$prop_seed"
echo "tracer_seed=$tracer_seed"
echo "stan_seed=$stan_seed"
echo "noises=$noises"

# Run R with command args
module load singularity
singularity exec -B .:/home/Project mixmustr_new.sif Rscript "$RSCRIPT" \
  --args \
  --tracers "$tracers" \
  --source_names "$source_names" \
  --degree "$degree" \
  --sim_id "$sim_id" \
  --model_id "$model_id" \
  --truth_stream "$truth_stream" \
  --G "$G" \
  --M "$M" \
  --input_file "$input_file" \
  --prop_seed "$prop_seed" \
  --tracer_seed "$tracer_seed" \
  --stan_seed "$stan_seed" \
  --noises "$noises"
