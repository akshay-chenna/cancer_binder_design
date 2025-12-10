export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

SCRIPT_TO_RUN="af2-1_2.sh"
TIME_LIMIT="60m"
GRACE_PERIOD="10s"

echo "=================================================="
echo "▶️ Starting loop for ${SCRIPT_TO_RUN}."
echo "The script will run in cycles of max ${TIME_LIMIT} until it completes."
echo "=================================================="

# Loop until the script finishes successfully.
while true; do
    echo "--> Starting a new ${TIME_LIMIT} attempt at $(date)"

    timeout -k ${GRACE_PERIOD} ${TIME_LIMIT} bash "${SCRIPT_TO_RUN}"
    EXIT_STATUS=$?

    if [ ${EXIT_STATUS} -eq 0 ]; then
        # Success! The script finished on its own.
        echo "✅ SUCCESS: Script completed within the time limit."
        break # Exit the loop.

    elif [ ${EXIT_STATUS} -eq 124 ]; then
        # Timed out. Loop will continue.
        echo "🔄 TIMEOUT: Script was killed. Restarting..."

    else
        # An unexpected error occurred.
        echo "⚠️ ERROR: Script failed with exit status ${EXIT_STATUS}. Stopping."
        break # Exit the loop to avoid repeating an error.
    fi
done

echo "=================================================="
echo "Process finished."
echo "=================================================="
