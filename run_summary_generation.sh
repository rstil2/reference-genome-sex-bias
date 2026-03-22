#!/bin/zsh

# This script is designed to be a foolproof way to generate the summary results
# for Project 33, using absolute paths to avoid any ambiguity about the
# execution context.

set -e # Exit immediately if a command exits with a non-zero status.
set -x # Print commands and their arguments as they are executed.

# 1. Define the absolute path for the project's root directory.
PROJECT_DIR="/Users/stillwell/Documents/Google Drive/Project 33 - Bias in Reference Genomes/Study_v2_Real_Data/"
echo "--- Project Directory: $PROJECT_DIR"

# 2. Change to the project directory. This is critical.
cd "$PROJECT_DIR"
echo "--- Changed current directory to: $(pwd)"

# 3. Define the absolute path for the Python script to be executed.
PYTHON_SCRIPT_PATH="$PROJECT_DIR/Code/12_generate_summary_absolute.py"
echo "--- Python script to execute: $PYTHON_SCRIPT_PATH"

# 4. Verify that the Python script actually exists before trying to run it.
if [ ! -f "$PYTHON_SCRIPT_PATH" ]; then
    echo "--- ERROR: The script does not exist at the specified path: $PYTHON_SCRIPT_PATH"
    exit 1
fi
echo "--- Verified script exists."

# 5. Execute the Python script using its absolute path.
# The python script itself also uses absolute paths for all file operations.
echo "--- Executing Python script..."
python3 "$PYTHON_SCRIPT_PATH"
echo "--- Python script execution finished."

# 6. Define the absolute path for the output directory.
FIGURES_DIR="$PROJECT_DIR/Figures/"
echo "--- Output directory: $FIGURES_DIR"

# 7. List the contents of the Figures directory to verify the output files.
echo "--- Verifying output files in $FIGURES_DIR..."
ls -l "$FIGURES_DIR"
echo "--- Verification complete."
