#!/bin/bash

# Activate the conda environment if needed (uncomment if using conda environment)
# source activate your_env_name

# Navigate to the project directory (optional, if you are already in the project folder)
cd "$(dirname "$0")"

# Run the Python script using Python 3
python3 main.py

# Print a message indicating the script has finished running
echo "main.py has been executed successfully."
