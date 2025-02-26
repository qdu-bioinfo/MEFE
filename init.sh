#!/bin/bash

# Ensure that the script is being run in the project's root directory
cd "$(dirname "$0")"

# Check if the Python environment is already activated (for conda users)
# Uncomment the following line if using conda environment
# source activate your_env_name

# Unzip gg2_index.xz in the 'database' folder, try xz first, then tar, then 7z if both fail
if [ -f "database/gg2_index.txt.xz" ]; then
    echo "Attempting to extract gg2_index.txt.xz..."

    # First try using xz to extract
    if xz -d database/gg2_index.txt.xz; then
        echo "Successfully extracted using xz."
    else
        echo "xz extraction failed, trying tar..."

        # If xz fails, try using tar to extract
        if tar -xf database/gg2_index.txt.xz -C database/; then
            echo "Successfully extracted using tar."
        else
            echo "Tar extraction failed, trying 7z..."

            # If both xz and tar fail, try using 7z to extract
            if 7z x database/gg2_index.txt.xz -odatabase/; then
                echo "Successfully extracted using 7z."
            else
                echo "7z extraction failed. Please check the file and try again."
                exit 1
            fi
        fi
    fi
else
    echo "gg2_index.xz not found in the 'database' folder. Please ensure it exists."
    exit 1
fi

# Install the required Python packages listed in requirements.txt
if [ -f requirements.txt ]; then
    echo "Installing required packages from requirements.txt..."
    pip install -r requirements.txt
else
    echo "requirements.txt not found. Please ensure it exists in the current directory."
    exit 1
fi

# Print a message indicating the installation is complete
echo "Installation of required packages is complete."
