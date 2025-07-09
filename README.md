
# MEFE 

Here, we introduce the Microbiome Elastic Feature Extraction (MEFE) algorithm, which leverages phylogenetic relationships, sequence similarity, and functional homology to elastically retrieve microbial candidate biomarkers. Evaluation on both synthetic and real microbiome datasets demonstrates that MEFE significantly reduces false positives and negatives compared to widely used methods, thereby improving microbiome status discrimination and disease detection.

## Project Structure

The project is structured as follows:

```
project/
│
├── data/                       # Input data files
├── database/                   # Database files
├── output/                     # Output files
├── main.py                     # Main entry point for running the algorithm
├── init.sh                     # Script to install dependencies from requirements.txt
├── run.sh                      # Script to run the main Python script
├── requirements.txt            # Required Python packages
└── README.md                   # Project documentation
```

## Features

- **Elastic Extraction Mechanism**: Refines the microbiome feature matrix using local attention weights based on the similarity of adjacent OTUs (Operational Taxonomic Units).

- **Statistical Analysis**: Uses the Wilcoxon test to identify significant biomarkers.

- **Data Preprocessing**: Handles missing values and outliers in the data.

- **Flexible Input and Output**: Can be easily modified to accept different datasets and outputs results in text files.

## Requirements

This project requires Python 3.x and the following Python libraries:

- `numpy`
- `pandas`
- `scikit-learn`
- `rpy2`

The dependencies can be installed by running:

```bash
pip install -r requirements.txt
```

## Installation

1. Move to the main folder:

    ```bash
    cd MEFE
    ```

2. Unzip the database files and install the required Python packages:

    ```bash
    source init.sh
    ```

3. Ensure that all input data is placed in the `data` folder.

## Usage

Once the dependencies are installed, you can run the `main.py` script by using the `run.sh` script, which will execute the pipeline:

```bash
 source run.sh
```

This will process the data files in the `data` folder, apply the MEFE algorithm, and save the results to the `output` folder.


### Output Files

- **biomarker_matrix.txt**: Contains the selected biomarker features after applying the elastic extraction mechanism.
- **biomarkers.txt**: Contains the p-values for the significant biomarkers selected by the Wilcoxon test.

### Citation

Guosen HOU, Yuxiao FEI, Hao GAO, Xiaoquan SU. MEFE: microbiome signature identification based on elastic feature extraction. Front. Comput. Sci., 2026, 20(1): 2001901 https://doi.org/10.1007/s11704-025-50323-1