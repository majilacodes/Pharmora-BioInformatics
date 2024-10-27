# Pharmora - Advanced Bioinformatics Research Platform

Welcome to **Pharmora**, an all-in-one bioinformatics research platform designed to accelerate drug discovery and molecular research. Pharmora provides a comprehensive set of tools for researchers and scientists working in the fields of cheminformatics, bioinformatics, and pharmaceutical research.

## Overview

Pharmora simplifies the drug discovery process by offering cutting-edge models that enable the analysis, prediction, and screening of drug-like molecules. By leveraging computational techniques, Pharmora helps users optimize lead compounds, predict molecular properties, and explore novel drug candidates, all within a user-friendly interface.

## Key Features and Models

Pharmora integrates **seven specialized models** to assist with various aspects of bioinformatics and drug discovery:

### 1. **Bioactivity Predictor (pIC50)**  
   This model predicts the bioactivity of small molecules against specific biological targets. Using molecular descriptors and machine learning algorithms, it calculates the half-maximal inhibitory concentration (pIC50) values. This helps researchers prioritize molecules with higher predicted biological activity.
   
   **Input:** Chemical structure (SMILES format or molecular file)  
   **Output:** Predicted pIC50 value
   
   - *Use Case:* Drug discovery campaigns seeking compounds with strong binding affinities to target proteins.

### 2. **New Molecule Discovery for Specific Targets**  
   This tool allows users to discover new molecules that are likely to interact with specific biological targets. By analyzing large compound libraries and applying similarity metrics, this model identifies novel candidates that may have therapeutic potential.

   **Input:** Target protein identifier or structure  
   **Output:** List of potential drug-like molecules
   
   - *Use Case:* Early-stage drug discovery focusing on generating new leads for target proteins.

### 3. **Solubility Predictor**  
   The solubility of a drug is crucial for its absorption and distribution in the body. Pharmora’s Solubility Predictor estimates the water solubility of compounds based on their molecular structure. This helps researchers assess whether a molecule has favorable drug-like properties.

   **Input:** Chemical structure  
   **Output:** Predicted solubility value (in mol/L)
   
   - *Use Case:* Screening drug candidates for formulation and bioavailability studies.

### 4. **DNA Nucleotide Counter**  
   This model analyzes DNA sequences and provides detailed nucleotide composition, helping researchers with genetic studies, molecular cloning, and sequence analysis.

   **Input:** DNA sequence  
   **Output:** Nucleotide count (A, T, C, G)
   
   - *Use Case:* Genetic analysis and sequence optimization in synthetic biology projects.

### 5. **Antimicrobial Activity Predictor for Peptides**  
   With the rise of antimicrobial resistance, peptides have emerged as promising therapeutic agents. This model predicts the antimicrobial activity of peptides, helping researchers identify potential candidates for therapeutic development.

   **Input:** Peptide sequence  
   **Output:** Predicted antimicrobial activity score
   
   - *Use Case:* Designing antimicrobial peptides for use in infections, wound healing, or novel therapeutics.

### 6. **Molecular Descriptor Calculator**  
   Molecular descriptors are vital for cheminformatics and QSAR (Quantitative Structure-Activity Relationship) studies. Pharmora's Molecular Descriptor Calculator generates a wide array of descriptors such as molecular weight, logP, polar surface area, and more.

   **Input:** Chemical structure  
   **Output:** Comprehensive list of molecular descriptors
   
   - *Use Case:* Data-driven drug discovery and chemical library analysis for QSAR modeling.

### 7. **Lipinski’s Rule of Five Filter**  
   Lipinski's Rule of Five is a widely accepted guideline to evaluate drug-likeness. This model automatically filters compounds to ensure they meet these criteria, such as molecular weight, hydrogen bond donors/acceptors, and logP, ensuring that the compounds are likely to be orally active.

   **Input:** Chemical structure  
   **Output:** Pass/Fail on Lipinski's Rule of Five
   
   - *Use Case:* Early-stage filtering of chemical libraries to focus on drug-like molecules.

## Intel® Extension for Scikit-learn

In this project, we utilized the **Intel® Extension for Scikit-learn** to enhance the performance of machine learning models. This extension optimizes scikit-learn for Intel® architectures, allowing us to speed up the training and inference of our models without changing the codebase.

### Why Intel® Extension for Scikit-learn?

- **Performance Boost:** The extension allows us to achieve higher computational efficiency, particularly for large datasets or complex models. We experienced faster training and inference times, enabling quicker iterations during the drug discovery process.
  
- **Seamless Integration:** The extension integrates smoothly into the existing scikit-learn ecosystem, allowing us to leverage its optimizations with minimal code modifications.

### How We Used Intel’s scikit-learn

We applied Intel® Extension for Scikit-learn in several models, including our **Bioactivity Predictor** and **Molecular Descriptor Calculator**, to accelerate computations involving large molecular datasets.

Here’s an example of how we patched scikit-learn in our code:

python
import numpy as np
import dpctl
from sklearnex import patch_sklearn, config_context
patch_sklearn()

# Example of using patched scikit-learn models
from sklearn.ensemble import RandomForestRegressor

# Using Intel-optimized RandomForestRegressor
model = RandomForestRegressor()
model.fit(X_train, y_train)
predictions = model.predict(X_test)


## Installation

To use Pharmora on your local machine, follow the steps below to install the required dependencies and run the platform.

### Prerequisites

- Python 3.8 or above
- Pip (Python package installer)

### Installation Steps

1. Clone the Pharmora repository to your local machine:

    ```bash
    git clone https://github.com/your-username/pharmora.git
    ```

2. Navigate to the project directory:

    ```bash
    cd pharmora
    ```

3. Install the required Python libraries:

    ```bash
    pip install -r requirements.txt
    ```

4. Launch the Streamlit app to interact with the platform:

    ```bash
    streamlit run app.py
    ```

Once the app is running, open your browser and navigate to the provided local address (typically `http://localhost:8501`) to start using the Pharmora interface.

## Usage

### Interface Overview
- **Home**: Get an overview of the platform and access quick links to the different tools.
- **Molecular Tools**: Select from a range of drug discovery models to start analysis or prediction on your compounds.
- **Upload Data**: Upload molecular files (in formats like SMILES, FASTA, or CSV) for analysis.
- **Visualize Results**: View the predicted results, bioactivity, molecular descriptors, and more in a tabular or graphical format.

### Example Usage

To predict the pIC50 value of a molecule, follow these steps:

1. Select **Bioactivity Predictor** from the sidebar.
2. Upload the chemical structure file (SMILES format or molecular file).
3. Click on **Predict** to get the pIC50 value for the compound.
4. View the result in the results section.

## Contributing

We are constantly improving Pharmora and welcome contributions from the community. Whether it's fixing bugs, adding new features, or improving documentation, your help is appreciated!

### How to Contribute:

1. Fork the repository on GitHub.
2. Create a new branch for your feature or fix.
3. Submit a pull request with a detailed explanation of your changes.

## Acknowledgments

Pharmora was inspired by the need for advanced, user-friendly bioinformatics tools that streamline drug discovery and molecular research. We thank the open-source community and bioinformatics researchers for their contributions to this field.

---

Thank you for choosing Pharmora – Your Comprehensive Partner in Drug Discovery and Molecular Research.
