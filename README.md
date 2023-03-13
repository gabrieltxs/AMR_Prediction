##AMR prediction with Random Forest and LGBMClassifier models
This repository contains a Python script for predicting antimicrobial resistance (AMR) using the Random Forest and LGBMClassifier models. The script reads input datasets from a directory, applies feature extraction techniques to obtain k-mer features, trains and tests the models using cross-validation, and outputs the results in text files.

Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

Prerequisites
This script requires the following Python libraries:

pandas
scikit-learn
numpy
tqdm
lightgbm
hyperopt
joblib
bayesian-optimization
skopt
Installing
Clone the repository to your local machine and install the required libraries:

shell
Copy code
$ git clone https://github.com/username/repo.git
$ cd repo
$ pip install -r requirements.txt
Usage
To use the script, execute the following command:

css
Copy code
$ python main.py
Code Structure
The main script consists of several sections:

Import necessary libraries
Set seed for reproducibility
Define function to get list of models to evaluate
Load list of selected samples
Call function to get list of models
Initialize KFold cross-validation
Iterate over values of k to read the corresponding k-mer feature dataset
Iterate over the models list
Write results to text file
Data Description
The input datasets are CSV files containing bacterial genomic sequences and their corresponding resistance profiles for selected antibiotics. The script reads these files from a directory and applies k-mer feature extraction techniques to obtain numerical feature vectors.

Models
The script uses two models for AMR prediction: Random Forest and LGBMClassifier.

Output
The script outputs the results of each model to a text file in the specified output directory. The results include accuracy, precision, recall, F1 score, and area under the ROC curve.

Authors
John Doe - Initial work - johndoe
License
This project is licensed under the MIT License - see the LICENSE.md file for details.
