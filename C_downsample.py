# Importing the libraries
import numpy as np
import pandas as pd
import math 
from sklearn.utils import resample
# Suppressing warnings because of skopt verbosity
import warnings
warnings.filterwarnings("ignore")
import glob, os
import matplotlib.pyplot as plt

# Importing the dataset
# Get the current working directory
cwd = os.getcwd()

# Construct paths to input and output files relative to the current working directory
input_file_path = os.path.join(cwd, 'Python\\libFiles\\processed')
output_file_path = os.path.join(cwd, 'Python\\libFiles\\downsample')
input_kmer_path = os.path.join(cwd, 'Dataset')
output_results_path = os.path.join(cwd, 'Resultados')

# load list of selected samples
#antibiotic_dfs = {}
antibiotics = ['ceftazidime', 'ciprofloxacin', 'meropenem', 'tobramycin']

for antibiotic in antibiotics:
    # Load the antibiotic dataframe
    df = pd.read_csv(os.path.join(input_file_path, antibiotic+'_AMR.csv'), index_col=False, sep=';', header=0)
    print(antibiotic, df.iloc[:,1].value_counts())
    
    # Separate majority and minority classes
    df_ones = df[df.iloc[:,1]==1]
    df_zeros = df[df.iloc[:,1]==0]

    if len(df_ones)>= len(df_zeros):
        df_majority = df_ones
        df_minority = df_zeros
        maj = 'one'
    else:
        df_majority = df_zeros
        df_minority = df_ones
        maj = 'zero'
    
    # Calculate the proportion of the majority class
    proportion_majority = len(df_majority) / len(df)
    
    if maj == 'one' and proportion_majority < 0.6:
        # Downsample majority class
        #srt = math.sqrt(    4*(     (len(df_majority)/0.55)  - len(df)**2  )   )
        
        #x = int(( (-2*len(df) + 1/0.55) + srt )/2)
        df_minority_downsampled = resample(df_minority, replace=False, n_samples=int(len(df_minority)*0.5), random_state=42)
        df_downsampled = pd.concat([df_minority_downsampled, df_majority])
        print(df_downsampled.iloc[:,1].value_counts())
    elif maj == 'zero' :
        df_majority_downsampled = resample(df_majority, replace=False, n_samples=int(len(df_minority)*0.5), random_state=42)
        df_downsampled = pd.concat([df_majority_downsampled, df_minority])

        #df_downsampled = df
        print(df_downsampled.iloc[:,1].value_counts())
    else:
        df_downsampled = df

    df_downsampled = df_downsampled.sort_index(ascending=True)

    xis = pd.read_csv(os.path.join(input_kmer_path,'kmer3'+'.csv'), sep=';', header = 0)






    # Save antibiotic_df to a CSV file
    #antibiotic_df.reset_index(inplace=True)
    df_downsampled.to_csv(os.path.join(output_file_path, f'{antibiotic}_AMR.csv'), index=False, sep = ';')

    # EDA for the data
    counts = df_downsampled.iloc[:,1].value_counts()
    total_counts = counts.sum()

    # Set custom color scheme
    colors = ['#8c510a', '#01665e']

    # Set font sizes
    SMALL_SIZE = 14
    MEDIUM_SIZE = 16
    BIGGER_SIZE = 18

    # Set font family
    plt.rcParams["font.family"] = "Arial"

    fig, ax = plt.subplots(figsize=(8,6))
    ax.bar(['1','0'], counts.values, color=colors)

    # Add annotations to the bars
    for i, count in enumerate(counts.values):
        ax.annotate(f"{count}\n({count/total_counts*100:.1f}%)", xy=(i, count), ha='center', va='bottom', fontsize=SMALL_SIZE)

    ax.set_title(f'{antibiotic} Resistance (n={total_counts})', fontsize=BIGGER_SIZE, pad=20)
    ax.set_xlabel('Resistance Category', fontsize=MEDIUM_SIZE)
    ax.set_ylabel('Count', fontsize=MEDIUM_SIZE)
    ax.set_xticks([0, 1])
    #ax.set_xticklabels(['Resistant', 'Susceptible'], fontsize=SMALL_SIZE)

    # Add captions to the labels
    ax.text(0, -60, 'Resistant', fontsize=SMALL_SIZE, ha='center')
    ax.text(1, -60, 'Susceptible', fontsize=SMALL_SIZE, ha='center')

    # Set tick parameters
    ax.tick_params(axis='both', which='major', labelsize=SMALL_SIZE, length=8, width=2)

    # Set spine parameters
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Save the plot to the output file path
    plt.savefig(os.path.join(output_file_path, f'{antibiotic}_AMR.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    