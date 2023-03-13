# Importing necessary libraries
import numpy as np # for numerical operations
import pandas as pd # for data processing and manipulation
import glob, os # for handling files and directories
from tqdm import tqdm # for progress bar visualization

# Defining a function to generate k-mers
def cMer(Tseq, k):
    # Creating a dictionary to store the frequency of k-mers
    kFreq = {}
    
    # Iterating through the entire sequence of nucleotides and generating k-mers
    for i in range(0, len(Tseq) - k + 1):
        kmer = Tseq[i:i +k]
        
        # Updating the frequency count for each k-mer in the dictionary
        if kmer in kFreq:
            kFreq[kmer] += 1
        else:
            kFreq[kmer] = 1
                
    # Returning the dictionary of k-mer frequencies
    return kFreq   

# Setting the current directory to the path of the dataset folder
os.chdir("D:\\OneDrive\\Documentos\\OneDrive\\Documentos\\Trabalho de Conclusão de Curso\\DATASET\\scaffold")

# Initializing an empty list to store the names of the files in the dataset folder
listaN = []

# Looping through the range of values [2, 3, 4] to generate 2-mer, 3-mer, and 4-mer frequencies
for i in range(2,5):
    # Initializing an empty pandas DataFrame to store the k-mer frequencies of all the files
    fullDF = pd.DataFrame()
    
    # Initializing a counter variable to track the number of files processed
    cont = 0
    
    # Looping through all the files in the dataset folder
    for file in tqdm(glob.glob("*")):
        
        # Printing the name of the current file
        print(file)
        
        # Opening the current file and reading its contents
        TargetFile = open(file)
        
        # Appending the name of the current file to the list of file names
        listaN.append(file)
        
        # Creating a new row in the pandas DataFrame to store the k-mer frequencies of the current file
        dfRow = pd.DataFrame(index = listaN)
        
        ################################################################
        # Generating the k-mer frequencies of the current file
        next(TargetFile)
        readTF = TargetFile.read()
        Tseq = "".join(readTF.split())
        rf = cMer(Tseq, i)
        rf = dict(sorted(rf.items()))
        ################################################################
        
        # Checking if the current k-mer is already present in the pandas DataFrame
        chave = fullDF.keys()
        kmerEX = [] 
        for k in rf.keys():
            if k not in chave:
                kmerEX.append(k)
        
        # Adding a new column for the current k-mer in the pandas DataFrame
        df = pd.DataFrame(columns = kmerEX)
        fullDF = pd.concat([fullDF,df], axis = 1)
        fullDF = pd.concat([fullDF, dfRow], axis=0)
        
        # Adding the k-mer frequencies of the current file to the pandas DataFrame
        for j in rf.keys():
            fullDF.at[file,j] = int(rf[j])
        ################################################################
        
        # Closing the current file
        TargetFile.close()
        
        # Incrementing the counter variable and removing the name of the current file from the list of file names
        cont = cont +1
        listaN.pop()    
        
    # Saving the pandas DataFrame as a CSV file for the current value of k
    
    fullDF.to_csv("D:\\OneDrive\\Documentos\\OneDrive\\Documentos\\Trabalho de Conclusão de Curso\\DATASET\\results\\kmer"+str(i)+".csv", sep = ';')

