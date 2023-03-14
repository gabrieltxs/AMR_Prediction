import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import glob, os
from tqdm import tqdm


filename = '401-isolates.fasta'
path = 'D:\\OneDrive\\Documentos\\OneDrive\\Documentos\\Trabalho de Conclusão de Curso\\DATASET\\' + filename

fasta = open(path)
f= open('scaffold\\' + 'Pseudomonas aeruginosa strain CH4443.fna', 'w')
for item in fasta:
    if item[0] =='>':
        nContig = item.split('contig_')[1].split()[0]
        #print(nContig)
        #print(item)
        if nContig == '1':
            f.close()
            
            fnaName = item.split('[')[1].split('|')[0]
            f= open('scaffold\\' + fnaName.strip() + '.fna', 'w')
    else:
        f.write(item)