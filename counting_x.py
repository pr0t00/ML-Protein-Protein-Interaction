# %%

from Bio import SeqIO  # https://biopython.org/wiki/SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import os, sys
import numpy as np
from scipy import stats
import seaborn as sns
sns.set(color_codes=True)
import Bio  # https://biopython.org/wiki/SeqIO


# this_dir = os.path.dirname(this_file_path)
dir_path = "/Users/universal/Google Drive/GOOGLE DRIVE/UNIVERSITY/BIOINFORMATICS/BIOINFO - PBL/ML-Protein-Protein-Interaction/"

print(dir_path)
seq_path = os.path.join(dir_path, "data", "PP_step1_trn.fas")
label_path = os.path.join(dir_path, "data", "PP_step1_trn.class")

# Importing the amino acids
prot_list = []

x = 0

for record in SeqIO.parse(seq_path, "fasta"):
    seq = str(record.seq)
    for char in seq:
        if char == "X":
            x += 1
print(x)