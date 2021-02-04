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


# get the right path
import inspect
# this_file_path = os.path.abspath(inspect.getframeinfo(inspect.currentframe()).filename)
# this_dir = os.path.dirname(this_file_path)
dir_path = "/Users/universal/Google Drive/GOOGLE DRIVE/UNIVERSITY/BIOINFORMATICS/BIOINFO - PBL/ML-Protein-Protein-Interaction/"

print(dir_path)
seq_path = os.path.join(dir_path, "data", "PP_step1_trn.fas")
label_path = os.path.join(dir_path, "data", "PP_step1_trn.class")

# Importing the amino acids
prot_list = []


for record in SeqIO.parse(seq_path, "fasta"):
    seq = str(record.seq)
    prot_list.append(seq)
prot_series = pd.Series(prot_list)


seq_len = []
for element in prot_list:
    seq_len.append(len(element))
len_series = pd.Series(seq_len)

# Importing the labels
label_list = []
for record in SeqIO.parse(label_path, "fasta"):
    seq = str(record.seq)
    label_list.append(seq)
label_series = pd.Series(label_list)



# Combining length series and labels into one data frame
data = {"Labels" : label_list, "Length" : seq_len, "Sequence" : prot_list}
prot_df = pd.DataFrame(data, columns = ["Labels", "Length", "Sequence"])

# changing the bind/nonbind labels to numerical values
# prot_df["Labels"] = prot_df["Labels"].apply(lambda x: 1 if x == "Bind" else 0)

print(prot_df.count())

# Sorting out all the duplicate entries
prot_df = prot_df.drop_duplicates().reset_index()
prot_df = prot_df.replace("X", "G")


print(prot_df.count())


# print(label_series.describe()) # Non_Bind: 651; Bind: 536



# %%
#print(os.chdir("/Users/universal/Google\ Drive/GOOGLE\ DRIVE/UNIVERSITY/BIOINFORMATICS/BIOINFO\ -\ PBL/ML-Protein-Protein-Interaction/plots"))


# Hypothesis 1: protein length is in relation to binding properties
#prot_df["Labels"] = prot_df["Labels"].apply(lambda x: "Interact" if x == "Bind" else "Non-interact")
plt.figure
ax = sns.violinplot(x="Labels", y="Length", data=prot_df)
plt.savefig("violinplot_len-and-labels", dpi=300)
plt.show()

# %% calculate avg length of binding and non binding prot
# prot_df

# interacting
bind = prot_df.loc[prot_df["Labels"] == "Bind"]
print(bind["Length"].mean())

# non interacting
non_bind = prot_df.loc[prot_df["Labels"] == "Non_Bind"]
print(non_bind["Length"].mean())

# %%
# Hypothesis 2: frequency of amino acids influences binding properties

dict_list = []


for i in range(len(prot_df)):
    seq = prot_df.iat[i,3]
    prot_dict = {"C": 0, "D": 0, "S": 0, "Q": 0, "K": 0, "I": 0, "P": 0, "T": 0, "F": 0, "N": 0, "G": 0, "H": 0, "L": 0,
                 "R": 0, "W": 0, "A": 0, "V": 0, "E": 0, "Y": 0, "M": 0, "X": 0, "Z": 0}
    for char in seq:
        prot_dict[char] += 1
    prot_dict["G"] += prot_dict["X"]
    prot_dict["X"] = 0
    prot_dict["G"] += prot_dict["Z"]
    prot_dict["Z"] = 0
    dict_list.append(prot_dict)

dict_list


# creating a new df for the freq_dicts to then append it to the master df
freq_df = pd.DataFrame()

for dict_dict in dict_list:
    freq_df = freq_df.append(dict_dict, ignore_index=True)

freq_df = freq_df.drop(["X"], axis=1)
freq_df = freq_df.drop(["Z"], axis=1)

print(freq_df)


prot_df = pd.concat([prot_df, freq_df], axis=1, join_axes=[prot_df.index])
print(prot_df)



# for dict_dict in dict_list:
#     for k, v in dict_dict.items():
#         if (k == "X" and v > 0):
#             print(v)


# %%


# %% Melting the dataframe to be able to plot things

temp_df = prot_df.drop(["index", "Length", "Sequence"], axis=1)

df = pd.melt(temp_df, id_vars="Labels")
df
df = df.dropna(axis=0, how="any", thresh=None, subset=None, inplace=False)

print(len(df))
print(df.iat[18099,0])

df


# %%

sns.set(rc={'figure.figsize':(20,8.27)})
plt.rcParams["axes.labelsize"] = 20
sns.set(font_scale = 2)

df["Labels"] = df["Labels"].apply(lambda x: "Interact" if x == "Bind" else "Non-interact")

ax = sns.violinplot(width=0.9, gridsize=100, x="variable", y="value", hue="Labels", data=df, palette="muted", split=True)
ax.set(xlabel='Amino-acids', ylabel='Frequency')
plt.setp(ax.get_legend().get_texts(), fontsize='17')
plt.savefig("amino-acids-frequency.png", dpi=300)
plt.show()


# %%


x_names = ["C", "D", "S", "Q", "K", "I", "P", "T", "F", "N", "G", "H", "L", "R", "W", "A", "V", "E", "Y", "M"]

#ax = sns.violinplot(x=df.index, y=prot_df[["A", "B"]], hue="Labels", data=prot_df, palette="muted", split=True, size=10)

plt.show()
print("here")

# %%

df = pd.DataFrame({'X_Axis':[1,3,5,7,10,20],
                   'col_2':[.4,.5,.4,.5,.5,.4],
                   'col_3':[.7,.8,.9,.4,.2,.3],
                   'col_4':[.1,.3,.5,.7,.1,.0],
                   'col_5':[.5,.3,.6,.9,.2,.4]})

print (df)

df = df.melt('X_Axis', var_name='cols',  value_name='vals')

print(df)
g = sns.catplot(x="X_Axis", y="vals", hue='cols', data=df)



# %%
# frequency_lists = [["C"], ["D"], ["S"], ["Q"], ["K"], ["I"], ["P"], ["T"], ["F"], ["N"], ["G"], ["H"], ["L"], ["R"], ["W"], ["A"], ["V"], ["E"], ["Y"], ["M"], ["X"], ["Z"]]
frequency_lists = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
for dict in frequency_dicts:
    list = []
    counter = 0
    for value in dict.values():
        frequency_lists[counter].append(value)
        counter += 1

print(frequency_lists)
# %%


# Combining frequency lists into one data frame
data = {"C Frequency" : frequency_lists[0], "D Frequency" : frequency_lists[1], "S Frequency" : frequency_lists[2], "Q Frequency" : frequency_lists[3], "K Frequency" : frequency_lists[4], "I Frequency" : frequency_lists[5], "P Frequency" : frequency_lists[6], "T Frequency" : frequency_lists[7],
        "F Frequency" : frequency_lists[8], "N Frequency" : frequency_lists[9], "G Frequency" : frequency_lists[10], "H Frequency" : frequency_lists[11], "L Frequency" : frequency_lists[12], "R Frequency" : frequency_lists[13], "W Frequency" : frequency_lists[14],
        "A Frequency" : frequency_lists[15], "V Frequency" : frequency_lists[16], "E Frequency" : frequency_lists[17], "Y Frequency" : frequency_lists[18], "M Frequency" : frequency_lists[19], "X Frequency" : frequency_lists[20], "Z Frequency" : frequency_lists[21], "Length" : seq_len, "Labels" : label_list}

freq_df = pd.DataFrame(data, columns = ["C Frequency", "D Frequency", "S Frequency", "Q Frequency", "K Frequency", "I Frequency", "P Frequency", "T Frequency", "F Frequency", "N Frequency", "G Frequency", "H Frequency", "L Frequency", "R Frequency", "W Frequency", "A Frequency", "V Frequency", "E Frequency", "Y Frequency", "M Frequency", "X Frequency", "Z Frequency", "Length", "Labels"])




# Trying the Boxplot with the different frequencies
# Draw a nested boxplot to show bills by day and time

ax = sns.violinplot(x="Labels", y="C Frequency", data=freq_df)
plt.savefig('C.png', dpi = 100)
plt.show()


ax = sns.violinplot(x="Labels", y="D Frequency", data=freq_df)
plt.savefig('D.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="S Frequency", data=freq_df)
plt.savefig('S.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="Q Frequency", data=freq_df)
plt.show()

ax = sns.violinplot(x="Labels", y="K Frequency", data=freq_df)
plt.savefig('K.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="I Frequency", data=freq_df)
plt.savefig('I.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="P Frequency", data=freq_df)
plt.savefig('P.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="T Frequency", data=freq_df)
plt.savefig('T.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="F Frequency", data=freq_df)
plt.savefig('F.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="N Frequency", data=freq_df)
plt.savefig('N.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="G Frequency", data=freq_df)
plt.savefig('G.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="H Frequency", data=freq_df)
plt.savefig('H.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="L Frequency", data=freq_df)
plt.savefig('L.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="R Frequency", data=freq_df)
plt.savefig('R.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="W Frequency", data=freq_df)
plt.savefig('W.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="A Frequency", data=freq_df)
plt.savefig('A.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="V Frequency", data=freq_df)
plt.savefig('V.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="E Frequency", data=freq_df)
plt.savefig('E.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="Y Frequency", data=freq_df)
plt.savefig('Y.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="M Frequency", data=freq_df)
plt.savefig('M.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="X Frequency", data=freq_df)
plt.savefig('X.png', dpi = 100)
plt.show()

ax = sns.violinplot(x="Labels", y="Z Frequency", data=freq_df)
plt.savefig('Z.png', dpi = 100)
plt.show()

# TODO: consider group box plot as well


#%% # Plotting the average distance between amino acids

for element in prot_series:
    element


prot_labels = {"C_avg_d", "D_avg_d", "S_avg_d", "Q_avg_d", "K_avg_d", "I_avg_d", "P_avg_d", "T_avg_d", "F_avg_d", "N_avg_d", "G_avg_d", "H_avg_d", "L_avg_d", "R_avg_d", "W_avg_d", "A_avg_d", "V_avg_d", "E_avg_d", "Y_avg_d", "M_avg_d","X_avg_d", "Z_avg_d"}
prot_list = {"C", "D", "S", "Q", "K", "I", "P", "T", "F", "N", "G", "H", "L", "R", "W", "A", "V", "E", "Y", "M","X", "Z"}

avg_dist_df = pd.DataFrame(columns=prot_labels)
avg_dist_df

part_to_add = "_avg_d"

first = 0

for i in range(len(prot_df)):
    seq = prot_df.iat[i, 2]
    prot_dict = {"C_avg_d": 0, "D_avg_d": 0, "S_avg_d": 0, "Q_avg_d": 0, "K_avg_d": 0, "I_avg_d": 0, "P_avg_d": 0,
                 "T_avg_d": 0, "F_avg_d": 0, "N_avg_d": 0, "G_avg_d": 0, "H_avg_d": 0, "L_avg_d": 0, "R_avg_d": 0,
                 "W_avg_d": 0, "A_avg_d": 0, "V_avg_d": 0, "E_avg_d": 0, "Y_avg_d": 0, "M_avg_d": 0, "X_avg_d": 0,
                 "Z_avg_d": 0}
    for prot in prot_list:
        avg_dist = 0
        freq_count = 0
        dist_count = 0


        for char in seq:
            if char == prot:
                freq_count += 1

            if freq_count >= 1:
                dist_count += 1
        if freq_count >= 1:
            avg_dist = freq_count / dist_count
            prot_add = prot + part_to_add
            prot_dict[prot_add] = avg_dist
            if first == 0:
                print(prot_dict)
                first = 1
        avg_dist_df = avg_dist_df.append(prot_dict, ignore_index=True)

        # TODO fix the wrong value error, first sequence definitely contains a K


print(prot_df.iat[0, 2])
