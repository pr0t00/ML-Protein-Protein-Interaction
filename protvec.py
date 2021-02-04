import pandas as pd


# read protvec file
protvec_file = "protVec_100d_3grams.csv"
protvec_data = pd.read_csv(protvec_file, sep="\t", header=0)

# dictionary to hold all 3-grams
protvec = {}


# number of rows
rows_protvec = protvec_data.shape[0]

# tree = {"hi" : [1,3]}
# tree["hi"][0]


for row in range(0,rows_protvec):

    dict_key = protvec_data.iloc[row][0] # [row][column]

    dict_values = []
    for column in range(1,101):
        value = protvec_data.iloc[0][column]
        dict_values.append(value)
    protvec[dict_key] = dict_values
    if row % 500 == 0:
        print(str(row) + " / 9048")


# TODO: check whether dictionary is created in the correct way and how I can retrieve specific list positions


# test sequence
seq = "MSRFSIEGKSLKLDAITTEDEKSVFAVLLEDDSVKEIVLSGNTIGTEAARWLSENIASKKDLEIAEFSDIFTGRVKDEIPEALRLLLQALLKCPKLHTVRLSDNAFGPTAQEPLIDFLSKHTPLEHLYLHNNGLGPQAGAKIARALQELAVNKKAKNAPPLRSIICGRNRLENGSMKEWAKTFQSHRLLHTVKMVQNGIRPEGIEHLLLEGLAYCQELKVLDLQDNTFTHLGSSALAIALKSWPNLRELGLNDCLLSARGAAAVVDAFSKLENIGLQTLRLQYNEIELDAVRTLKTVIDEKMPDLLFLELNGNRFSEEDDVVDEIREVFSTRGRGELDELDDMEELTDEEEEDEEEEAESQSPEPETSEEEKEDKELADELSKAHI"

# creating a list which is then filled with 100 zero values
sequence_vec = []
for i in range(0,100):
    sequence_vec.append(0)


# iterate overall overlapping 3-grams (e.g. MAK, AKL, KLP, LPQ,...)
for i in range(0, len(seq)-2, 1):
    threegram = seq[i:i+3]
    print(seq[i:i+3])
    #prot_vec_3_gram = protvec[threegram]
    sequence_vec = sequence_vec #+ prot_vec_3_gram

print(sequence_vec)




# playground

#list = [[[1,2],2],[3,4]]
#dic = {"three" : [1,2]}
#dic["three"]
#print(list[[0]])