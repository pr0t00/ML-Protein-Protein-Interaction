def parse_fasta_to_df(file_path, content_alias, id_alias = "ProtID"):
    from Bio import SeqIO
    import pandas as pd
    with open(file_path) as fasta_file:
        ids = []
        contents = []
        for record in SeqIO.parse(fasta_file, 'fasta'):
            ids.append(record.id)
            contents.append(record.seq)
        df = pd.DataFrame(data = {id_alias: ids, content_alias: contents})
        return df.set_index(id_alias)
    
def filter_unique_rownames(df):
    return df[~df.index.duplicated(keep='first')]


def import_protvec(filepath, namescol = "words"):
    """
    Import data frame of ProtVec 3-grams. 
    
    :param filepath: path to a TSV.
    :param namescol: name of a column with row names 
    :return: pandas dataframe with 3-grams as rownames.
    """
    import pandas as pd
    protvec_df = pd.read_csv(filepath, sep = "\t", header = 0)
    protvec_df_3gramidx = protvec_df.set_index(namescol)
    return(protvec_df_3gramidx)

def get3gramvec(threegr_df, threegr_name, as_list = False):
    if not threegr_name in threegr_df.index:
        raise ValueError(''.join(["The supplied ProtVec dataset is not trained for the threegram: ", threegr_name]))
    vec = threegr_df.loc[threegr_name].values
    if (as_list):
        vec = vec.tolist()
    return vec
    
def convert_seq_to_protvec(seq, threegr_df, substitute_any_with="G"):
    """
    Get ProtVec representation of a given sequence
    """
    import numpy as np
    protvec = np.zeros(100)
    for i in range(0, len(seq) - 3):
        this3gram = str(seq[i:i+3])
        this3gram = this3gram.replace("X", substitute_any_with)
        if not this3gram in threegr_df.index:
            # skip untrained 3grams
            continue
        this3gramvec = get3gramvec(threegr_df, this3gram)
        protvec = np.add(protvec, this3gramvec)
    return protvec
