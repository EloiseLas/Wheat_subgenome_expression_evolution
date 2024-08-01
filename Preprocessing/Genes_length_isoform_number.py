from Bio import SeqIO
import pandas as pd
import numpy as np

def average(Length_list):
    # if np.std(Length_list)<100:
    return np.mean(Length_list)
    # else:
    #     return "NA"

def std_flag(Length_list):
    return np.std(Length_list)>100


fasta_sequences = SeqIO.parse(open("Data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_cds_from_genomic.fna"),'fasta')
ids=[]
lengths=[]

for fasta in fasta_sequences:
    name, sequence = fasta.description, str(fasta.seq)
    ID = name.split("[gene=",1)[1].split("]",1)[0]
    seq_len = len(sequence)
    ids.append(ID)
    lengths.append(seq_len)

df = pd.DataFrame({'ID': ids, 'Length': lengths})

df["Isoform_number"]=df.groupby(["ID"])["Length"].transform("count")

# df = df.drop_duplicates()
# print(df)

ids = df["ID"]
dupli=df[df.duplicated(subset=['ID'],keep=False)]
dupli["Length_list"]=dupli.groupby("ID")["Length"].transform(lambda x : [x.tolist()]*len(x))
dupli["Length"]=dupli.apply(lambda row : average(row["Length_list"]), axis = 1)
dupli["High_std"]=dupli.apply(lambda row : std_flag(row["Length_list"]), axis = 1)
dupli=dupli.drop(["Length_list"], axis = 1)
dupli = dupli.drop_duplicates()
print(dupli)


notdupli=df[~df.duplicated(subset=['ID'],keep=False)]
# print(notdupli)

data = pd.concat([dupli,notdupli],axis=0)
data["High_std"] = data["High_std"].fillna(False)
data_length=data[["ID","Length","High_std"]]
data_iso=data[["ID","Isoform_number"]]

print(data_length)

counts=data_length.groupby(["High_std"]).count()
print(counts)

data_length.to_csv('Data/Genes_length.csv')
data_iso.to_csv('Data/Genes_iso.csv')