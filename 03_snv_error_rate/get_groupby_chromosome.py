import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import collections
import numpy as np
import re
import time

def prob_to_phred(prob: float):
    """
    Convert probability of base call being correct into phred score
    Values are clipped to stay within 0 to 60 phred range
    Args:
        prob  (float): probability of base call being correct
    Returns:
        phred_score (byte)
    """
    return np.rint(-10 * np.log10(np.clip(1-prob, 1-0.999999, 0.999999))).astype('B')




def prob_to_phred(prob: float):
    """
    Convert probability of base call being correct into phred score
    Values are clipped to stay within 0 to 60 phred range
    Args:
        prob  (float): probability of base call being correct
    Returns:
        phred_score (byte)
    """
    return np.rint(-10 * np.log10(np.clip(1-prob, 1-0.999999, 0.999999))).astype('B')



def custom_sort(df):
    # Extract relevant sorting information
    method = df['method'].iloc[0]
    name = df['name'].iloc[0]

    # Create sorting keys
    if method == 'NANO':
        order = 1
    elif method == 'NOVA':
        order = 2
    elif method == 'NOVA-ecco':
        order = 3
    elif method == 'RCS':
        order = 4
    else:
        order = 5  # in case there's an unexpected method

    # Extract the sample number (HCxx)
    sample_number = int(re.search(r'HC(\d+)', name).group(1))

    return (order, sample_number)


def synonym_name(x):
    if x == "RCS":
        return "Consensus\nNanoRCS"
    elif x == "NOVA":
        return "NovaSeq"
    elif x == "NOVA-ecco":
        return "NovaSeq\npaired\nend"

    elif x == "NANO":
        return "Raw NanoRCS"
    else:
        raise "method is not NOVA, NANO or CYC."

## New files: 

names = []
methods = []
unsort_dfs = []



# READ RAW NANOPORE SEQUENCING
in_path = "/Users/liting/00_projects/NanoRCS/output/processed_data/Fig2A_pickles/"
#in_path = "/Users/liting/00_projects/genome_wide_cyclomics_project/Figure1_Qscore/data/02_exclude_overlap_pickle/NANO"



for a_file in os.listdir(in_path):
    if not a_file.endswith("pickle.gz"):
        continue

    if "NANO" in a_file:
        print("Nanopore sequencing", a_file)

        name = a_file.split("_overlap")[0]
        sample = name.split("_")[0].split("-")[0]
        df = pd.read_pickle(os.path.join(in_path, a_file))
        df['method'] = "NANO"
        df['name'] = sample
        names.append(sample)
        methods.append("NANO")
        unsort_dfs.append(df)

## READ ILLUMINA NOVASEQ without error correction
#in_path = "/Users/liting/01_data/MANUSCRIPT_DATA/Figure1_Errorrate/output3_select_pickle_for_plotting"


for a_file in os.listdir(in_path):
    if not a_file.endswith("pickle.gz"):
        continue
    if "ecco_NOVA" in a_file:
        continue
    if "NOVA" in a_file:


        name = a_file.split("_overlap")[0]

        sample = name.split("_")[0]
        print('NOVA sample', sample)
        df = pd.read_pickle(os.path.join(in_path, a_file))
        df['method'] = "NOVA"
        df['name'] = sample

        names.append(sample)
        methods.append("NOVA")
        unsort_dfs.append(df)


## READ ILLUMINA NOVASEQ ERROR CORRECTION



for a_file in os.listdir(in_path):
    if not a_file.endswith("pickle.gz"):
        continue
    if "ecco_NOVA" in a_file:
        print("NovaSeq Error correction", a_file)

        name = a_file.split("_overlap")[0]
        sample = name.split("_")[0].split("-")[0]
        df = pd.read_pickle(os.path.join(in_path, a_file))
        df['method'] = "NOVA-ecco"
        df['name'] = sample

        names.append(sample)
        methods.append("NOVA-ecco")
        unsort_dfs.append(df)

        
## READ CYCLOMICSSEQ 

for a_file in os.listdir(in_path):
    if not a_file.endswith("pickle.gz"):
        continue

    if "RCS" in a_file:
        print("NanoRCS", a_file)
        name = a_file.split("_overlap")[0]
        method = name.split("_")[1]
        sample = name.split("_")[0]
#         print(name, method,sample)
        df = pd.read_pickle(os.path.join(in_path, a_file))
        df['method'] = "RCS"
        df['name'] = sample
        names.append(sample)
        methods.append("RCS")
        unsort_dfs.append(df)

        

for df in unsort_dfs:
    if method == "NANO":
        df['Sequencing Method'] = "Raw NanoRCS"
    elif method == "NOVA":
        df['Sequencing Method'] = "NovaSeq"
    elif method == "NOVA-ecco":
        df['Sequencing Method'] = "NovaSeq paired end"
    elif method == "CYC":
        df['Sequencing Method'] = 'Consensus NanoRCS'
        

sorted_indexes = sorted(range(len(unsort_dfs)), key=lambda i: custom_sort(unsort_dfs[i]))

sorted_names = [names[i] for i in sorted_indexes]
sorted_methods = [methods[i] for i in sorted_indexes]
dfs = [unsort_dfs[i] for i in sorted_indexes]



# Long plotting names

        
source_data_out = []
for name, method,  df in zip(names, methods, dfs):
    df_sum_by_chrom = df.groupby(['chromosome']).sum()
    df_sum_by_chrom['mismatch_rate'] = df_sum_by_chrom.apply(lambda x:  x['mismatches'] / x['inferred_read_length'], axis = 1 )
    df_sum_by_chrom['chromosome'] = df_sum_by_chrom.index
    construct_sorted_values = []
    for x in range(1,23):
        mismatch_rate = df_sum_by_chrom[df_sum_by_chrom['chromosome'] == str(x)]['mismatch_rate'].values[0]
        construct_sorted_values.append(mismatch_rate)
        source_data_out.append([name, method, synonym_name(method), x, mismatch_rate])

        
# Read in source data
out = pd.DataFrame(source_data_out, columns = ['sample','method','method_plot_name','chromosome','error_rate_in_1M_reads'])
out.to_csv(f"../source_data/{filename}.csv")

