import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import collections
import numpy as np

# Phred score definition
def prob_to_phred(prob):
    """
    Convert probability of base call being correct into phred score
    Values are clipped to stay within 0 to 60 phred range
    Args:
        prob  (float): probability of base call being correct
    Returns:
        phred_score (byte)
    """
    return np.rint(-10 * np.log10(np.clip(1-prob, 1-0.999999, 0.999999))).astype('B')




# Long plotting names
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


def synonym(x):
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



if __name__ == '__main__':
    # Collect data from pickle files
    dfs_dict = collections.defaultdict(dict)
    in_path = "../output/processed_data/Fig2A_pickles"
    # Read raw NanoRCS sequencing
    for a_file in os.listdir(in_path):
        if not a_file.endswith("pickle.gz"):
            continue
        if "NANO" in a_file:
#            print("Raw NanoRCS file collected:", a_file)
            name = a_file.split("_overlap")[0]
            sample = name.split("_")[0].split("-")[0]
            df = pd.read_pickle(os.path.join(in_path, a_file))
            df['method'] = "NANO"
            df['Sequencing Method'] = "Raw NanoRCS"
            df['name'] = sample
            dfs_dict["NANO"][sample] = df
    ## READ ILLUMINA NOVASEQ without error correction
#    in_path = "/Users/liting/01_data/MANUSCRIPT_DATA/Figure1_Errorrate/output3_select_pickle_for_plotting"
    for a_file in os.listdir(in_path):
        if not a_file.endswith("pickle.gz"):
            continue
        if "ecco_NOVA" in a_file:
            continue
        if "NOVA" in a_file:

            name = a_file.split("_overlap")[0]
            sample = name.split("_")[0]
            # print("NovaSeq file collected:", a_file)
            df = pd.read_pickle(os.path.join(in_path, a_file))
            df['method'] = "NOVA"
            df['Sequencing Method'] = "NovaSeq"
            df['name'] = sample
            dfs_dict["NOVA"][sample] = df
    ## READ ILLUMINA NOVASEQ ERROR CORRECTION
    for a_file in os.listdir(in_path):
        if not a_file.endswith("pickle.gz"):
            continue
        if "ecco_NOVA" in a_file:
#            print("NovaSeq Error correction file collected: ", a_file)
            name = a_file.split("_overlap")[0]
            sample = name.split("_")[0].split("-")[0]
            df = pd.read_pickle(os.path.join(in_path, a_file))
            df['method'] = "NOVA-ecco"
            df['Sequencing Method'] = "NovaSeq paired end"
            df['name'] = sample
            dfs_dict["NOVA-ecco"][sample] = df
    ## READ NanoRCS
    for a_file in os.listdir(in_path):
        if not a_file.endswith("pickle.gz"):
            continue
        if "RCS" in a_file:
 #           print("Consensus NanoRCS file collected:", a_file)
            name = a_file.split("_overlap")[0]
            method = name.split("_")[1]
            sample = name.split("_")[0]
            df = pd.read_pickle(os.path.join(in_path, a_file))
            df['method'] = "RCS"
            df['Sequencing Method'] = 'Consensus NanoRCS'
            df['name'] = sample
            dfs_dict["RCS"][sample] = df

    # Done collecting data from pickle files, start generating combined statistics.
    ## Fig2A
    df_concat = []
    for method in ['NANO','NOVA','NOVA-ecco','RCS']:
        for name in ['HC01','HC02','HC03']:
            df = dfs_dict[method][name]
            if df.shape[0] == 0:
                print(name, method, "Skipped")
                continue
            else:
                # python 3.8
                df_sum_by_chrom = df.groupby(['batch'],).sum( numeric_only=True)
                df_sum_by_chrom['mismatch_rate'] = df_sum_by_chrom.apply(
                    lambda x: x['mismatches'] / x['inferred_read_length'], axis=1)
                df_sum_by_chrom['Quality Score (phred)'] = df_sum_by_chrom['mismatch_rate'].apply(
                    lambda x: prob_to_phred(1 - x))
                df_sum_by_chrom['Method'] = method
                df_sum_by_chrom['sample'] = name
            df_concat.append(df_sum_by_chrom)
    df_concated = pd.concat(df_concat)
    df_concated['Sequencing Method'] = df_concated['Method'].apply(synonym)
    df_concated['SNV Error rate'] = df_concated['mismatch_rate'].copy()
    df_concated = df_concated.reset_index()
    df_concated = df_concated.reset_index(drop=True)
    Fig2A = df_concated.groupby(by=['Method', 'sample']).sum(numeric_only=True)[['inferred_read_length', 'nm_tag', 'indels', 'mismatches']]
    Fig2A['SNV error rate'] = Fig2A.apply(lambda x: x['mismatches'] / x['inferred_read_length'], axis=1)
    Fig2A['Method'] = [x for x, y in Fig2A.index]
    Fig2A['Sequencing Method'] = Fig2A['Method'].apply(synonym)
    Fig2A.to_csv("../10_figures/source_data/Fig2A.csv")

    ## Main text "per method error rate":
    for method_print in df_concated['Sequencing Method'].unique():
        mean_rate = df_concated[df_concated['Sequencing Method'] == method_print]['mismatch_rate'].mean()
        print(method_print.replace('/n',' '))
        print("Error rate:", round(mean_rate, 9))
        print("Phred score (Q score):", prob_to_phred(1 - mean_rate))
        print("One in how many SNV is an error:", round(1 / mean_rate))

    ## SupplFig2A
    source_data_out = []
    for method in ['NANO','NOVA','NOVA-ecco','RCS']:
        for name in ['HC01','HC02','HC03']:
            df = dfs_dict[method][name]
            if df.shape[0] == 0:
                print(name, method, "Skipped")
                continue
            else:
                df_sum_by_chrom = df.groupby(['chromosome']).sum(numeric_only=True)
                df_sum_by_chrom['mismatch_rate'] = df_sum_by_chrom.apply(
                    lambda x: x['mismatches'] / x['inferred_read_length'], axis=1)
                df_sum_by_chrom['chromosome'] = df_sum_by_chrom.index
                construct_sorted_values = []
                for x in range(1, 23):
                    mismatch_rate = df_sum_by_chrom[df_sum_by_chrom['chromosome'] == str(x)]['mismatch_rate'].values[0]
                    source_data_out.append([name, method, synonym_name(method), x, mismatch_rate])
    out = pd.DataFrame(source_data_out,
                       columns=['sample', 'method', 'method_plot_name', 'chromosome', 'error_rate_in_1M_reads'])
    out.to_csv("../10_figures/source_data/SupplFig2A.csv")


    ## SupplFig2B
    fraction_q_score = []
    for method in ['NANO','NOVA','NOVA-ecco','RCS']:
        for name in ['HC01','HC02','HC03']:
            df = dfs_dict[method][name]
            if df.shape[0] == 0:
                print(name, method, "Skipped")
                continue
            else:
                df['Q_score_per_read'] = df.apply(
                    lambda row: 10 * -np.log10(row['mismatches'] / row['inferred_read_length']) if row['mismatches'] > 0 \
                        else np.inf , axis = 1)
                for q_threshold in [10, 20, 30, 40, ]:
                    fraction_above =  df[df['Q_score_per_read'] >= q_threshold].shape[0] / df.shape[0]
                    avg_read_length = df[df['Q_score_per_read'] >= q_threshold]['inferred_read_length'].mean()
                    amount_read = df[df['Q_score_per_read'] >= q_threshold].shape[0]
                    fraction_q_score.append(
                        [method, name, q_threshold, fraction_above, avg_read_length, amount_read])
    fraction_above_q_df = pd.DataFrame(fraction_q_score, columns=["method",
                                                                  "sample",
                                                                  "q_score_threshold",
                                                                  "fraction_reads_above_threshold",
                                                                  "avg_read_length",
                                                                  "amount_read"
                                                                  ]
                       )
    fraction_above_q_df.to_csv("../10_figures/source_data/SupplFig2B.csv")

