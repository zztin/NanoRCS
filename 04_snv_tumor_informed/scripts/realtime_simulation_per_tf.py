import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np


def get_average_read_length(read_length_table,
                            max_read_length = 1000):
    read_lengths = pd.read_csv(read_length_table, header=1)
    nanorcs_read_length = read_lengths[read_lengths['Technique'] == 'NanoRCS']
    max_read_length = max_read_length - 31
    nanorcs_average_read_length_density = nanorcs_read_length.sum(axis=0, numeric_only=True)[:max_read_length]
    return nanorcs_average_read_length_density

def get_throughput(a_folder, ):

    ## 1. Get column "Basecalled Reads Passed" in all csv sequence throughput files in a folder.
    # Collect dfs (TODO: trim down to only valid columns if necessary)
    dfs = []
    run_names = []
    for a_file in os.listdir(a_folder):
        if a_file.endswith('.csv'):
            print(a_file)
            df1 = pd.read_csv(a_folder + a_file)
            run_name = a_file.split('/')[-1].split('.')[0]
            dfs.append(df1)
            run_names.append(run_name)
    # Collect useful columns
    for i, df_x in enumerate(dfs):
        if i == 0:
            df_final = df_x
        else:
            df_final = df_final.join(df_x.loc[:, df1.columns.str.contains('Basecalled Reads Passed')], how='outer',
                                     rsuffix=run_names[i])

    mean_read_throughput_per_mins = df_final.loc[:, df_final.columns.str.contains('Basecalled Reads Passed')].mean(
        axis=1)[::EVERY_X_MINS]
    mean_read_throughput_per_mins_std = df_final.loc[:, df_final.columns.str.contains('Basecalled Reads Passed')].std(
        axis=1)[::EVERY_X_MINS]
    mean_read_throughput_per_mins_cv = mean_read_throughput_per_mins_std / mean_read_throughput_per_mins

    return mean_read_throughput_per_mins, mean_read_throughput_per_mins_cv


def get_total_cfDNA_bases(nanorcs_average_read_length_density,
                          mean_read_throughput_per_mins,
                          mean_read_throughput_per_mins_cv,
                          max_read_length = 1000):

    values = np.arange(0, max_read_length-31)  # from 0 to X inclusive
    probabilities = nanorcs_average_read_length_density / np.sum(nanorcs_average_read_length_density)
    # Sampling based on known distribution in all samples (exported from Suppl Table 4)
    nested_list_for_summary_df = []
    for j, (num_samples, cv) in enumerate(
            zip(mean_read_throughput_per_mins.values, mean_read_throughput_per_mins_cv.values)):
        sampled_values = np.random.choice(values, size=int(num_samples), p=probabilities)
        total_covered_consensus_bases_per_x_mins = sampled_values.sum()
        # Print the sampled values
        # print("Consensus reads:", num_samples, "Consensus bases:",total_covered_consensus_bases_per_x_mins, "Average read length:", total_covered_consensus_bases_per_x_mins/num_samples )
        nested_list_for_summary_df.append([j * EVERY_X_MINS, num_samples, total_covered_consensus_bases_per_x_mins, cv])

    df_throughput = pd.DataFrame(nested_list_for_summary_df,
                                 columns=['Experiment Time (minutes)', 'Basecalled Reads Passed',
                                          'Total cfDNA bases (NanoRCS consensus)', 'Coefficient of variation']
                                 )
    return df_throughput


if __name__ == '__main__':
    ## Variables
    amount_of_tumor_specific_snv = 5500
    TFs = [0.001, 0.005, 0.01, 0.1, 0.2, 0.4, 0.6]

    ## FIXED PARAMETERS
    genome_size = 3_117_275_501
    coverage = [0.04, 0.25, 0.8]
    error_rates = [0.00744, 0.00082]
    nanorcs_error_rate = 0.00082
    tumor_fraction_min = 0.0001
    tumor_fraction_max = 1
    VAF = 0.5
    EVERY_X_MINS = 5
    max_read_length = 1000
    ## PATHS
    a_folder = '/path/to/sequencing_summary_files'
    read_length_table = '/path/to/suppl_table_4.csv'

    ## Get read length distribution in cfDNA NanoRCS reads
    nanorcs_average_read_length_density = get_average_read_length(read_length_table,
                                                                  max_read_length = max_read_length)
    plt.plot(nanorcs_average_read_length_density)
    plt.show()


    ## Get throughput from sequencing summary files
    mean_read_throughput_per_mins, mean_read_throughput_per_mins_cv = get_throughput(a_folder)

    ## Get Total cfDNA base counts by combining throughput and read length
    df_throughput = get_total_cfDNA_bases(nanorcs_average_read_length_density,
                          mean_read_throughput_per_mins,
                          mean_read_throughput_per_mins_cv,
                          max_read_length = max_read_length)


    ## 3.Plot a curve of throughput. (TODO: add CV)
    sns.scatterplot(x='Experiment Time (minutes)',y= 'Total cfDNA bases (NanoRCS consensus)', data=df_throughput,
                    label='PromethION NanoRCS consensus bases output')
    plt.legend()
    plt.show()

    ## 4. Estimate when will discovery of tumor specific SNVs is higher than background
    ## in different TFs & tumor SNV count?
    ## TODO: Binomial trial instead of constant rate.
    for TF in TFs:
        df_throughput[f'Tumor SNV detected at {TF}'] = np.floor(
            df_throughput['Total cfDNA bases (NanoRCS consensus)'] * (
                        amount_of_tumor_specific_snv / genome_size) * TF * VAF)
    df_throughput[f'Tumor SNV detected at BACKGROUND'] = np.floor(
        df_throughput['Total cfDNA bases (NanoRCS consensus)'] * (
                    amount_of_tumor_specific_snv / genome_size) * nanorcs_error_rate)

    ## 5. Plot discovery
    for TF in TFs:
        sns.scatterplot(x='Experiment Time (minutes)', y=f'Tumor SNV detected at {TF}', data=df_throughput,
                        label=f'TF = {TF}')
    sns.scatterplot(x='Experiment Time (minutes)', y=f'Tumor SNV detected at BACKGROUND', data=df_throughput,
                    label=f'Background error rate')
    plt.title(f'Tumor SNV = {amount_of_tumor_specific_snv}')
    plt.legend()
    plt.ylabel("Tumor SNV detected (n)")
    plt.show()


    ## Plot - Zoom in to specific time period
    time_tuple = 550, 1500

    for TF in TFs:
        sns.scatterplot(x='Experiment Time (minutes)', y=f'Tumor SNV detected at {TF}', data=df_throughput,
                        label=f'TF = {TF}')
    sns.scatterplot(x='Experiment Time (minutes)', y=f'Tumor SNV detected at BACKGROUND', data=df_throughput,
                    label=f'Background error rate')
    plt.title(f'Tumor SNV = {amount_of_tumor_specific_snv}')
    plt.legend(bbox_to_anchor=(1, 1))
    plt.ylabel("Tumor SNV detected (n)")
    plt.ylim(-3, 20)
    plt.xlim(time_tuple)
    # plt.xticks(range(0, 250, 10))
    plt.show()



