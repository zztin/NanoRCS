# Import relevant packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set font size, style
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['font.family'] = ['Arial']
SMALL_SIZE = 5
MEDIUM_SIZE = 6
BIGGER_SIZE = 7
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
sns.set(rc={'font.family': 'Arial'}, font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid': False})

if __name__ == '__main__':
    # Load data base
    df3 = pd.read_csv("./source_data/SupplFig2A_2.csv", index_col = 0)
    # Plotting
    fig, ax = plt.subplots(figsize = (6,6))
    sns.boxplot(x = 'q_score_threshold',
                y = 'fraction_reads_above_threshold',
                hue = 'method',
                data = df3,
                palette={"NANO": "#d8e6f6", "NOVA": "#d4d4d4", "NOVA-ecco":"#606680","RCS":"#60ba6c"},
                ax = ax)
    plt.ylabel("Fraction of reads above SNV error rate")
    plt.xlabel("SNV Error rate (Phred score, -log10(E)*10)")
    plt.savefig("../output/output_figures/SupplFig2A_2.pdf")
    plt.show()
