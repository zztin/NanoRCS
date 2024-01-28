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
    # Import dataframe
    df = pd.read_csv("./source_data/Fig2A.csv", index_col=[0,1])
    # Plotting snv error rate per sequencing and analysis method
    p1, ax = plt.subplots(figsize=(5*1.1, 4*1.1))
    p1 = sns.barplot(data=df,
                     orient="h",
                     x='SNV error rate',
                     y="Sequencing Method",
                     order=["Raw NanoRCS", "Consensus\nNanoRCS", "NovaSeq", "NovaSeq\npaired\nend", ],
                     palette=['lightgrey', 'grey', 'lightgrey', 'grey', ],
                     )
    ax.set_ylabel('')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks(np.arange(0.000, 0.009, 0.001))
    ax.tick_params(axis='x', which='major', length=4, grid_color='grey')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("../output/output_figures/Fig2A.pdf")
    plt.show()