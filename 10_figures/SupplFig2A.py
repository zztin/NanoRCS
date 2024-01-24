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
    # Get long plotting names
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
            raise "method is not NOVA, NANO or RCS."

    # Plotting parameters:
    color_type_dict = {"RCS":plt.cm.get_cmap('Greens'),
                       "NOVA":plt.cm.get_cmap('Greys'),
                       "NOVA-ecco":plt.cm.get_cmap('bone'),
                       "NANO": plt.cm.get_cmap('Blues')}
    marker_type_dict =  {"HC01":"o","HC02":"s", "HC03":"^"}

    # Load data
    df2 = pd.read_csv("./source_data/SupplFig2A.csv", index_col=[0])
    names = list(df2['sample'].unique())*4
    methods = list(np.repeat(list(df2['method'].unique()),3))
    # Plot
    plt.subplots(figsize=(7.93, 3.12))
    i = 0
    for name, method in zip(names, methods):
        i += 1
        error_rate_by_chrom = df2[(df2['sample'] == name) & (df2['method'] == method)].sort_values(by='chromosome')['error_rate_in_1M_reads'].values
        marker_type = marker_type_dict[name]
        color_type = color_type_dict[method]
        plt.plot([str(i) for i in range(1,23)],
                 error_rate_by_chrom,
                 marker=marker_type,
                 linestyle="None",
                 label= f"{name}_{method}",
                 color = color_type((i+3)*10)
                )
    plt.legend(bbox_to_anchor= (1,1))
    plt.ylabel("SNV Error rate")
    plt.xlabel("Chromosome")
    plt.tight_layout()
    plt.savefig("../output/output_figures/SupplFig2A.pdf")
    plt.show()
