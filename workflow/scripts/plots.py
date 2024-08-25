import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3
import warnings

def npg_palette():
    palette = ["#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
               "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF"]
    return sns.color_palette(palette, len(palette))

def plot_venn(sig_dict, tools, metrics, ax = None, pretty_print=None):

    if not isinstance(tools, list): 
        tools = [tools]
    if not isinstance(metrics, list): 
        metrics = [metrics]
        
    n_sets = len(tools) * len(metrics)
    if n_sets not in [2,3]:
        raise Exception(f"Venn diagram with {n_sets} sets not supported") # TO DO: venn4
    
    #idx = pd.IndexSlice
    #summary_df.loc[:, idx[tools, metrics, "qvalue"]]

    toolIsTopLevel = len(tools) < len(metrics)

    sets, labels = [], []
    for tool in tools:
        for metric in metrics:
            s = sig_dict[tool][metric]
            sets.append(s)
            labels.append(metric if toolIsTopLevel else tool)
            if pretty_print: # pretty print labels
                if labels[-1] in pretty_print:
                    labels[-1] = pretty_print[labels[-1]]

    if not ax:
        fig, ax = plt.subplots(1,1, figsize=(4,4))

    with warnings.catch_warnings(action="ignore"):
        if len(sets) == 2:
            venn2(sets,set_labels=labels, ax=ax)
        elif len(sets) == 3:
            venn3(sets,set_labels=labels, ax=ax)
    
    plt.title(tool if toolIsTopLevel else metric)
    
    if not ax:
        return fig