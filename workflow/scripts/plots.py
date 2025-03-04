import os
import sys
import warnings
import yaml
import pickle
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import seaborn as sns

from upsetplot import from_memberships
from upsetplot import UpSet

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from scripts.explore_results import get_sig_dict


def npg_palette():
    palette = [
        "#E64B35FF",
        "#4DBBD5FF",
        "#00A087FF",
        "#3C5488FF",
        "#F39B7FFF",
        "#8491B4FF",
        "#91D1C2FF",
        "#DC0000FF",
        "#7E6148FF",
        "#B09C85FF",
    ]
    return sns.color_palette(palette, len(palette))


def plot_venn(sig_dict, tools, metrics, ax=None, pretty_print=None):
    if not isinstance(tools, list):
        tools = [tools]
    if not isinstance(metrics, list):
        metrics = [metrics]

    n_sets = len(tools) * len(metrics)
    if n_sets not in [2, 3]:
        raise Exception(f"Venn diagram with {n_sets} sets not supported")  # TO DO: venn4

    # idx = pd.IndexSlice
    # summary_df.loc[:, idx[tools, metrics, "qvalue"]]

    toolIsTopLevel = len(tools) < len(metrics)

    sets, labels = [], []
    for tool in tools:
        for metric in metrics:
            s = sig_dict[tool][metric]
            sets.append(s)
            labels.append(metric if toolIsTopLevel else tool)
            if pretty_print:  # pretty print labels
                if labels[-1] in pretty_print:
                    labels[-1] = pretty_print[labels[-1]]

    if not ax:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    with warnings.catch_warnings(action="ignore"):
        if len(sets) == 2:
            venn2(sets, set_labels=labels, ax=ax)
        elif len(sets) == 3:
            venn3(sets, set_labels=labels, ax=ax)

    plt.title(tool if toolIsTopLevel else metric)

    if not ax:
        return fig


def make_bar_plots(
    summary_dict: Dict,
    figpath: str,
    project_name: str,
    lib_names: Dict,
    pretty_print: Dict,
    qval: float = 0.05,
    palette=npg_palette(),
    max_depth: int = 0,
    ext: str = "pdf",
):
    sns.set_theme(font_scale=1.2)

    nlib = len(lib_names.keys())

    with sns.axes_style("whitegrid"):
        fig, axes = plt.subplots(2, nlib, figsize=(nlib * 5, 10))
    axes = axes.flatten()

    ### Intersection depth

    for ax, lib in zip(axes[:nlib], lib_names.keys()):
        depth_df = summary_dict[lib]["depth_df"]
        if len(depth_df) < 1:
            print(f"No terms found for {lib}")
            continue
        h = sns.histplot(
            depth_df["Depth"],
            bins=depth_df["Depth"].max() - depth_df["Depth"].min() + 1,
            discrete=True,
            ax=ax,
            alpha=1,
            color=palette[3],
        )
        h.bar_label(h.containers[0])
        ax.set(title=lib, xlabel="Robustness", ylabel="Terms")
        xmax = max_depth if max_depth else depth_df["Depth"].max()
        ax.set_xticks(range(1, xmax + 1))
        ax.set_xlim(0.25, xmax + 0.75)
        ax.grid(axis="x")

    ### Number of terms

    for ax, lib in zip(axes[nlib:], lib_names.keys()):
        summary_df = summary_dict[lib]["summary_df"]
        qv = summary_df.drop(["Combined", "nan"], axis=1, level=0)
        qv = qv.xs("qvalue", level=2, axis=1)
        qv = qv.replace(np.nan, 1)
        qv = qv < qval

        qqv = qv.sum().reset_index()
        if pretty_print:
            qqv.replace({"Tool": pretty_print, "Metric": pretty_print}, inplace=True)
        qqv.index = qqv["Tool"] + "\n" + qqv["Metric"]
        qqv = qqv.drop(["Tool", "Metric"], axis=1)
        qqv = qqv.sort_values(by=qqv.columns[0], ascending=False)

        if ax == axes[nlib]:
            hue_order = {qqv.index[i]: palette[i] for i in range(len(qqv))}

        qqv["hue"] = qqv.index
        b = sns.barplot(
            data=qqv, x=qqv.index, y=qqv.iloc[:, 0], ax=ax, alpha=1, hue="hue", palette=palette, hue_order=hue_order
        )
        for i in b.containers:
            b.bar_label(
                i,
            )

        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=60, ha="right")
        ax.set(title=lib, ylabel="Terms", xlabel=None)

        # apply offset transform to all x ticklabels.
        offset = matplotlib.transforms.ScaledTranslation(12 / 72.0, 3 / 72.0, fig.dpi_scale_trans)
        for label in ax.xaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)
    for i in range(len(axes)):
        axes[i].annotate(
            chr(ord("A") + i),
            xy=(-0.08, 1.04),
            xycoords="axes fraction",
            weight="bold",
            va="center",
            ha="center",
            fontsize=18,
        )

    fig.tight_layout()
    fig.savefig(f"{figpath}/bars.{project_name}.{ext}")


def make_venn_plots(
    summary_dict,
    figpath: str,
    project_name: str,
    lib_names: Dict,
    metrics: List,
    tools: List,
    pretty_print: Dict = {},
    qval=0.05,
    ext: str = "pdf",
):
    if len(tools) == 1:
        save_empty(f"{figpath}/venn.methodcomp.{project_name}.{ext}", "")
    else:
        fig, ax = plt.subplots(len(lib_names), len(metrics), figsize=(len(metrics) * 4, len(lib_names.keys()) * 4))

        if len(lib_names) == 1:
            ax = np.expand_dims(ax, axis=0)  # Make it a 2D array with one row
        if len(metrics) == 1:
            ax = np.expand_dims(ax, axis=1)

        for i, lib in enumerate(lib_names.keys()):
            summary_df = summary_dict[lib]["summary_df"]
            sig_dict = get_sig_dict(summary_df, tools, metrics, qval=qval)
            for j, metric in enumerate(metrics):
                print(tools, metric)
                plot_venn(sig_dict, tools, metric, ax[i][j], pretty_print)
                ax[i][j].set_title(f"{pretty_print[metric] if pretty_print else metric} ({lib})", fontweight="bold")

        fig.tight_layout()
        fig.savefig(f"{figpath}/venn.methodcomp.{project_name}.{ext}")

    if len(metrics) == 1:
        save_empty(f"{figpath}/venn.metriccomp.{project_name}.{ext}", "")
    else:
        fig, ax = plt.subplots(len(lib_names), len(tools), figsize=(len(tools) * 4, len(lib_names.keys()) * 4))

        if len(lib_names) == 1:
            ax = np.expand_dims(ax, axis=0).T

        for i, lib in enumerate(lib_names.keys()):
            summary_df = summary_dict[lib]["summary_df"]
            sig_dict = get_sig_dict(summary_df, tools, metrics, qval=0.05)
            for j, tool in enumerate(tools):
                plot_venn(sig_dict, tool, metrics, ax[i][j], pretty_print)
                ax[i][j].set_title(f"{pretty_print[tool] if pretty_print else tool} ({lib})", fontweight="bold")

        fig.tight_layout()
        fig.savefig(f"{figpath}/venn.metriccomp.{project_name}.{ext}")


def make_upset_plots(
    summary_dict: Dict, lib_names: Dict, figpath: str, project_name: str, pretty_print: Dict = {}, ext: str = "pdf"
):
    for lib in lib_names.keys():
        outfile = f"{figpath}/upset.{lib}.{project_name}.{ext}"
        depth_df = summary_dict[lib]["depth_df"]
        if len(depth_df) < 1:
            print(f"No terms found for {lib}")
            save_empty(outfile, lib)
            continue
        memberships = depth_df["Configurations"]
        memberships_list = [categories.split("; ") for categories in memberships.values]
        upset_ready = from_memberships(memberships_list)
        upset_ready.index.names = [
            " ".join([pretty_print[i] if pretty_print else i for i in u.split(".")]) for u in upset_ready.index.names
        ]  # pretty print

        pd.options.mode.copy_on_write = False
        print(upset_ready.head())
        try:
            UpSet(upset_ready, subset_size="count", sort_by="cardinality", show_counts="{:,}").plot()
        except AttributeError:
            save_empty(outfile, lib)
        pd.options.mode.copy_on_write = True
        plt.savefig(outfile)


def lollipop_plots(
    df,
    lib,
    figpath,
    project_name,
    max_depth=0,
    ext="pdf",
    title=None,
    suffix="",
    x_val="SignedDepth",
    hue_subontology=False,
):
    outfile = f"{figpath}/lollipop{'.' + suffix if suffix != '' else ''}.{lib}.{project_name}.{ext}"
    if len(df) < 1:
        save_empty(outfile, lib)
        return

    df["logFDR"] = -np.log10(df["Combined FDR"])

    if x_val == "SignedDepth":
        df["SignedDepth"] = df["Depth"] * df["Direction"].apply(lambda x: 1 if x == "Up" else 0 if x == "Both" else -1)
        hue = "logFDR"
    elif x_val == "NegSignedlogFDR":
        df["NegSignedlogFDR"] = df["logFDR"] * df["Direction"].apply(
            lambda x: 1 if x == "Up" else 0 if x == "Both" else -1
        )
        max_depth = df["NegSignedlogFDR"].abs().max()
        if "Genes" in df:
            df["n_genes"] = df["Genes"].str.split(";").apply(lambda x: len(x))
            hue = "n_genes"
        else:
            hue = "Depth"
    else:
        raise Exception(f"Unknown x_val: {x_val}")

    if "ONTOLOGY" not in df:
        hue_subontology = False

    if hue_subontology:
        hue = "ONTOLOGY"

    ordered_df = df.sort_values(by=x_val)

    my_range = range(1, len(df.index) + 1)
    max_label_length = max(len(label) for label in ordered_df["Description"])

    with sns.axes_style("ticks"):
        fig_width = 4 + max_label_length * 0.08
        fig_height = max(3.6, len(df) // 3)
        fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))

    npg = npg = npg_palette()
    # hue subontology
    palette = {"BP": npg[0], "CC": npg[1], "MF": npg[2]}

    markers = {"BP": "o", "CC": "D", "MF": "s"}

    ax.hlines(y=my_range, xmin=0, xmax=ordered_df[x_val], zorder=98, color="grey")
    sns.scatterplot(
        data=ordered_df,
        x=x_val,
        y=range(1, 1 + len(ordered_df)),
        hue=hue,
        ax=ax,
        zorder=99,
        s=100,
        palette=palette if hue_subontology else None,
        style=hue if hue_subontology else None,
        markers=markers if hue_subontology else None,
    )

    ax.set_yticks(my_range, ordered_df["Description"])

    ax.set(title=f"Top {lib} terms {project_name}" if title is None else title)

    if hue_subontology:
        ax.legend(facecolor="white", bbox_to_anchor=(1.01, 1), loc="upper left")

    else:
        ### COLOR BAR
        cmap = sns.cubehelix_palette(as_cmap=True)
        norm = plt.Normalize(ordered_df[hue].min(), ordered_df[hue].max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        # Remove the legend and add a colorbar
        ax.get_legend().remove()
        fig.colorbar(sm, cax=cax)
        # clb.ax.set_title('This is a title')
        fig.axes[1].set(title=hue, xlabel="", ylabel="")

    # some extra spacing top and bottom
    ax.set_ylim(my_range[0] - 1, my_range[-1] + 1)

    if max_depth:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        # ax.set_xticks(range(-max_depth, max_depth+1))

    if max_depth > 0:
        low, high = ax.get_xlim()
        if ordered_df[x_val].max() > 0:
            high = max_depth + 0.5
        if ordered_df[x_val].min() < 0:
            low = -max_depth - 0.5
        ax.set_xlim(low, high)
        if low == -max_depth - 0.5 and high == max_depth + 0.5:
            ax.axvline(0, ls="--", color="grey")

    fig.tight_layout()
    fig.savefig(outfile)


def make_lollipop_plots(
    summary_dict: Dict,
    lib_names: Dict,
    figpath: str,
    project_name: str,
    top_terms: int = 30,
    qval: float = 0.05,
    depth_cutoff: int = 1,
    ext: str = "pdf",
    split_by_subontology=False,
    **kwargs,
):
    sns.set_theme(font_scale=1)

    for lib in lib_names.keys():
        d = summary_dict[lib]["depth_df"]

        if "ONTOLOGY" not in d:
            split_by_subontology = False

        if split_by_subontology:
            subonts = ["BP", "CC", "MF", "ALL"]
        else:
            subonts = [""]

        for ont in subonts:
            suffix = f"{ont}." if ont != "" else ""

            do = d[d["ONTOLOGY"] == ont] if ont not in ["", "ALL"] else d

            if len(do) < 1:
                print(f"No terms found for {lib}{suffix}")
                save_empty(f"{figpath}/lollipop.{suffix}{lib}.{project_name}.{ext}", lib)
                continue

            dd = do[(do["Depth"] > depth_cutoff) & (do["Combined FDR"] < qval)]
            if len(dd) < 1:
                print(f"No terms found for {lib}{suffix}")
                save_empty(f"{figpath}/lollipop.{suffix}{lib}.{project_name}.{ext}", lib)
                continue

            title = f"Top {lib}\n{project_name}\nDepth>{depth_cutoff}"
            lollipop_plots(
                dd.iloc[:top_terms], lib, figpath, project_name, ext=ext, title=title, suffix=f"{suffix}", **kwargs
            )

            # subset plots to enrichr terms
            if "Enrichr" in dd:
                title = f"Top {lib}\nEnrichr filtered\n{project_name}\nDepth>{depth_cutoff}"
                dd = dd[dd["Enrichr"]]
                if len(dd) > 1:
                    lollipop_plots(
                        dd.iloc[:top_terms],
                        lib,
                        figpath,
                        project_name,
                        ext=ext,
                        title=title,
                        suffix=f"{suffix}enrichr",
                        **kwargs,
                    )

            # subset plots to goesmsim fitlered terms
            if "Top_GO_Cluster" in dd:
                title = f"Top {lib}\nGOSemSim filtered\n{project_name}\nDepth>{depth_cutoff}"
                dd = dd[dd["Top_GO_Cluster"]]
                if len(dd) > 1:
                    lollipop_plots(
                        dd.iloc[:top_terms],
                        lib,
                        figpath,
                        project_name,
                        ext=ext,
                        title=title,
                        suffix=f"{suffix}gosemsim",
                        **kwargs,
                    )


# Save dummy figs to prevent SnakeMake jobs from failing when no terms are found
def save_empty(outfile, lib=""):
    fig, ax = plt.subplots(1, 1)
    if lib != "":
        ax.set_title(f"No terms found for {lib}")
    else:
        ax.set_title(f"No plot generated")
    fig.savefig(outfile)


if __name__ == "__main__":
    config_file_path = os.path.join("config", "config.yaml")

    with open(config_file_path, "r") as stream:
        config = yaml.safe_load(stream)

    libraries = config.get("libraries", [])
    lib_names = {os.path.splitext(os.path.basename(lib))[0]: lib for lib in libraries}

    tools = config.get("tools", [])
    metrics = config.get("metrics", [])
    pretty_print = config.get("pretty_print", {})
    project_name = config.get("project_name")
    qval = config.get("qval")
    max_depth = len(metrics) * len(tools)
    depth_cutoff = config.get("depth_cutoff_lollipop")
    x_val = config.get("x_val_lollipop")
    figpath = os.path.join("results", project_name, "figures")

    fig_formats = config.get("fig_formats", [])
    if isinstance(fig_formats, str):
        fig_formats = [fig_formats]

    summary_dict_file = os.path.join("results", project_name, "combined", f"syn.summary_dict.{project_name}.txt")

    with open(summary_dict_file, "rb") as f:
        summary_dict = pickle.load(f)

    for ext in fig_formats:
        make_bar_plots(
            summary_dict, figpath, project_name, lib_names, pretty_print, qval=qval, max_depth=max_depth, ext=ext
        )
        make_venn_plots(
            summary_dict, figpath, project_name, lib_names, metrics, tools, pretty_print, qval=qval, ext=ext
        )
        make_upset_plots(summary_dict, lib_names, figpath, project_name, pretty_print, ext=ext)
        make_lollipop_plots(
            summary_dict,
            lib_names,
            figpath,
            project_name,
            top_terms=30,
            qval=qval,
            depth_cutoff=depth_cutoff,
            max_depth=max_depth,
            ext=ext,
            x_val=x_val,
        )
