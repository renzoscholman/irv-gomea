import math
import os
import sys

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import seaborn as sn
import numpy as np
import pandas as pd
import scienceplots
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import mannwhitneyu

from rvgomea.defaults import DEFAULT_MAX_NUM_EVALUATIONS

# plt.style.use('science')

# Prevent scienceplots from being purged as import
# scienceplots.listdir(".")

cmap = matplotlib.colormaps['tab10']

LABELS = {
    # "univariate": "Univariate",
    # "full": "Full",
    # "uni-hg-gbo-without_clique_seeding-conditional": "Static-UCond-HG",
    # "mp-hg-gbo-without_clique_seeding-conditional": "Static-MCond-HG",
    "mp-hg-gbo-with_clique_seeding-conditional": "Static-CS-LM", #"Static-MCond-HG-CS",
    # "lt-fb-online-pruned": "FB-LT",
    # "uni-hg-fb_no_order-without_clique_seeding-conditional": "FB-UCond-HG",
    # "mp-hg-fb_no_order-without_clique_seeding-conditional": "FB-MCond-HG",
    "mp-hg-fb_no_order-with_clique_seeding-conditional": "FB-CS-LM", #"FB-MCond-HG-CS",
    "vkd-cma": "VkD-CMA-ES",
}

PROBLEMS = {
    "sphere": "Sphere",
    "rosenbrock": "Rosenbrock",
    "osoreb-big-strong": "OSoREB",
    "reb-grid": "REBGrid",
    "reb5-disjoint-pairs": "REB 5 Disjoint Pairs",
    "reb2-chain-weak": "REB 2 Weak",
    "reb2-chain-strong": "REB 2 Strong",
    "reb2-chain-alternating": "REB 2 Alternating",
    "reb5-no-overlap": "REB 5 No Overlap",
    "reb5-small-overlap": "REB 5 Small Overlap",
    "reb5-large-overlap": "REB 5 Large Overlap",
    "reb5-small-overlap-alternating": "REB 5 Alternating",
    "reb10-no-overlap": "REB 10 No Overlap",
    "reb10-small-overlap": "REB 10 Small Overlap",
    "reb10-large-overlap": "REB 10 Large Overlap",
    "reb10-small-overlap-alternating": "REB 10 Alternating",
}

MARKERS = {
    "univariate": "+",
    "full": "*",
    "uni-hg-gbo-without_clique_seeding-conditional": ">",
    "mp-hg-gbo-without_clique_seeding-conditional": "<",
    "mp-hg-gbo-with_clique_seeding-conditional": "d",
    "lt-fb-online-pruned": "s",
    "uni-hg-fb_no_order-without_clique_seeding-conditional": "2",
    "mp-hg-fb_no_order-without_clique_seeding-conditional": "x",
    "mp-hg-fb_no_order-with_clique_seeding-conditional": "X",
    "vkd-cma": "o",
}

COLOR_ORDER = {
    "univariate": 8,
    "full": 7,
    "uni-hg-gbo-without_clique_seeding-conditional": 6,
    "mp-hg-gbo-without_clique_seeding-conditional": 5,
    "mp-hg-gbo-with_clique_seeding-conditional": 4,
    "lt-fb-online-pruned": 2,
    "uni-hg-fb_no_order-without_clique_seeding-conditional": 1,
    "mp-hg-fb_no_order-without_clique_seeding-conditional": 0,
    "mp-hg-fb_no_order-with_clique_seeding-conditional": 9,
    "vkd-cma": 3,
}

COLORS = {}
for lm, color_index in COLOR_ORDER.items():
    COLORS[lm] = cmap(color_index / (len(COLOR_ORDER) - 1))

def extract_dfs(data_dir, start_pattern, filename, incremental, lm_override=False):
    dirs = [(i.replace(start_pattern, ""), os.path.join(data_dir, i, filename)) for i in os.listdir(data_dir) if i.startswith(start_pattern)]
    dfs = []
    for (dir, file) in dirs:
        if not os.path.exists(file):
            print(f"File {file} does not exist")
            continue
        if dir not in PROBLEMS:
            continue
        df = pd.read_csv(file)
        df["problem"] = PROBLEMS[dir]
        df["incremental"] = incremental
        for lm in df["linkage_model"].unique():
            if lm_override:
                if lm not in lm_override:
                    df = df.loc[df["linkage_model"] != lm]
                    continue
                df.loc[df["linkage_model"] == lm, "linkage_model"] = lm_override[lm]
            else:
                if lm not in LABELS:
                    df = df.loc[df["linkage_model"] != lm]
                    continue
                df.loc[df["linkage_model"] == lm, "linkage_model"] = LABELS[lm]
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True).drop(columns=["Unnamed: 0"])

if __name__ == '__main__':
    # Read results of iRV-GOMEA
    final_df = extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_dmd_halved_gom_only/", "scalability-bisection-incremental-", "aggregated_results.csv", True)
    extrapolated_dfs = [final_df]
    for i in range(5):
        extrapolated_dfs.append(extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_dmd_halved_gom_only/", "scalability-bisection-incremental-", f"extrapolated_results_{i}.csv", True))
    final_df = pd.concat(extrapolated_dfs, ignore_index=True)
    
    # Read results of RV-GOMEA
    original_df = extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/population_bisections/georgios/scalability-aggregated/", "scalability-bisection-", "aggregated_results.csv", False)
    original_df_reb10 = extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_rvg_gom_only/", "scalability-bisection-", "aggregated_results.csv", False)
    extrapolated_orig_dfs = [original_df, original_df_reb10]
    for i in range(5):
        extrapolated_orig_dfs.append(extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/population_bisections/georgios/", "scalability-bisection-", f"extrapolated_results_{i}.csv", False))
    for i in range(5):
        extrapolated_orig_dfs.append(extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_rvg_gom_only/", "scalability-bisection-", f"extrapolated_results_{i}.csv", False))
    original_df = pd.concat(extrapolated_orig_dfs, ignore_index=True)
    # print(final_df.loc[final_df["problem"] == "REB 10 Small Overlap"].groupby(["dimensionality", "linkage_model"])["population_size"].median())
    # print(original_df.loc[original_df["problem"] == "REB 10 Small Overlap"].groupby(["dimensionality", "linkage_model"])["corrected_num_evaluations"].median())
    # print(final_df.loc[final_df["problem"] == "REB 10 Small Overlap"].groupby(["dimensionality", "linkage_model"])["corrected_num_evaluations"].median())
    # exit(0)
    
    # Read results of VkD-CMA
    cma_df = original_df[original_df["linkage_model"] == "VkD-CMA-ES"]
    original_df = original_df[original_df["linkage_model"] != "VkD-CMA-ES"]
    original_df_reb10 = extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_rvg/", "scalability-bisection-", "aggregated_results.csv", False, {'vkd-cma': 'VkD-CMA-ES'})
    extrapolated_vkd_dfs = [cma_df, original_df_reb10]
    for i in range(5):
        extrapolated_vkd_dfs.append(extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_extrapolated_rvg/", "scalability-bisection-", f"extrapolated_results_{i}.csv", False, {'vkd-cma': 'VkD-CMA-ES'}))
    for i in range(5):
        reb10_df = extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_vkd/", "scalability-bisection-", f"extrapolated_results_{i}.csv", False, {'vkd-cma': 'VkD-CMA-ES'})
        reb10_df = reb10_df[reb10_df["problem"].str.contains("REB 10")]
        if len(reb10_df):
            extrapolated_vkd_dfs.append(reb10_df)
    cma_df = pd.concat(extrapolated_vkd_dfs, ignore_index=True)
    
    # cma_df = extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/population_bisections/georgios/scalability-aggregated/", "scalability-bisection-", "aggregated_results.csv", False, {'vkd-cma': 'VkD-CMA-ES'})
    # cma_georgios = original_df[original_df["linkage_model"] == "VkD-CMA-ES"]
    # cma_df = original_df[original_df["linkage_model"] == "VkD-CMA-ES"]
    # original_df = original_df[original_df["linkage_model"] != "VkD-CMA-ES"]
    # # cma_df = pd.concat([cma_df, cma_georgios], ignore_index=True)
    
    # extrapolated_orig_dfs = []
    # for i in range(5):
    #     extrapolated_orig_dfs.append(extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/population_bisections/georgios/", "scalability-bisection-", f"extrapolated_results_{i}.csv", False))
    # extrapolated_orig_dfs = pd.concat(extrapolated_orig_dfs, ignore_index=True)
    
    # extrapolated_vkd_dfs = []
    # extrapolated_rvg_dfs = []
    # for i in range(5):
    #     extrapolated_vkd_dfs.append(extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/population_bisections/georgios/", "scalability-bisection-", f"extrapolated_results_{i}.csv", False, {'vkd-cma': 'VkD-CMA-ES'}))
    # for i in range(5):
    #     extrapolated_rvg_dfs.append(extract_dfs("/mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_extrapolated_rvg/", "scalability-bisection-", f"extrapolated_results_{i}.csv", False))
    # extrapolated_vkd_dfs = pd.concat(extrapolated_vkd_dfs, ignore_index=True)
    # extrapolated_rvg_dfs = pd.concat(extrapolated_rvg_dfs, ignore_index=True)
    # extrapolated_rvg_dfs = extrapolated_rvg_dfs[extrapolated_rvg_dfs["linkage_model"] == "VkD-CMA-ES"]
    # cma_df = pd.concat([cma_df, extrapolated_vkd_dfs, extrapolated_rvg_dfs], ignore_index=True)
    # print(extrapolated_orig_dfs.loc[extrapolated_orig_dfs["problem"] == "Rosenbrock"].groupby(["linkage_model", "dimensionality"])["population_size"].median())
    # print(original_df.loc[original_df["problem"] == "Rosenbrock"].groupby(["linkage_model", "dimensionality"])["population_size"].median())
    # exit(0)
    # print(extrapolated_dfs)
    # exit(0)
    
    # print(original_df["problem"].unique())
    # print(cma_df["problem"].unique())
    # exit(0)
    
    # print(final_df)
    # print(cma_df)
    # exit(0)
    
    # Remove outliers
    final_df = final_df[final_df["corrected_num_evaluations"] < 1e8]
    original_df = original_df[original_df["corrected_num_evaluations"] < 1e8]
    cma_df = cma_df[cma_df["corrected_num_evaluations"] < 1e8]
    # print(cma_df.loc[cma_df["problem"] == "REBGrid"])
    # exit(0)
    
    # Build ratio DF
    ratio_df = pd.DataFrame(columns=["problem", "linkage_model", "dimensionality", "ratio"])
    eval_df = pd.DataFrame(columns=["problem", "linkage_model", "dimensionality", "evaluations"])
    for key, df in final_df.groupby(["problem", "linkage_model", "dimensionality"]):
        new_eval_df = pd.DataFrame({"problem": key[0], "linkage_model": "iRV-GOMEA "+key[1], "dimensionality": key[2], "evaluations": df["corrected_num_evaluations"].median()}, index=[0])
        eval_df = pd.concat([eval_df, new_eval_df], ignore_index=True)
        
        df_original = original_df[(original_df["problem"] == key[0]) & (original_df["linkage_model"] == key[1]) & (original_df["dimensionality"] == key[2])]
        if len(df_original) == 0:
            continue
        ratio = df["corrected_num_evaluations"].median() / df_original["corrected_num_evaluations"].median()
        new_df = pd.DataFrame({"problem": key[0], "linkage_model": key[1], "dimensionality": key[2], "ratio": ratio}, index=[0])
        ratio_df = pd.concat([ratio_df, new_df], ignore_index=True)
        
        new_eval_df = pd.DataFrame({"problem": key[0], "linkage_model": "RV-GOMEA "+key[1], "dimensionality": key[2], "evaluations": df_original["corrected_num_evaluations"].median()}, index=[0])
        eval_df = pd.concat([eval_df, new_eval_df], ignore_index=True)
        
        if key[1].startswith("FB-"):
            df_cma = cma_df[(cma_df["problem"] == key[0]) & (cma_df["dimensionality"] == key[2])]
            if len(df_cma) == 0:
                print(f"No CMA data for {key[0]} {key[2]}")
                continue
            ratio = df["corrected_num_evaluations"].median() / df_cma["corrected_num_evaluations"].median()
            new_df = pd.DataFrame({"problem": key[0], "linkage_model": "VkD-CMA-ES", "dimensionality": key[2], "ratio": ratio}, index=[0])
            ratio_df = pd.concat([ratio_df, new_df], ignore_index=True)
        
            new_eval_df = pd.DataFrame({"problem": key[0], "linkage_model": "VkD-CMA-ES", "dimensionality": key[2], "evaluations": df_cma["corrected_num_evaluations"].median()}, index=[0])
            eval_df = pd.concat([eval_df, new_eval_df], ignore_index=True)
            
    # for key, df in extrapolated_dfs.groupby(["problem", "linkage_model", "dimensionality"]):
    #     new_eval_df = pd.DataFrame({"problem": key[0], "linkage_model": "iRV-GOMEA "+key[1], "dimensionality": key[2], "evaluations": df["corrected_num_evaluations"].median()}, index=[0])
    #     eval_df = pd.concat([eval_df, new_eval_df], ignore_index=True)
        
    #     df_original_extrapolated = extrapolated_orig_dfs[(extrapolated_orig_dfs["problem"] == key[0]) & (extrapolated_orig_dfs["linkage_model"] == key[1]) & (extrapolated_orig_dfs["dimensionality"] == key[2])]
    #     if len(df_original_extrapolated) == 0:
    #         continue
    #     ratio = df["corrected_num_evaluations"].median() / df_original_extrapolated["corrected_num_evaluations"].median()
    #     new_df = pd.DataFrame({"problem": key[0], "linkage_model": key[1], "dimensionality": key[2], "ratio": ratio}, index=[0])
    #     ratio_df = pd.concat([ratio_df, new_df], ignore_index=True)
    
    #     new_eval_df = pd.DataFrame({"problem": key[0], "linkage_model": "RV-GOMEA "+key[1], "dimensionality": key[2], "evaluations": df_original_extrapolated["corrected_num_evaluations"].median()}, index=[0])
    #     eval_df = pd.concat([eval_df, new_eval_df], ignore_index=True)
        
    #     if key[1].startswith("FB-"):
    #         df_cma = cma_df[(cma_df["problem"] == key[0]) & (cma_df["dimensionality"] == key[2])]
    #         if len(df_cma) == 0:
    #             print(f"No CMA data for {key[0]} {key[2]}")
    #             continue
    #         ratio = df["corrected_num_evaluations"].median() / df_cma["corrected_num_evaluations"].median()
    #         new_df = pd.DataFrame({"problem": key[0], "linkage_model": "VkD-CMA-ES", "dimensionality": key[2], "ratio": ratio}, index=[0])
    #         ratio_df = pd.concat([ratio_df, new_df], ignore_index=True)
        
    #         new_eval_df = pd.DataFrame({"problem": key[0], "linkage_model": "VkD-CMA-ES", "dimensionality": key[2], "evaluations": df_cma["corrected_num_evaluations"].median()}, index=[0])
    #         eval_df = pd.concat([eval_df, new_eval_df], ignore_index=True)
    
    output_dir = "figures/final_bisections_gom_only/"
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    col_order = PROBLEMS.values()
    hue_order = sorted(ratio_df["linkage_model"].unique())
    # sn.set(rc={"xtick.bottom" : True, "xtick.top" : True, "ytick.left" : True, "ytick.right" : True})
    fg = sn.relplot(
        data=ratio_df, x="dimensionality", y="ratio", col="problem", palette="bright",
        hue="linkage_model", style="linkage_model", kind="line", dashes=False, markers=True,
        col_wrap=4, col_order=col_order, hue_order=hue_order, height=1.6, aspect=2**0.5, facet_kws={"sharey": False}
    )
    # fg = sn.FacetGrid(ratio_df, col="problem", hue="linkage_model", col_wrap=4, col_order=col_order, hue_order=hue_order, height=2**0.5, aspect=2**0.5, sharey=False)
    # fg.map(sn.lineplot, "dimensionality", "ratio", errorbar=None, markers=True)
    sn.move_legend(fg, "upper center", bbox_to_anchor=(0.45, 0.04), ncol=3, fontsize="small", columnspacing=1.0, title="")
    # fg.add_legend(loc='upper center',bbox_to_anchor=(0.35, 0.025), ncol=3)
    fg.set(xscale='log')
    # fg.set(yticks = [0.1, 1.0], yticklabels = [0.1, 1.0])
    fg.refline(y=1.0)
    fg.set_titles(col_template="{col_name}")
    fg.set_axis_labels("$\\ell$", "Ratio")
    for ax in fg.axes:
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.tick_params(axis='y', which='both', labelright=False, right=True, labelleft=True, left=True)
        ax.tick_params(axis='x', which='both', labeltop=False, top=True, labelbottom=False, bottom=True)
    for ax in fg.axes[-4:]:
        ax.tick_params(axis='y', which='both', labelright=False, right=True, labelleft=True, left=True)
        ax.tick_params(axis='x', which='both', labeltop=False, top=True, labelbottom=True, bottom=True)
    # plt.tight_layout()
    plt.savefig(output_dir+"gecco_ratios_cond.svg", bbox_inches='tight')
    plt.close()
    
    # exit(0)
    
    col_order = PROBLEMS.values()
    hue_order = sorted(eval_df["linkage_model"].unique())
    fg = sn.relplot(
        data=eval_df, x="dimensionality", y="evaluations", col="problem", palette="bright",
        hue="linkage_model", style="linkage_model", kind="line", dashes=False, markers=True,
        col_wrap=4, col_order=col_order, hue_order=hue_order, height=1.6, aspect=2**0.5, facet_kws={"sharey": False}
    )
    # fg = sn.FacetGrid(eval_df, col="problem", hue="linkage_model", col_wrap=4, col_order=col_order, hue_order=hue_order, height=2**0.5, aspect=2**0.5, sharey=False)
    # fg.map(sn.lineplot, "dimensionality", "evaluations", errorbar=None, markers='o')
    sn.move_legend(fg, "upper center", bbox_to_anchor=(0.4, 0.04), ncol=5, fontsize="small", columnspacing=1.0, title="")
    fg.set(yscale='log', xscale='log')
    # fg.set(yticks = [0.1, 1.0], yticklabels = [0.1, 1.0])
    fg.set_titles(col_template="{col_name}")
    fg.set_axis_labels("Dimensionality $\\ell$", "Evaluations")
    fg.tight_layout()
    for ax in fg.axes:
        # ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.tick_params(axis='y', which='both', labelright=False, right=True, labelleft=True, left=True)
        ax.tick_params(axis='x', which='both', labeltop=False, top=True, labelbottom=False, bottom=True)
    for ax in fg.axes[-4:]:
        ax.tick_params(axis='y', which='both', labelright=False, right=True, labelleft=True, left=True)
        ax.tick_params(axis='x', which='both', labeltop=False, top=True, labelbottom=True, bottom=True)
    plt.savefig(output_dir+"evals.svg", bbox_inches='tight')
    plt.close()
    
    original_df["linkage_model"] = original_df["linkage_model"].apply(lambda x: "RV-GOMEA "+x)
    final_df["linkage_model"] = final_df["linkage_model"].apply(lambda x: "iRV-GOMEA "+x)
    # extrapolated_dfs["linkage_model"] = extrapolated_dfs["linkage_model"].apply(lambda x: "iRV-GOMEA "+x)
    # extrapolated_orig_dfs["linkage_model"] = extrapolated_orig_dfs["linkage_model"].apply(lambda x: "RV-GOMEA "+x if x != "VkD-CMA" else x)
    pop_df = pd.concat([final_df, original_df, cma_df], ignore_index=True)
    
    col_order = PROBLEMS.values()
    hue_order = sorted(pop_df["linkage_model"].unique())
    fg = sn.relplot(
        data=pop_df, x="dimensionality", y="population_size", col="problem", palette="bright",
        hue="linkage_model", style="linkage_model", kind="line", dashes=False, markers=True,
        col_wrap=4, col_order=col_order, hue_order=hue_order, height=1.6, aspect=2**0.5, 
        errorbar=None, legend="full", facet_kws={"legend_out":True}
    )
    # fg = sn.FacetGrid(eval_df, col="problem", hue="linkage_model", col_wrap=4, col_order=col_order, hue_order=hue_order, height=2**0.5, aspect=2**0.5, sharey=False)
    # fg.map(sn.lineplot, "dimensionality", "evaluations", errorbar=None, markers='o')
    # fg.add_legend(loc='upper center',bbox_to_anchor=(0.35, 0.025), ncol=5)
    sn.move_legend(fg, "upper center", bbox_to_anchor=(0.4, 0.04), ncol=5, fontsize="small", columnspacing=1.0, title="")
    # leg = fg._legend
    # leg.set_bbox_to_anchor([0.5, 1.0])
    fg.set(yscale='log', xscale='log')
    # fg.set(yticks = [0.1, 1.0], yticklabels = [0.1, 1.0])
    fg.set_titles(col_template="{col_name}")
    fg.set_axis_labels("$\\ell$", "Population Size")
    for ax in fg.axes:
        # ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.tick_params(axis='y', which='both', labelright=False, right=True, labelleft=True, left=True)
        ax.tick_params(axis='x', which='both', labeltop=False, top=True, labelbottom=False, bottom=True)
    for ax in fg.axes[-4:]:
        ax.tick_params(axis='y', which='both', labelright=False, right=True, labelleft=True, left=True)
        ax.tick_params(axis='x', which='both', labeltop=False, top=True, labelbottom=True, bottom=True)
    fg.tight_layout()
    # plt.tick_params(axis='y', which='both', labelright='off', labelleft='on')
    # plt.tick_params(axis='x', which='both', labeltop='off', labelbottom='on')
    plt.savefig(output_dir+"gecco_popsize.svg", bbox_inches='tight')
    plt.close()
    
    # fg = sn.FacetGrid(final_df, col="problem", hue="linkage_model", col_wrap=4, col_order=col_order, hue_order=hue_order)
    # fg.map(sn.lineplot, "dimensionality", "corrected_num_evaluations", errorbar=None).set(yscale = 'log')
    # fg.add_legend()
    # plt.savefig(output_dir+"incremental_evals.svg")
    # plt.close()
    
    # fg = sn.FacetGrid(original_df, col="problem", hue="linkage_model", col_wrap=4, col_order=col_order, hue_order=hue_order)
    # fg.map(sn.lineplot, "dimensionality", "corrected_num_evaluations", errorbar=None).set(yscale = 'log')
    # fg.add_legend()
    # plt.savefig(output_dir+"non_incremental_evals.svg")
    # plt.close()
    
    # fg = sn.FacetGrid(final_df, col="problem", hue="linkage_model", col_wrap=4, col_order=col_order, hue_order=hue_order)
    # fg.map(sn.lineplot, "dimensionality", "population_size", errorbar=None).set(yscale = 'log')
    # fg.add_legend()
    # plt.savefig(output_dir+"incremental_popsize.svg")
    # plt.close()
    
    # fg = sn.FacetGrid(original_df, col="problem", hue="linkage_model", col_wrap=4, col_order=col_order, hue_order=hue_order)
    # fg.map(sn.lineplot, "dimensionality", "population_size", errorbar=None).set(yscale = 'log')
    # fg.add_legend()
    # plt.savefig(output_dir+"non_incremental_popsize.svg")
    # plt.close()
    
    # # figsize = (10, 6)
    # figsize = (20, 12)
    # fig, axs = plt.subplots(3, 4, figsize=figsize, sharex=True, sharey=True)

    # for i, (problem_id, problem_label) in enumerate(zip(problem_ids, problem_labels)):
    #     ax = axs[i // 4, i % 4]
    #     line_style = "--" if metric == "population_size" else "-"
    #     ax.plot(extra_dimensions, extra_values, linestyle=line_style, color=COLORS[lm], marker=MARKERS[lm],
    #             label="_" + LABELS[lm])

    #     handles, labels = make_one_plot(ax, base_directory + problem_id, problem_id, linkage_models, metric, problem_label)
    
    #     if metric == "corrected_num_evaluations":
    #         ax.set_ylim(1e2, 1e7)

    # fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.5, 1.01), ncol=math.ceil(len(linkage_models) / 2))

    # # Global labels
    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.xlabel('Dimensionality')
    # plt.ylabel(metric_label)

    # if len(problem_ids) == 1:
    #     prefix = "set_cover_scalability"
    # else:
    #     prefix = "scalability"

    # plt.tight_layout()
    # plt.savefig(os.path.join(plot_directory, f'{prefix}_{metric}.png'), bbox_inches='tight')
    # # plt.savefig(os.path.join(plot_directory, f'{prefix}_{metric}.pdf'), bbox_inches='tight')
    # plt.clf()
    # plt.close(fig)
    # main(sys.argv[1], sys.argv[2].split(","), sys.argv[3].split(","), sys.argv[4].split(","))
