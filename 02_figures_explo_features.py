#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures to explore the features for calibration

@author: LucieBourreau
@date: 2025/02/17
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import numpy as np

df = pd.read_csv("./merged_LOKI2013_ecotaxa_masks_features.csv")

summary_table = df.pivot_table(
    index='object_date',            
    columns='object_annotation_category',        
    values='object_id',             
    aggfunc='count',                
    fill_value=0                    
)

summary_table.to_csv("summary_table_LOKI2013.csv")

#### Needed for figures

colors = {
    'Calanus hyperboreus': '#7fcdbb',
    'Calanus glacialis': '#feb24c',
}

species_list = ["Calanus hyperboreus", "Calanus glacialis"]

lipids_min = df['total_lipids_ugC'].min()
lipids_max = df['total_lipids_ugC'].max() + 50

 
def plot_lipid_vs_fullness_hist(ax_main, ax_histx, ax_histy, df, stage, lipids_min, lipids_max, equation_ref):
    """
    Plot lipid content in carbon against fullness ratio in a scatter plot with the associated
    histograms of each trait. 

    Parameters
    ----------
    ax_main : fig.subplot
        Axis for the scatter plot.
    ax_histx : fig.subplot
        Axis for the histogram associated to the x axis (lipids).
    ax_histy : fig.subplot
        Axis for the histogram associated to the y axis (fullness).
    df : data.frame
        Dataset with the traits.
    stage : str
        Specify the stage ('civstage', 'cvstage' or 'female').
    lipids_min : float
        Min value of lipids (for axis limits).
    lipids_max : float
         Max value of lipids (for axis limits).
    equation_ref : str
        With which method the fullness has been computed ('Carbon Area' or 'Carbon Fullness').

    Returns
    -------
    None.

    """
    data_at_stage = df[df['object_annotation_category'].str.contains(stage, case=False, na=False)]
    
    data_hyperboreus = data_at_stage[data_at_stage['object_annotation_category'].str.contains('hyperboreus', case=False, na=False)]
    data_glacialis = data_at_stage[data_at_stage['object_annotation_category'].str.contains('glacialis', case=False, na=False)]
        
    if equation_ref == "Carbon Area":
         
        hyp_total_lipids_ugC = data_hyperboreus['total_lipids_ugC']
        gla_total_lipids_ugC = data_glacialis['total_lipids_ugC']
         
        hyp_fullness_ratio = data_hyperboreus['fullness_ratio_carbon_area']
        gla_fullness_ratio = data_glacialis['fullness_ratio_carbon_area']
        
    if equation_ref == "Carbon Volume":
          
        hyp_total_lipids_ugC = data_hyperboreus['total_lipids_ugC']
        gla_total_lipids_ugC = data_glacialis['total_lipids_ugC']
          
        hyp_fullness_ratio = data_hyperboreus['fullness_ratio_carbon_volume']
        gla_fullness_ratio = data_glacialis['fullness_ratio_carbon_volume']

    # Histogramme des lipides (en haut)
    ax_histx.hist([hyp_total_lipids_ugC,
                   gla_total_lipids_ugC],
                   bins=20,
                   range=(lipids_min, lipids_max),
                   color=[colors['Calanus hyperboreus'], colors['Calanus glacialis']],
                   #stacked=True, 
                   alpha=0.7)
                   
    ax_histx.spines['top'].set_visible(False)
    ax_histx.spines['right'].set_visible(False)
    ax_histx.tick_params(axis='x', labelbottom=False)
    ax_histx.tick_params(axis='y', left=False)

    # Scatter plot (nuage de points)
    ax_main.scatter(hyp_total_lipids_ugC, 
                    hyp_fullness_ratio,
                    color=colors['Calanus hyperboreus'], 
                    label='C. hyperboreus', 
                    s=20, 
                    alpha=0.7)
    
    ax_main.scatter(gla_total_lipids_ugC, 
                    gla_fullness_ratio,
                    color=colors['Calanus glacialis'], 
                    label='C. glacialis', 
                    s=20, 
                    alpha=0.7)
    
    ax_main.set_xlabel("Total Lipids (ugC)")
    ax_main.set_ylabel(f"Fullness Ratio ({equation_ref})")
        
    ax_main.set_xlim([lipids_min, lipids_max])
    ax_main.set_ylim([0, 1])
    
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)

    # Histogramme de la fullness (à droite)
    ax_histy.hist([hyp_fullness_ratio,
                   gla_fullness_ratio],
                   bins=20, 
                   range=(0, 1),
                   color=[colors['Calanus hyperboreus'], colors['Calanus glacialis']],
                   #stacked=True, 
                   alpha=0.7, 
                   orientation='horizontal')
    
    ax_histy.spines['top'].set_visible(False)
    ax_histy.spines['right'].set_visible(False)
    ax_histy.tick_params(axis='y', labelleft=False)
    ax_histy.tick_params(axis='x', bottom=False)



#### Figure 1: scatter plot and histograms per stage CV and Females (Carbon Area)

fig = plt.figure(figsize=(14, 6)) 
gs = gridspec.GridSpec(3, 6, figure=fig, width_ratios=[4, 1.5, 0.2, 4, 1.5, 0.5], height_ratios=[0.5, 1.5, 4], wspace=0.3, hspace=0.3)

# **Titles**
ax_title_cv = fig.add_subplot(gs[0, 0:2], frameon=False)
ax_title_fem = fig.add_subplot(gs[0, 2:5], frameon=False)

for ax in [ax_title_cv, ax_title_fem]:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

ax_title_cv.text(0.5, 0.1, "C5 stage - August 2013", ha='center', va='center', fontsize=14, fontweight='bold')
ax_title_fem.text(0.5, 0.1, "Females - August 2013", ha='center', va='center', fontsize=14, fontweight='bold')

# **Figure for CV stage**
ax_histx_cv = fig.add_subplot(gs[1, 0])
ax_main_cv = fig.add_subplot(gs[2, 0])
ax_histy_cv = fig.add_subplot(gs[2, 1], sharey=ax_main_cv)

plot_lipid_vs_fullness_hist(ax_main_cv, ax_histx_cv, ax_histy_cv, df, 'cvstage', lipids_min, lipids_max, "Carbon Area")

# **Figure for Females**
ax_histx_fem = fig.add_subplot(gs[1, 3])
ax_main_fem = fig.add_subplot(gs[2, 3])
ax_histy_fem = fig.add_subplot(gs[2, 4], sharey=ax_main_fem)

plot_lipid_vs_fullness_hist(ax_main_fem, ax_histx_fem, ax_histy_fem, df, 'female', lipids_min, lipids_max, "Carbon Area")

# **Add legend**
legend_ax = fig.add_subplot(gs[:, 5], frameon=False)
legend_ax.set_xticks([])
legend_ax.set_yticks([])
legend_ax.set_frame_on(False)

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#7fcdbb', markersize=8, label='C. hyperboreus'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#feb24c', markersize=8, label='C. glacialis'),
]

legend_ax.legend(handles=legend_elements, loc='center', frameon=False, fontsize=10)

plt.tight_layout()

plt.savefig("Lipids_against_fullness_hist_CV_Females_August_LOKI2013_Carbon_Area.png")

plt.show()


#### Figure 2: scatter plot and histograms per stage CV and Females (Carbon Volume)

fig = plt.figure(figsize=(14, 6)) 
gs = gridspec.GridSpec(3, 6, figure=fig, width_ratios=[4, 1.5, 0.2, 4, 1.5, 0.5], height_ratios=[0.5, 1.5, 4], wspace=0.3, hspace=0.3)

# **Titles**
ax_title_cv = fig.add_subplot(gs[0, 0:2], frameon=False)
ax_title_fem = fig.add_subplot(gs[0, 2:5], frameon=False)

for ax in [ax_title_cv, ax_title_fem]:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

ax_title_cv.text(0.5, 0.1, "C5 stage - August 2013", ha='center', va='center', fontsize=14, fontweight='bold')
ax_title_fem.text(0.5, 0.1, "Females - August 2013", ha='center', va='center', fontsize=14, fontweight='bold')

# **Figure for CV stage**
ax_histx_cv = fig.add_subplot(gs[1, 0])
ax_main_cv = fig.add_subplot(gs[2, 0])
ax_histy_cv = fig.add_subplot(gs[2, 1], sharey=ax_main_cv)

plot_lipid_vs_fullness_hist(ax_main_cv, ax_histx_cv, ax_histy_cv, df, 'cvstage', lipids_min, lipids_max, "Carbon Volume")

# **Figure for Females**
ax_histx_fem = fig.add_subplot(gs[1, 3])
ax_main_fem = fig.add_subplot(gs[2, 3])
ax_histy_fem = fig.add_subplot(gs[2, 4], sharey=ax_main_fem)

plot_lipid_vs_fullness_hist(ax_main_fem, ax_histx_fem, ax_histy_fem, df, 'female', lipids_min, lipids_max, "Carbon Volume")

ax_main_cv.set_ylim(0,1)
ax_main_fem.set_ylim(0,1)

# **Add legend**
legend_ax = fig.add_subplot(gs[:, 5], frameon=False)
legend_ax.set_xticks([])
legend_ax.set_yticks([])
legend_ax.set_frame_on(False)

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#7fcdbb', markersize=8, label='C. hyperboreus'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#feb24c', markersize=8, label='C. glacialis'),
]

legend_ax.legend(handles=legend_elements, loc='center', frameon=False, fontsize=10)

plt.tight_layout()

plt.savefig("Lipids_against_fullness_hist_CV_Females_August_LOKI2013_Carbon_Volume.png")

plt.show()
