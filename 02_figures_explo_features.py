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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend_handler import HandlerTuple
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


def plot_lipid_vs_fullness_hist(ax_main, ax_histx, ax_histy, df, stage, lipids_min, lipids_max):
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
        
          
    hyp_total_lipids_ugC = data_hyperboreus['total_lipids_ugC']
    gla_total_lipids_ugC = data_glacialis['total_lipids_ugC']
      
    hyp_fullness_ratio = data_hyperboreus['fullness_ratio_carbon']
    gla_fullness_ratio = data_glacialis['fullness_ratio_carbon']

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
    ax_main.set_ylabel("Fullness Ratio")
        
    ax_main.set_xlim([0, lipids_max+20])
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


#### Figure 1: scatter plot and histograms per stage CV and Females

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

plot_lipid_vs_fullness_hist(ax_main_cv, ax_histx_cv, ax_histy_cv, df, 'cvstage', lipids_min, lipids_max)

# **Figure for Females**
ax_histx_fem = fig.add_subplot(gs[1, 3])
ax_main_fem = fig.add_subplot(gs[2, 3])
ax_histy_fem = fig.add_subplot(gs[2, 4], sharey=ax_main_fem)

plot_lipid_vs_fullness_hist(ax_main_fem, ax_histx_fem, ax_histy_fem, df, 'female', lipids_min, lipids_max)

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

plt.savefig("Lipids_against_fullness_hist_CV_Females_August_LOKI2013.png")

plt.show()


#### Figure 2: scatter plot and histograms per stage CIV, CV and Females

fig = plt.figure(figsize=(10, 14)) 
gs = gridspec.GridSpec(11, 3, figure=fig, width_ratios=[4, 1.5, 0.5], height_ratios=[1, 1.5, 4, 0.2, 1, 1.5, 4, 0.2, 1, 1.5, 4], wspace=0.3, hspace=0.3)


# **Titles**
ax_title_civ = fig.add_subplot(gs[0, 0:2], frameon=False)
ax_title_cv = fig.add_subplot(gs[4, 0:2], frameon=False)
ax_title_fem = fig.add_subplot(gs[8, 0:2], frameon=False)

for ax in [ax_title_civ, ax_title_cv, ax_title_fem]:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

ax_title_civ.text(0.5, 0.1, "C4 stage - August 2013", ha='center', va='center', fontsize=12, fontweight='bold')
ax_title_cv.text(0.5, 0.1, "C5 stage - August 2013", ha='center', va='center', fontsize=12, fontweight='bold')
ax_title_fem.text(0.5, 0.1, "Females - August 2013", ha='center', va='center', fontsize=12, fontweight='bold')

# **Figure for CIV stage**
ax_histx_civ = fig.add_subplot(gs[1, 0])
ax_main_civ = fig.add_subplot(gs[2, 0])
ax_histy_civ = fig.add_subplot(gs[2, 1], sharey=ax_main_civ)

plot_lipid_vs_fullness_hist(ax_main_civ, ax_histx_civ, ax_histy_civ, df, 'civstage', lipids_min, lipids_max)


# **Figure for CV stage**
ax_histx_cv = fig.add_subplot(gs[5, 0])
ax_main_cv = fig.add_subplot(gs[6, 0])
ax_histy_cv = fig.add_subplot(gs[6, 1], sharey=ax_main_cv)

plot_lipid_vs_fullness_hist(ax_main_cv, ax_histx_cv, ax_histy_cv, df, 'cvstage', lipids_min, lipids_max)

# **Figure for Females**
ax_histx_fem = fig.add_subplot(gs[9, 0])
ax_main_fem = fig.add_subplot(gs[10, 0])
ax_histy_fem = fig.add_subplot(gs[10, 1], sharey=ax_main_fem)

plot_lipid_vs_fullness_hist(ax_main_fem, ax_histx_fem, ax_histy_fem, df, 'female', lipids_min, lipids_max)

ax_main_civ.set_ylim(0,1)
ax_main_cv.set_ylim(0,1)
ax_main_fem.set_ylim(0,1)

# **Add legend**
legend_ax = fig.add_subplot(gs[:, 2], frameon=False)
legend_ax.set_xticks([])
legend_ax.set_yticks([])
legend_ax.set_frame_on(False)

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#7fcdbb', markersize=8, label='C. hyperboreus'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#feb24c', markersize=8, label='C. glacialis'),
]

legend_ax.legend(handles=legend_elements, loc='center', frameon=False, fontsize=10)

fig.subplots_adjust(top=1, bottom=0.05) 
plt.tight_layout()

plt.savefig("Lipids_against_fullness_hist_CIV_CV_Females_August_LOKI2013.png")

plt.show()

#### Figure 3: Fullness ratio (area) against prosome length 

df['stage'] = df['object_annotation_category'].str.extract(r'^(.*?)\+')[0]
df['species'] = df['object_annotation_category'].str.extract(r'<([^<]+)$')[0]

smooth_green_blue = LinearSegmentedColormap.from_list(
    "smooth_green_blue",
    ['azure', '#a0ddc0', '#5fcab4', '#006eca']
)
colors_hyp = [smooth_green_blue(i) for i in np.linspace(0.0, 0.9, 7)]

smooth_orange_red = LinearSegmentedColormap.from_list(
    "smooth_orange_red",
    ['#fffff0', '#FFEE74', '#feb24c', '#e34a33']
)
colors_gla = [smooth_orange_red(i) for i in np.linspace(0.0, 0.9, 7)]

# Take only the 3 last values for the 3 older stages
colors_hyp_stage = colors_hyp[-3:]
colors_gla_stage = colors_gla[-3:]

stage_labels = ['C4', 'C5', 'AF']

stage_to_idx = {
    'civstage': 0,  #C4
    'cvstage': 1,   #C5
    'female': 2     #AF
}

markers = {
    'civstage': 'D',
    'cvstage': '^',
    'female': 'o'
}

fig, ax = plt.subplots(figsize=(8, 6))
for species in ['Calanus hyperboreus', 'Calanus glacialis']:
    for stage, group in df[df['species'] == species].groupby('stage'):
        idx = stage_to_idx.get(stage, None)
        if idx is None:
            color = 'grey'
        else:
            if species == 'Calanus hyperboreus':
                color = colors_hyp_stage[idx]
            elif species == 'Calanus glacialis':
                color = colors_gla_stage[idx]
            else:
                color = 'grey'
            marker = markers.get(stage, 'o')

        ax.scatter(
            group['prosome_major_axis_mm'],
            group['fullness_ratio_area'],
            color=color,
            # facecolor=color,
            edgecolor='none',
            marker=marker,
            label=f"{stage} - {species}",
            alpha=0.8,
            s=30
        )

# Legend and labels

legend_handles = []
stage_keys = ['civstage', 'cvstage', 'female'] 
    
for stage in stage_keys:
    idx = stage_to_idx[stage]
    # Patch Hyperboreus 
    handle_hyp = Line2D(
        [0], [0],
        marker=markers[stage],
        markerfacecolor=colors_hyp_stage[idx],
        markeredgecolor='none',
        markersize=7,
        linestyle='None'
    )
    # Patch Glacialis
    handle_gla = Line2D(
        [0], [0],
        marker=markers[stage],
        markerfacecolor=colors_gla_stage[idx],
        markeredgecolor='none',
        markersize=7,
        linestyle='None'
    )
    legend_handles.append((handle_hyp, handle_gla))
    
fig.legend(
    legend_handles,
    stage_labels,
    handler_map={tuple: HandlerTuple(ndivide=None)},
    loc='center right',
    bbox_to_anchor=(1, 0.5),
    frameon=False,
)

fig.text(
    0.91,    
    0.56,    
    r'$C.\ hyperboreus$', 
    rotation=90, 
    fontsize=10, 
    va='bottom', 
    ha='left',
    style='italic'
)

fig.text(
    0.935,    
    0.56,    
    r'$C.\ glacialis$', 
    rotation=90, 
    fontsize=10, 
    va='bottom', 
    ha='left',
    style='italic'
)

ax.spines['top'].set_visible(False) 
ax.spines['right'].set_visible(False)
ax.set_xlim([2,8])
ax.set_ylim([0,0.8])
ax.set_xlabel('Prosome length (mm)')
ax.set_ylabel('Lipid fullness (LSA:PA)')

plt.savefig("FullnessRatio_against_PL_CIV_CV_Females_August_LOKI2013.png", dpi=600)
plt.show()


#### Figure 4: lipid sac volume against prosome length

fig, ax = plt.subplots(figsize=(8, 6))
for species in ['Calanus hyperboreus', 'Calanus glacialis']:
    for stage, group in df[df['species'] == species].groupby('stage'):
        idx = stage_to_idx.get(stage, None)
        if idx is None:
            color = 'grey'
        else:
            if species == 'Calanus hyperboreus':
                color = colors_hyp_stage[idx]
            elif species == 'Calanus glacialis':
                color = colors_gla_stage[idx]
            else:
                color = 'grey'
            marker = markers.get(stage, 'o')

        ax.scatter(
            group['prosome_major_axis_mm'],
            group['lipid_volume_mm3'],
            color=color,
            edgecolor='none',
            marker=marker,
            label=f"{stage} - {species}",
            alpha=0.8,
            s=30
        )

# Legend and labels

legend_handles = []
stage_keys = ['civstage', 'cvstage', 'female']
    
for stage in stage_keys:
    idx = stage_to_idx[stage]
    # Patch Hyperboreus 
    handle_hyp = Line2D(
        [0], [0],
        marker=markers[stage],
        markerfacecolor=colors_hyp_stage[idx],
        markeredgecolor='none',
        markersize=7,
        linestyle='None'
    )
    # Patch Glacialis
    handle_gla = Line2D(
        [0], [0],
        marker=markers[stage],
        markerfacecolor=colors_gla_stage[idx],
        markeredgecolor='none',
        markersize=7,
        linestyle='None'
    )
    legend_handles.append((handle_hyp, handle_gla))
    
fig.legend(
    legend_handles,
    stage_labels,
    handler_map={tuple: HandlerTuple(ndivide=None)},
    loc='center right',
    bbox_to_anchor=(1, 0.5),
    frameon=False,
)

fig.text(
    0.91,    
    0.56,    
    r'$C.\ hyperboreus$', 
    rotation=90, 
    fontsize=10, 
    va='bottom', 
    ha='left',
    style='italic'
)

fig.text(
    0.935,    
    0.56,    
    r'$C.\ glacialis$', 
    rotation=90, 
    fontsize=10, 
    va='bottom', 
    ha='left',
    style='italic'
)

ax.spines['top'].set_visible(False) 
ax.spines['right'].set_visible(False)
ax.set_xlim([2,8])
ax.set_ylim([0,10])
ax.set_xlabel('Prosome length (mm)')
ax.set_ylabel('Lipid sac volume (mm3)')

plt.savefig("LipidSacVolume_against_PL_CIV_CV_Females_August_LOKI2013.png", dpi=600)
plt.show()


#### Figure 5: Lipid sac volume from ellispe VS cylinder

df['lipid_volume_mm3_cylinder'] = (np.pi * df['lipid_area_mm2']**2) / (4 * df['lipid_major_axis_mm'])


plt.scatter(df['lipid_volume_mm3_cylinder'], df['lipid_volume_mm3'], edgecolor='black', facecolor='none')

min_val = min(df['lipid_volume_mm3_cylinder'].min(), df['lipid_volume_mm3'].min())
max_val = max(df['lipid_volume_mm3_cylinder'].max(), df['lipid_volume_mm3'].max())
lims = [min_val, max_val+0.5]

plt.plot(lims, lims, 'k--')

plt.xlabel("Lipid sac volume (mm3) - Cylinder")
plt.ylabel("Lipid sac volume (mm3) - Ellipse")

# Mettre un aspect égal
plt.gca().set_aspect('equal', adjustable='box')

# Limites égales sur x et y
plt.xlim(lims)
plt.ylim(lims)

plt.savefig("LipidSacVolume_from_cylinderVSellipse_CIV_CV_Females_August_LOKI2013.png", dpi=600)
plt.show()


