#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot all prosome and lipid masks with major and minor axis

@author: LucieBourreau
@date: 2025/03/10
"""


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import json
from skimage.draw import polygon
import os
import ast

def polygon_to_mask(polygon_points, image_width, image_height):
    """
    Create a mask from polygon points

    Parameters
    ----------
    polygon_points : list
        in pixel.
    image_width : float
        in pixel.
    image_height : float
        in pixel.

    Returns
    -------
    mask : np.array
        associated mask (bool).

    """
    
    # Zero mask from the image
    mask = np.zeros((image_height, image_width), dtype=np.uint8)
    
    # Put to 1 pixel within the polygon
    rr, cc = polygon([point[1] for point in polygon_points], [point[0] for point in polygon_points], shape=mask.shape)
    
    mask[rr, cc] = 1
    
    return mask


def create_transparent_colormap(base_cmap):
    cmap = plt.cm.get_cmap(base_cmap)  # base map color map
    cmap_colors = cmap(np.arange(cmap.N))  # Isolate the colors
    cmap_colors[0, -1] = 0  # Put the first color transparent
    return ListedColormap(cmap_colors)

transparent_reds = create_transparent_colormap("Reds")
transparent_blues = create_transparent_colormap("Blues")


### Plot and save all msks with major and minor axis into a folder

output_dir = "./Masks_with_major_minor_axis"
os.makedirs(output_dir, exist_ok=True)  # Create folder if it does not exists

data = pd.read_csv("./merged_LOKI2013_ecotaxa_masks_features.csv")

pixel_size_mm = 0.023

for index, row in data.head(10).iterrows():
    object_id = row["object_id"]
    
    # Polygon data pour prosome et lipid
    prosome_polygon_points = row['prosome_polygon_points_px']
    prosome_polygon_width = json.loads(row['prosome_polygon'])[0]["original_width"]
    prosome_polygon_height = json.loads(row['prosome_polygon'])[0]["original_height"]
    

    lipid_polygon_points = json.loads(row['lipid_polygon'])[0]
    lipid_polygon_width = row["shape_x_lipid_polygon"]
    lipid_polygon_height = row["shape_y_lipid_polygon"]

    # Compute the mask for the figure
    prosome_mask = polygon_to_mask(ast.literal_eval(prosome_polygon_points), prosome_polygon_width, prosome_polygon_height)
    lipid_mask = polygon_to_mask(lipid_polygon_points, lipid_polygon_width, lipid_polygon_height)
    
    # Found centroid of each mask
    prosome_centroid = row['prosome_centroid']
    lipid_centroid = row['lipid_centroid']
    
    # Create the figure
    fig, ax = plt.subplots()
    ax.imshow(prosome_mask, cmap=transparent_reds, alpha=0.4)
    ax.imshow(lipid_mask, cmap=transparent_blues, alpha=0.4)

    prosome_major_end_1 = ast.literal_eval(row['prosome_major_end_1']) if isinstance(row['prosome_major_end_1'], str) else row['prosome_major_end_1']
    prosome_major_end_2 = ast.literal_eval(row['prosome_major_end_2']) if isinstance(row['prosome_major_end_2'], str) else row['prosome_major_end_2']
    prosome_minor_end_1 = ast.literal_eval(row['prosome_minor_end_1']) if isinstance(row['prosome_minor_end_1'], str) else row['prosome_minor_end_1']
    prosome_minor_end_2 = ast.literal_eval(row['prosome_minor_end_2']) if isinstance(row['prosome_minor_end_2'], str) else row['prosome_minor_end_2']

    lipid_major_end_1 = ast.literal_eval(row['lipid_major_end_1']) if isinstance(row['lipid_major_end_1'], str) else row['lipid_major_end_1']
    lipid_major_end_2 = ast.literal_eval(row['lipid_major_end_2']) if isinstance(row['lipid_major_end_2'], str) else row['lipid_major_end_2']
    lipid_minor_end_1 = ast.literal_eval(row['lipid_minor_end_1']) if isinstance(row['lipid_minor_end_1'], str) else row['lipid_minor_end_1']
    lipid_minor_end_2 = ast.literal_eval(row['lipid_minor_end_2']) if isinstance(row['lipid_minor_end_2'], str) else row['lipid_minor_end_2']


    # Plot the axis
    ax.plot([prosome_major_end_1[0], prosome_major_end_2[0]], 
            [prosome_major_end_1[1], prosome_major_end_2[1]], 
            '--r', 
            linewidth=1)  # Prosome major axis
    ax.plot([prosome_minor_end_1[0], prosome_minor_end_2[0]], 
            [prosome_minor_end_1[1], prosome_minor_end_2[1]], 
            '--r', 
            linewidth=1)  # Prosome minor axis
    
    ax.plot([lipid_major_end_1[0], lipid_major_end_2[0]], 
            [lipid_major_end_1[1], lipid_major_end_2[1]], 
            '--b', 
            linewidth=1)  # Prosome major axis
    ax.plot([lipid_minor_end_1[0], lipid_minor_end_2[0]], 
            [lipid_minor_end_1[1], lipid_minor_end_2[1]], 
            '--b', 
            linewidth=1)  # Prosome minor axis

    
    ax.set_title(f"id: {object_id}")
    
    # Save plot
    file_path = os.path.join(output_dir, f"{object_id}.png")
    plt.savefig(file_path, dpi=300)
    plt.close(fig)


