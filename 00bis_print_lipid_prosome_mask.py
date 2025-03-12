#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the prosome and lipid masks on the original image and save them in a folder

@author: LucieBourreau
@date: 2025/02/17
"""

from PIL import Image, ImageDraw, ImageFont
import ast
import pathlib
import numpy as np
import pandas as pd
import json

def draw_masks_from_polygon_data(data, root_folderpath: pathlib.Path, original_images_folder: pathlib.Path):
    """
    Print original LOKI images with prosome and lipid masks.

    Parameters
    ----------
    data : pd.DataFrame
        Data from merged prosome and lipid masks and Ecotataxa classif.
    root_folderpath : pathlib.Path
        Folder path where to store images with masks.
    original_images_folder : pathlib.Path
        Folder path to original LOKI images.

    Returns
    -------
    None.

    """
    
    #### LIPIDS
    
    # Create a black mask for lipids 
    lipid_width = data["shape_x_lipid_polygon"]
    lipid_height = data["shape_y_lipid_polygon"]
    
    # Extract lipid points from polygon
    lipid_scaled_points = [
        [int(x), int(y)] for x, y in ast.literal_eval(data["lipid_polygon"])[0]
    ]
    
    #  Convert to integer coordinates for PIL
    lipid_polygon_points = [(x, y) for x, y in lipid_scaled_points]
    
    # Create lipid mask (boolean)
    lipid_mask_image = Image.new("RGBA", (lipid_width, lipid_height), (0, 0, 0, 0))  
    lipid_draw = ImageDraw.Draw(lipid_mask_image)
    lipid_draw.polygon(lipid_polygon_points, outline=(0, 0, 255, 90), fill=(0, 0, 255, 90)) # transparent blue

    #### PROSOME
    prosome_mask = json.loads(data["prosome_polygon"])
    
    prosome_width = prosome_mask[0]["original_width"]
    prosome_height = prosome_mask[0]["original_height"]
    
    # Necessary because Label Studio gives polygon coordinates in pixel percentages rather than pixels.
    prosome_scaled_points = [
        [int(x * prosome_width / 100), int(y * prosome_height / 100)] for x, y in prosome_mask[0]["points"]
    ]
    
    prosome_polygon_points = [(x, y) for x, y in prosome_scaled_points]
    
    prosome_mask_image = Image.new("RGBA", (prosome_width, prosome_height), (0, 0, 0, 0))
    prosome_draw = ImageDraw.Draw(prosome_mask_image)
    prosome_draw.polygon(prosome_polygon_points, outline=(255, 0, 0, 90), fill=(255, 0, 0, 90))  # transparent red

    #### IMAGE

    # Load original LOKI image
    original_images_folder = pathlib.Path(original_images_folder)
    original_image_path = original_images_folder / f"{data['object_id']}.jpg"
    original_image = Image.open(original_image_path)

    # Ensure that the size of the original image is the same as that of the mask
    original_image = original_image.resize((lipid_width, lipid_height))

    # Convert original image to RGBA for transparent overlay
    original_image = original_image.convert("RGBA")

    # Overlay masks on the original image with transparency
    original_image.paste(prosome_mask_image, (0, 0), prosome_mask_image)       
    original_image.paste(lipid_mask_image, (0, 0), lipid_mask_image)
    
    # Add the title (object_id)
    draw = ImageDraw.Draw(original_image)
    font = ImageFont.load_default(size=8 ) 
    
    # Add text on top of the image
    title = data['object_annotation_category']
    draw.text((10, 10), title, fill="white", font=font)

    #### SAVE

    # Save image with masks
    root_folderpath = pathlib.Path(root_folderpath)
    output_path = root_folderpath / f"{data['object_id'].replace(' ', '_')}_masks.png"
    original_image.save(output_path)


df = pd.read_csv("./merged_LOKI2013_ecotaxa_masks.csv")


for i in range(0,len(df.head(10))):
    draw_masks_from_polygon_data(df.iloc[i], 
                                 root_folderpath="./Images_with_masks",
                                 original_images_folder="/Users/luciebourreau/Library/CloudStorage/OneDrive-UniversitéLaval/PhD_ULaval/Data/LOKI/Export_Ecotaxa_ValidatedAll_LOKI2013/export_2331_20250212_1631")
