#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Link Ecotaxa classification with Lipids and Prosome masks.
Here lipid sac polygons were done by one annotator while prosome polygons were 
done by two other annotators on Label Studio. 
The lipid sac and prosome polygons are then of different types.

Only Calanus in lateral position are keep here.

@author: LucieBourreau
@date: 2025/02/12
"""

import pandas as pd

#### Load datasets

classification = pd.read_csv("./ecotaxa_export_2331_20250212_1631.tsv",
                             sep="\t")

lipid_masks = pd.read_csv("./lipid_segmentation.csv", sep=";")

# Prosome masks were done by different annotator so we have more than one file
prosome_masks_1 = pd.read_csv("./prosome_segmentation.csv")
prosome_masks_2 = pd.read_csv("./project-14-at-2024-10-18-16-12-341c8e96.csv")
prosome_masks_3 = pd.read_csv("./project-10-at-2025-02-24-13-03-e3aca7e5.csv")

# Merge the 3 prosome masks df
prosome_masks = pd.concat([prosome_masks_1, prosome_masks_2, prosome_masks_3], ignore_index=True)

#### Isolate ID of calanus + lateral +  C4 - C5 - Females

filtered_ids = classification[classification["object_annotation_category"].str.contains("calanus", case=False, na=False) & 
                              classification["object_annotation_category"].str.contains("lateral", case=False, na=False) &
                              -classification["object_annotation_category"].str.contains("ciiistage", case=False, na=False)]["object_id"]

filtered_ids = filtered_ids.tolist()

#### Create a merged df with calanus + lateral classification and associated
#### lipid and prosome masks

def extract_object_id(image_path):
    """
    Replace the image name from label studio to the same object id as in Ecotaxa
    """
    image_name = image_path.split("/")[-1]  # Keep only file name
    image_name = image_name.replace(".bmp", "")  # Remove ".bmp"
    image_name = image_name.replace(".jpg", "")  # Remove ".jpg"
    
    if "-" in image_name:
        image_name = image_name.split("-")[-1]  # Remove everything before '-'
    
    image_name = image_name.replace("_", " ")  # Replace "_" by " "
    
    return image_name

# Add 'object_id' column
lipid_masks["object_id"] = lipid_masks["filename"].str.replace(".bmp", "", regex=False)
prosome_masks["object_id"] = prosome_masks["image"].apply(extract_object_id)

# Filter classification
filtered_classification = classification[classification["object_id"].isin(filtered_ids)]

# Construct the dataset with selected columns
df_merged = (
    filtered_classification[["object_id", "object_date", "object_annotation_category"]]  
    .merge(lipid_masks[["segmentation", "shape", "shape_y", "shape_x", "object_id"]], on="object_id", how="left") 
    .merge(prosome_masks[["label", "object_id"]], on="object_id", how="left")  
)

# Rename columns
df_final = df_merged.rename(columns={
    "segmentation": "lipid_polygon",
    "shape": "shape_lipid_polygon",
    "shape_x": "shape_x_lipid_polygon",
    "shape_y": "shape_y_lipid_polygon",
    "label": "prosome_polygon"
})

#### Export dataset
df_final.to_csv("./merged_LOKI2013_ecotaxa_masks.csv", index=False)
