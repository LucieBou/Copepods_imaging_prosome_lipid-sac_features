#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create the dataset for Coltrane calibration.
Simply to keep only the necessary variables in the dataframe.

@author: LucieBourreau
@date: 2025/03/12
"""

import pandas as pd

df = pd.read_csv("merged_LOKI2013_ecotaxa_masks_features.csv")

to_keep = ['object_annotation_category', 'total_lipids_ugC', 'fullness_ratio_carbon_volume']

df_final = df[to_keep]

df_final.to_csv("merged_LOKI2013_ecotaxa_masks_features_for_calibration.csv", index=False)