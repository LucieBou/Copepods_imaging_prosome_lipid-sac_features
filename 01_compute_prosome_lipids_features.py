#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute features (prosome length, lipids to carbon, fullness etc.) from polygons. 
Those features are necessary for Coltrane calibration (https://github.com/LucieBou/Coltrane_calibration.git).

@author: LucieBourreau
@date: 2025/02/17
"""

import pandas as pd
import numpy as np
import json
from itertools import combinations

def modify_prosome_perc_px_into_px(prosome_polygon):
    """
    Modify prosome polygon points that are in percentage of pixel when we extract
    them from Label Studio to be in pixel.

    Parameters
    ----------
    prosome_polygon : list
        polygon points in percentage of pixel.

    Returns
    -------
    prosome_polygon_points_px : list
        polygon points in pixel.

    """
    
    prosome_polygon_points = json.loads(prosome_polygon)[0]["points"]
    prosome_polygon_width = json.loads(prosome_polygon)[0]["original_width"]
    prosome_polygon_height = json.loads(prosome_polygon)[0]["original_height"]

    prosome_polygon_points_px = [
        [int(x * prosome_polygon_width / 100), int(y * prosome_polygon_height / 100)] for x, y in prosome_polygon_points
    ]
    
    return prosome_polygon_points_px

        
def polygon_area_mm2(polygon_points, pixel_size_mm):
    """
    Compute polygon area in mm2 from px2

    Parameters
    ----------
    polygon_points : LIST
        List containing list of coordinates in px. eg: [[x1,y1], [x2,y2], ...]

    Returns
    -------
    polygon_area_mm2 : LIST
        polygon area in mm2.

    """
    
    coords = np.array(polygon_points)
    
    x = coords[:, 0]
    y = coords[:, 1]

    polygon_area_px2 = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))) # Area in px2
    polygon_area_mm2 = polygon_area_px2 * pixel_size_mm**2  # Convert into mm2
    return polygon_area_mm2
    

def lipid_mm2_to_mg(lipid_area):
    """
    Convert lipid area in mm2 to mg of lipids (Vogedes et al., 2010)

    Parameters
    ----------
    lipid_area : FLOAT
        Lipid sac area in mm2.

    Returns
    -------
    total_lipids_mg : FLOAT
        Total lipid content in mg.

    """
    total_lipids_mg = 0.167 * lipid_area**(1.42)
    return total_lipids_mg


def major_axis(points, pixel_size_mm):
    """
    Compute major axis of a polygon and center it on the centroid.

    Parameters
    ----------
    points : list
        Polygon points in pixel.
    pixel_size_mm : float
        Pixel size in mm.

    Returns
    -------
    major_axis_endpoints : tuple
        Coordinates of the centered major axis.
    max_dist_mm : float
        Major axis length in mm.
    centroid: tuple
        Coordinates of the centroid.

    """
    max_dist = 0
    major_axis_endpoints = None

    for p1, p2 in combinations(points, 2):
        p1, p2 = np.array(p1), np.array(p2)
        dist = np.linalg.norm(p1 - p2)
        if dist > max_dist:
            max_dist = dist
            major_axis_endpoints = (tuple(p1), tuple(p2))
            
    max_dist_mm = max_dist*pixel_size_mm
    
    # Compute the centroid
    points_array = np.array(points)
    centroid = np.mean(points_array, axis=0)

    # Direction vector of major axis
    v = np.array(major_axis_endpoints[1]) - np.array(major_axis_endpoints[0])
    v_unit = v / np.linalg.norm(v)  # Normalize the vector
    
    # Compute new coordinates centered on the centroid of the polygon
    new_p1 = centroid - (v_unit * max_dist / 2)
    new_p2 = centroid + (v_unit * max_dist / 2)    
    
    return tuple(new_p1), tuple(new_p2), max_dist_mm, tuple(centroid)


def minor_axis_from_area(area, major_axis_length_mm, major_axis_endpoints_1, major_axis_endpoints_2):
    """
    Compute the minor axis of a polygon by considering the area of the polygon and its 
    major axis using the equation for Ellipse area. 
    Then compute the associated coordinates based on the center of the major axis 
    and a 90° angle. This last part is not accurate as the minor axis is not necessarily 
    in the middle of the major axis but it permit to have easy coordinates value for 
    visual validation. 

    Parameters
    ----------
    area : float
        Polygon area in mm2.
    major_axis_length_mm : float
        Major axis length in mm.
    major_axis_endpoints_1 : tuple
        First coordinates of the major axis.
    major_axis_endpoints_2 : tuple
        Second coordinates of the major axis.

    Returns
    -------
    minor_axis_length_mm : float
        Minor axis length in mm.
    minor_axis_endpoint_1: tuple
        First coordinates of the minor axis.
    minor_axis_endpoint_2: tuple
        Second coordinates of the minor axis.

    """
    minor_axis_length_mm = (4 * area) / (np.pi * major_axis_length_mm)
    
    (x1, y1) = major_axis_endpoints_1
    (x2, y2) = major_axis_endpoints_2
    
    center_x = (x1 + x2) / 2
    center_y = (y1 + y2) / 2
    
    # Compute major axis angle
    dx = x2 - x1
    dy = y2 - y1
    major_axis_angle = np.arctan2(dy, dx)
    
    # Minor axis angle is 90° of major axis angle
    minor_axis_angle = major_axis_angle + np.pi / 2
    
    half_length_minor = (minor_axis_length_mm/0.023) / 2 # in pixel
    
    x3 = center_x + half_length_minor * np.cos(minor_axis_angle)
    y3 = center_y + half_length_minor * np.sin(minor_axis_angle)
    
    x4 = center_x - half_length_minor * np.cos(minor_axis_angle)
    y4 = center_y - half_length_minor * np.sin(minor_axis_angle)

    return minor_axis_length_mm, (x3, y3), (x4, y4)


def polygon_volume(lipid_major_axis_length_mm, lipid_minor_axis_length_mm):
    """
    Compute polygon volume based on major and minor axis by using an Ellipse
    volume equation.

    Parameters
    ----------
    lipid_major_axis_length_mm : float
        Major axis length in mm.
    lipid_minor_axis_length_mm : float
        Minor axis length in mm.

    Returns
    -------
    volume : float
        Ellipse volume in mm3.

    """
    
    volume = (4/3) * np.pi * (lipid_major_axis_length_mm/2) * (lipid_minor_axis_length_mm/2)**2
    
    return volume



#### Load dataset
data = pd.read_csv("./merged_LOKI2013_ecotaxa_masks.csv")

pixel_size_mm = 0.023

#### Compute all features

## Convert prosome polygon points into pixel (rather than percentage of pixel originaly)
data['prosome_polygon_points_px'] = data['prosome_polygon'].apply(lambda x: modify_prosome_perc_px_into_px(x))

## lipid area
data['lipid_area_mm2'] = data['lipid_polygon'].apply(lambda x: polygon_area_mm2(json.loads(x)[0], pixel_size_mm))

## prosome area
data['prosome_area_mm2'] = data['prosome_polygon_points_px'].apply(lambda x: polygon_area_mm2(x, pixel_size_mm))

## lipid mg
data['total_lipids_mg'] = data['lipid_area_mm2'].apply(lambda x: lipid_mm2_to_mg(x))

## lipid ugC
data['total_lipids_ugC'] = data['total_lipids_mg'] * 0.79 * 1000 # Tarling et al., 2022

## Fullness ratio area
data['fullness_ratio_area'] = data['lipid_area_mm2'] / data['prosome_area_mm2']

## Fullness ratio carbon area
data['total_lipids_carbon_area'] = 0.79 * data['lipid_area_mm2'] # Tarling et al., 2022
data['prosome_carbon_area'] = (data["prosome_area_mm2"] - data['lipid_area_mm2']) * 0.2 * 0.45 # Ikeda & Skjoldal, 1989
data['fullness_ratio_carbon_area'] = data['total_lipids_carbon_area'] / (data['prosome_carbon_area'] + data['total_lipids_carbon_area'])

## Major and minor axis of lipid and prosome mask to compute their volume
data[['prosome_major_end_1', 'prosome_major_end_2', 'prosome_major_axis_mm', 'prosome_centroid']] = data.apply(
    lambda row: pd.Series(major_axis(row['prosome_polygon_points_px'], 0.023)),
    axis=1
)

data[['prosome_minor_axis_mm', 'prosome_minor_end_1', 'prosome_minor_end_2']] = data.apply(
    lambda row: pd.Series(minor_axis_from_area(
        row["prosome_area_mm2"],
        row['prosome_major_axis_mm'],
        row['prosome_major_end_1'],
        row['prosome_major_end_2']
    )),
    axis=1
)

data[['lipid_major_end_1', 'lipid_major_end_2', 'lipid_major_axis_mm', 'lipid_centroid']] = data.apply(
    lambda row: pd.Series(major_axis(json.loads(row['lipid_polygon'])[0], 0.023)),
    axis=1
)

data[['lipid_minor_axis_mm', 'lipid_minor_end_1', 'lipid_minor_end_2']] = data.apply(
    lambda row: pd.Series(minor_axis_from_area(
        row["lipid_area_mm2"],
        row['lipid_major_axis_mm'],
        row['lipid_major_end_1'],
        row['lipid_major_end_2']
    )),
    axis=1
)

## Compute prosome and lipid volume (Vogedes et al., 2010)
data['prosome_volume_mm3'] = data.apply(
    lambda row: polygon_volume(
        row['prosome_major_axis_mm'],
        row['prosome_minor_axis_mm'],
    ), 
    axis=1
)

data['lipid_volume_mm3'] = data.apply(
    lambda row: polygon_volume(
        row['lipid_major_axis_mm'],
        row['lipid_minor_axis_mm'],
    ), 
    axis=1
)

## Fullness ratio carbon volume
data['total_lipids_carbon_volume'] = 0.79 * data['lipid_volume_mm3'] # Tarling et al., 2022
data['prosome_carbon_volume'] = (data["prosome_volume_mm3"] - data['lipid_volume_mm3']) * 0.2 * 0.45 # Ikeda & Skjoldal., 1989
data['fullness_ratio_carbon_volume'] = data['total_lipids_carbon_volume'] / (data['prosome_carbon_volume'] + data['total_lipids_carbon_volume'])

## Remove ind from profiles 2013-08-16	2013-08-18 (only 13 individuals)
dates_to_remove = ["2013-08-16", "2013-08-18"]
data_filtered = data[~data["object_date"].isin(dates_to_remove)]

#### Save outputs
data.to_csv("merged_LOKI2013_ecotaxa_masks_features.csv", index=False)
data_filtered.to_csv("merged_LOKI2013_ecotaxa_masks_no2profiles_features.csv", index=False)

