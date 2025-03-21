#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 13:43:42 2025

This plots the simulation data and data from PM 


@author: Gajendra Gurung
"""
import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime

### _____ HERE IS THE SIMULATION PLOT ____ ####

spectra = pd.read_csv("/Users/Gajju/Downloads/pm_sim_v2/bi207spectra.txt",
                      sep='\s+', header=None)

# Assign column names
spectra.columns = ["BinValue", "outer_long", "inner_long_electrons", 
                   "inner_long_gammas", 
                   "inner_long", "outer_short", "inner_short_electrons", 
                   "inner_short_gammas", "inner_short"]

# Plot for long drift spectra
plt.figure(figsize=(10, 6))
for col in ["outer_long", "inner_long_electrons", "inner_long_gammas", 
            "inner_long"]:
    plt.step(spectra["BinValue"], spectra[col], label=col)
plt.xlabel("BinValue")
plt.ylabel("Count")
plt.title("Long PM Spectra")
plt.legend()
plt.show()

# Plot for short drift spectra
plt.figure(figsize=(10, 6))
for col in ["outer_short", "inner_short_electrons", "inner_short_gammas", 
            "inner_short"]:
    plt.step(spectra["BinValue"], spectra[col], label=col)
plt.xlabel("BinValue")
plt.ylabel("Count")
plt.title("Short PM Spectra")
plt.legend()
plt.show()




### _____ HERE IS THE DATA PLOT ____ ####

# Map textual month to a numeric string
MONTH_MAP = {
    'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
    'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
    'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'
}

spectra_mapping = {
    "F1.txt": "inner_short",
    "F2.txt": "inner_long",
    "F3.txt": "outer_short",
    "F4.txt": "outer_long"
}


# Set the base path (always in Year_month format)
base_path = "/Users/Gajju/NP02_activities/purityMonitor/np02data/stable_data/2024_Mar"

# List all subdirectories (including base_path itself)
subdirs = [root for root, dirs, files in os.walk(base_path)]

# Initialize an empty DataFrame for all the data
exp_data = pd.DataFrame()

for d in subdirs:
    for file_name in ["F1.txt", "F2.txt", "F3.txt", "F4.txt"]:
        file_path = os.path.join(d, file_name)
        if os.path.exists(file_path):
            # Read the CSV file.
            data = pd.read_csv(file_path, sep=",", header=0)
            # 
            
            # Extract timestamp from directory structure.
            # Expected structure: .../<Year_Month>/<Day>/<Hour>/<Minute>/
            path_parts = d.split(os.sep)
            if len(path_parts) >= 4:
                year_month = path_parts[-4]  # e.g. "2024_Dec"
                day = path_parts[-3]
                hour = path_parts[-2]
                minute = path_parts[-1]
                
                try:
                    year_str, month_word = year_month.split('_')
                    # Use the mapping to get the numeric month, defaulting to 'NA' if not found.
                    month_str = MONTH_MAP.get(month_word, 'NA')
                    if month_str == 'NA':
                        timestamp = "NA"
                    else:
                        timestamp_str = f"{year_str}-{month_str}-{day} {hour}:{minute}"
                        timestamp = datetime.strptime(timestamp_str, '%Y-%m-%d %H:%M')
                except (ValueError, IndexError):
                    timestamp = "NA"
            else:
                timestamp = "NA"
            
            data['timestamp'] = timestamp
            
            # Add the spectra column based on the file name.
            data['spectra'] = spectra_mapping.get(file_name, "NA")
            
            # Append to the aggregated DataFrame.
            exp_data = pd.concat([exp_data, data], ignore_index=True)



inner_long_data = exp_data[exp_data['spectra'] == "inner_long"]
outer_long_data = exp_data[exp_data['spectra'] == "outer_long"]

# Extract unique measurement times and sort them for inner_long only.
unique_ts = inner_long_data['timestamp'].dropna().unique()

plt.figure(figsize=(10, 6))
for ts in sorted(unique_ts):
    subset = inner_long_data[inner_long_data['timestamp'] == ts]
    # Plot each individual histogram with a faded grey line.
    plt.step(subset['BinCenter'], subset['Population'], linestyle='-', color='pink', alpha=0.3)

# Compute the average histogram by grouping by BinCenter.
avg_inner_long_data = inner_long_data.groupby('BinCenter', as_index=False)['Population'].mean()

# Overlay the average histogram with a bold line.
plt.step(avg_inner_long_data['BinCenter'], avg_inner_long_data['Population'], linestyle='-', color='red', linewidth=1, label='Average')

plt.xlabel("Voltage [V]")
plt.ylabel("count")
plt.title("Inner anode Spectrum of long PM (Stable data Mar 12 & 13)")
plt.legend(loc='best', fontsize='small')
plt.grid(True)
plt.show()


plt.figure(figsize=(10, 6))
for ts in sorted(unique_ts):
    subset = outer_long_data[outer_long_data['timestamp'] == ts]
    # Plot each individual histogram with a faded grey line.
    plt.step(subset['BinCenter'], subset['Population'], linestyle='-', color='pink', alpha=0.3)

# Compute the average histogram by grouping by BinCenter.
avg_outer_long_data = outer_long_data.groupby('BinCenter', as_index=False)['Population'].mean()

# Overlay the average histogram with a bold line.
plt.step(avg_outer_long_data['BinCenter'], avg_outer_long_data['Population'], linestyle='-', color='red', linewidth=1, label='Average')

plt.xlabel("Voltage [V]")
plt.ylabel("count")
plt.title("Outer anode Spectrum of long PM (Stable data Mar 12 & 13)")
plt.legend(loc='best', fontsize='small')
plt.grid(True)
plt.show()


inner_short_data = exp_data[exp_data['spectra'] == "inner_short"]
plt.figure(figsize=(10, 6))
for ts in sorted(unique_ts):
    subset = inner_short_data[inner_short_data['timestamp'] == ts]
    # Plot each individual histogram with a faded grey line.
    plt.step(subset['BinCenter'], subset['Population'], linestyle='-', color='pink', alpha=0.3)

# Compute the average histogram by grouping by BinCenter.
avg_inner_short_data = inner_short_data.groupby('BinCenter', as_index=False)['Population'].mean()

# Overlay the average histogram with a bold line.
plt.step(avg_inner_short_data['BinCenter'], avg_inner_short_data['Population'], linestyle='-', color='red', linewidth=1, label='Average')

plt.xlabel("Voltage [V]")
plt.ylabel("count")
plt.title("Inner anode Spectrum of Short PM (Stable data Mar 12 & 13)")
plt.legend(loc='best', fontsize='small')
plt.grid(True)
plt.show()



### _____ HERE IS THE DATA AND SIMULATION OVERLAY ____ ####

# Define scaling factors for simulation data.
x_scale_factor = 0.737  # Scale factor for the x-axis (BinValue)
y_scale_factor = 1.755  # Scale factor for the y-axis (inner_long)

plt.figure(figsize=(10, 6))

plt.step(
    spectra["BinValue"] * x_scale_factor, 
    spectra["inner_long"] * y_scale_factor, 
    linestyle='-', 
    color='red', 
    linewidth=1.5, 
    label="Inner long Simulation"
)

plt.step(
    spectra["BinValue"] * 0.8, 
    spectra["outer_long"] * 1.75, 
    linestyle='-', 
    color='blue', 
    linewidth=1.5, 
    label="Outer Long Simulation"
)

plt.step(
    avg_inner_long_data['BinCenter'], 
    avg_inner_long_data['Population'], 
    linestyle='-', 
    color='green', 
    linewidth=1.5, 
    label="Inner long Data"
)

plt.step(
    avg_outer_long_data['BinCenter'], 
    avg_outer_long_data['Population'], 
    linestyle='-', 
    color='darkorange', 
    linewidth=1.5, 
    label="Outer Long Data"
)

plt.xlabel("Voltage [V]")
plt.ylabel("Count")
plt.title("Overlayed Scaled Simulation vs Data for Inner anode for Long PM")
plt.legend(loc='best', fontsize='small')
plt.grid(True)
plt.show()



