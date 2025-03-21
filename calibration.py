#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 12:04:29 2025

For Calibration with Test pulse 

@author: Gajju
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.dates as mdates

# Define the Gaussian function.
def gaussian(x, A, mu, sigma):
    return A * np.exp(- (x - mu)**2 / (2 * sigma**2))

# Map textual month to a numeric string.
MONTH_MAP = {
    'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
    'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
    'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'
}

# Mapping for file names to spectra names.
spectra_mapping = {
    "F1.txt": "inner_short",
    "F2.txt": "inner_long",
    "F3.txt": "outer_short",
    "F4.txt": "outer_long"
}

# Set the base path.
base_path = "/Users/Gajju/NP02_activities/purityMonitor/np02data/calib_Mar20"

# List all subdirectories (including base_path itself).
subdirs = [root for root, dirs, files in os.walk(base_path)]

# Initialize an empty DataFrame for all the data.
exp_data = pd.DataFrame()

for d in subdirs:
    for file_name in ["F1.txt", "F2.txt", "F3.txt", "F4.txt"]:
        file_path = os.path.join(d, file_name)
        if os.path.exists(file_path):
            data = pd.read_csv(file_path, sep=",", header=0)
            # Extract timestamp from directory structure.
            path_parts = d.split(os.sep)
            if len(path_parts) >= 4:
                year_month = path_parts[-4]  # e.g. "2024_Dec"
                day = path_parts[-3]
                hour = path_parts[-2]
                minute = path_parts[-1]
                try:
                    year_str, month_word = year_month.split('_')
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
            # Add the spectra column.
            data['spectra'] = spectra_mapping.get(file_name, "NA")
            exp_data = pd.concat([exp_data, data], ignore_index=True)

# Define spectra types.
spectra_types = ["inner_short", "inner_long", "outer_short", "outer_long"]

# Dictionary to store Gaussian fit results (if needed for later).
fit_results = {spec: [] for spec in spectra_types}

# Create a 2x2 grid of subplots.
fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)

x_min, x_max = 1.325, 1.475

for ax, spec in zip(axes.flatten(), spectra_types):
    spec_data = exp_data[exp_data['spectra'] == spec]
    unique_ts = spec_data['timestamp'].dropna().unique()
    sorted_ts = sorted(unique_ts)
    cmap = plt.get_cmap("tab10")
    
    # If there are timestamps, use the date from the first one for the title.
    if sorted_ts and isinstance(sorted_ts[0], datetime):
        date_str = sorted_ts[0].strftime("%Y-%m-%d")
    else:
        date_str = ""
        
    # Loop over each timestamp.
    for i, ts in enumerate(sorted_ts):
        subset = spec_data[spec_data['timestamp'] == ts]
        color = cmap(i % cmap.N)
        # Use only time for label.
        time_label = ts.strftime("%H:%M") if isinstance(ts, datetime) else str(ts)
        
        # Plot the histogram (faded).
        ax.step(subset['BinCenter'], subset['Population'], linestyle='-', color=color, alpha=0.3)
        
        # Prepare data for Gaussian fit.
        xdata = subset['BinCenter'].values
        ydata = subset['Population'].values
        if len(xdata) < 3:
            continue
        
        # Improved initial guesses and bounds:
        A0 = np.max(ydata)
        mu0 = np.sum(xdata * ydata) / np.sum(ydata)
        sigma0 = np.sqrt(np.sum(ydata * (xdata - mu0)**2) / np.sum(ydata))
        try:
            popt, _ = curve_fit(gaussian, xdata, ydata, p0=[A0, mu0, sigma0],
                                bounds=([0, x_min, 0], [np.inf, x_max, np.inf]))
            A_fit, mu_fit, sigma_fit = popt
            fit_results[spec].append((ts, mu_fit, sigma_fit))
            x_fit = np.linspace(x_min, x_max, 500)
            y_fit = gaussian(x_fit, A_fit, mu_fit, sigma_fit)
            # Compute chi²/ndf.
            y_fit_data = gaussian(xdata, *popt)
            
            # Label: only time and chi²/ndf.
            gf_label = f"{time_label}"
            ax.plot(x_fit, y_fit, color=color, linewidth=1, label=gf_label)
        except Exception as e:
            print(f"Gaussian fit failed for {spec} at {ts}: {e}")
    
    ax.set_title(f"{spec} - {date_str}")
    ax.legend(loc='best', fontsize='small')
    

plt.xlim(x_min, x_max)
fig.text(0.5, 0.04, 'Voltage [V]', ha='center', fontsize=12)
fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical', fontsize=12)
fig.suptitle("Response for 3 mV Test Pulse for Calibration", fontsize=14)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()


fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

for spec in spectra_types:
    # Sort results by timestamp.
    results = sorted(fit_results[spec], key=lambda x: x[0])
    if not results:
        continue
    times = [r[0] for r in results]
    mus = [r[1] for r in results]
    sigmas = [r[2] for r in results]
    
    ax1.plot(times, mus, marker='o', linestyle='None', label=spec)
    ax2.plot(times, sigmas, marker='o', linestyle='None', label=spec)

ax1.set_ylabel("Mean [V]")
ax1.set_title("Gaussian Fit Mean vs Time (Test pulse = 3 mV)")
ax1.legend(loc='best', fontsize='small')

ax2.set_ylabel("Sigma [V]")
ax2.set_title("Gaussian Fit Sigma vs Time (Test pulse = 3 mV)")
ax2.legend(loc='best', fontsize='small')

# Format x-axis of ax2 with date labels.
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
for label in ax2.get_xticklabels():
    label.set_rotation(45)

ax2.set_xlabel("Time")
plt.tight_layout()
plt.show()


# Second figure: Plot mean with sigma as error bars.
fig3, ax3 = plt.subplots(figsize=(12, 8))

for spec in spectra_types:
    results = sorted(fit_results[spec], key=lambda x: x[0])
    if not results:
        continue
    times = [r[0] for r in results]
    mus = [r[1] for r in results]
    sigmas = [r[2] for r in results]
    
    ax3.errorbar(times, mus, yerr=sigmas, marker='o', linestyle='None', capsize=2, label=spec)

ax3.set_xlabel("Time")
ax3.set_ylabel("Mean [V]")
ax3.set_title("Gaussian Fit Mean vs Time with Sigma as Error Bars")
ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
for label in ax3.get_xticklabels():
    label.set_rotation(45)

ax3.legend(loc='best', fontsize='small')
plt.tight_layout()
plt.show()


fig, ax = plt.subplots(figsize=(12, 8))

for spec in spectra_types:
    # Sort the fit results by timestamp.
    results = sorted(fit_results[spec], key=lambda r: r[0])
    if not results:
        continue
    times = [r[0] for r in results]
    # Compute sigma/mean ratio for each timestamp.
    ratios = [2.35482*r[2] / r[1] if r[1] != 0 else np.nan for r in results]
    avg_ratio = np.nanmean(ratios)

# Create a label that includes the spectra name and its mean sigma/mean ratio.
    label_text = f"{spec} (mean: {avg_ratio:.3f})"
    
    ax.plot(times, ratios, marker='o', linestyle='None', label=label_text)

ax.set_xlabel("Time")
ax.set_ylabel("FWHM / Mean")
ax.set_title("FWHM / Mean Ratio vs Time")
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
for label in ax.get_xticklabels():
    label.set_rotation(45)
ax.legend(loc='best', fontsize='small')
plt.tight_layout()
plt.show()
