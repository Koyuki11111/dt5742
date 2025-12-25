#!/usr/bin/env python3

import argparse
import os
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def check_hdf5_qa(input_dir, output_dir, num_events=5):
    """
    Reads HDF5 files from input_dir and generates QA plots in output_dir.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Find HDF5 files
    h5_files = sorted(glob.glob(os.path.join(input_dir, "*.h5")))
    if not h5_files:
        print(f"No HDF5 files found in {input_dir}")
        return

    print(f"Found {len(h5_files)} HDF5 files: {[os.path.basename(f) for f in h5_files]}")

    for h5_file in h5_files:
        file_name = os.path.basename(h5_file)
        print(f"\nProcessing {file_name}...")
        
        try:
            with h5py.File(h5_file, 'r') as f:
                if 'Hits' not in f:
                    print(f"  WARNING: 'Hits' dataset not found in {file_name}. Skipping.")
                    continue

                data = f['Hits'][:]
                
                # Check if data is empty
                if len(data) == 0:
                    print(f"  WARNING: 'Hits' is empty in {file_name}. Skipping.")
                    continue

                # Extract fields
                # The structured array fields: event, channel, sensor_id, column_id, strip_id, ampMax, ...
                events = data['event']
                amp_max = data['ampMax']
                col_ids = data['column_id']
                strip_ids = data['strip_id']
                sensor_ids = data['sensor_id']
                
                unique_sensors = np.unique(sensor_ids)
                
                # 1. Summary Histogram of Amplitude
                plt.figure(figsize=(10, 6))
                plt.hist(amp_max, bins=100, range=(0, 4096), histtype='stepfilled', alpha=0.7, label='All Channels')
                plt.title(f"Amplitude Distribution - {file_name}")
                plt.xlabel("Amplitude (ADC)")
                plt.ylabel("Counts")
                plt.yscale('log')
                plt.grid(True, which="both", ls="-", alpha=0.5)
                plt.legend()
                
                summary_plot_name = f"summary_amp_{os.path.splitext(file_name)[0]}.png"
                plt.savefig(os.path.join(output_dir, summary_plot_name))
                plt.close()
                print(f"  Saved summary plot: {summary_plot_name}")

                # 2. Event-by-Event Amplitude Maps
                unique_events = np.unique(events)
                print(f"  Total events: {len(unique_events)}. Plotting first {num_events}...")

                for i, ev_id in enumerate(unique_events[:num_events]):
                    # Mask for current event
                    mask = (events == ev_id)
                    ev_data = data[mask]
                    
                    if len(ev_data) == 0:
                        continue

                    # Determine grid size for this event (and sensor)
                    # We assume one sensor per file usually, or split plots by sensor
                    for sid in unique_sensors:
                        sensor_mask = (ev_data['sensor_id'] == sid)
                        sensor_data = ev_data[sensor_mask]
                        
                        if len(sensor_data) == 0:
                            continue

                        s_cols = sensor_data['column_id']
                        s_strips = sensor_data['strip_id']
                        s_amps = sensor_data['ampMax']
                        
                        # Determine bounds
                        min_col, max_col = np.min(s_cols), np.max(s_cols)
                        min_strip, max_strip = np.min(s_strips), np.max(s_strips)
                        
                        # Create a dense grid for plotting (filled with 0 or NaN)
                        # We use a dictionary or directly scatter/hist2d
                        # hist2d is good for "maps"
                        
                        # Calculate bin edges to center the pixels
                        x_bins = np.arange(min_col - 0.5, max_col + 1.5, 1)
                        y_bins = np.arange(min_strip - 0.5, max_strip + 1.5, 1)
                        
                        plt.figure(figsize=(8, 6))
                        h = plt.hist2d(s_cols, s_strips, weights=s_amps, bins=[x_bins, y_bins], 
                                       cmin=0.1, cmap='viridis', vmin=0, vmax=300) # Fixed scale for comparison
                        plt.colorbar(h[3], label='Amplitude (ADC)')
                        
                        plt.title(f"Event {ev_id} - Sensor {sid}\n{file_name}")
                        plt.xlabel("Column ID")
                        plt.ylabel("Strip ID")
                        plt.xticks(np.arange(min_col, max_col + 1, 1))
                        plt.yticks(np.arange(min_strip, max_strip + 1, 1))
                        plt.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

                        # Annotate values
                        for sc, ss, sa in zip(s_cols, s_strips, s_amps):
                            plt.text(sc, ss, f"{sa:.0f}", ha='center', va='center', color='white', fontsize=8)

                        plot_name = f"event_{ev_id:06d}_sensor{sid}_{os.path.splitext(file_name)[0]}.png"
                        plt.savefig(os.path.join(output_dir, plot_name))
                        plt.close()
                        print(f"    Saved event map: {plot_name}")

        except Exception as e:
            print(f"  ERROR processing {file_name}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate QA plots from HDF5 exported files.")
    parser.add_argument("--input-dir", required=True, help="Directory containing HDF5 files")
    parser.add_argument("--output-dir", default="qa_plots", help="Directory to save plots")
    parser.add_argument("--num-events", type=int, default=5, help="Number of events to inspect")
    
    args = parser.parse_args()
    
    check_hdf5_qa(args.input_dir, args.output_dir, args.num_events)
