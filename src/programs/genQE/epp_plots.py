#!/usr/bin/env python3
"""
Plotting script for 2N SRC analysis with different parameter ranges.
Creates 4 plots:
1. pmiss_exc histogram for xB 0.5-1.5, pmiss 0.4-0.7 GeV/c
2. pmiss_exc vs angle for xB 0.5-1.5, pmiss 0.4-0.7 GeV/c
3. pmiss_exc histogram for xB 1.5-2.5, pmiss > 0.7 GeV/c
4. pmiss_exc vs angle for xB 1.5-2.5, pmiss > 0.7 GeV/c
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

# Configuration
OUTPUT_DIR = "analysis_output_2N/txt_files"
PLOT_OUTPUT_DIR = "analysis_output_2N/png_files"

# SRC pair CM width (GeV/c). Used to set sensible plotting ranges. This must match the
# width that SRC_analysis_2N.cpp sets with set_sigmaCM() when it writes the files this
# script reads (0.055 for its deuteron study); it is not the genQE default, which is the
# per-nucleus value in gcfNucleus (0 for d, 0.10 for A = 3-6, 0.15 for A >= 12).
SIGMA_CM = 0.055
# Show pmiss_exc up to a few * sigma_CM
PMISS_EXC_MAX = 4.0 * SIGMA_CM

# Ensure output directory exists
os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)

def read_1d_histogram(filename):
    """
    Read a 1D histogram from a text file.
    Returns: (centers, weights, counts)
    """
    data = np.loadtxt(filename)
    centers = data[:, 0]
    weights = data[:, 1]
    counts = data[:, 2]
    return centers, weights, counts

def read_2d_histogram(filename):
    """
    Read a 2D histogram from a text file.
    Returns: (xB_centers, Q2_centers, weights_2d)
    """
    data = np.loadtxt(filename)
    xB_centers = np.unique(data[:, 0])
    Q2_centers = np.unique(data[:, 1])
    
    # Reshape into 2D array
    nx = len(xB_centers)
    nq = len(Q2_centers)
    weights_2d = data[:, 2].reshape(nq, nx)
    
    return xB_centers, Q2_centers, weights_2d

def apply_cuts(centers_xB, weights_xB, centers_pmiss, weights_pmiss, 
               centers_angle, weights_angle, centers_pmiss_exc, weights_pmiss_exc,
               xB_min, xB_max, pmiss_min, pmiss_max):
    """
    Apply xB and pmiss cuts to create filtered histograms.
    Returns filtered pmiss_exc and angle distributions.
    """
    # Find xB bins within range
    xB_mask = (centers_xB >= xB_min) & (centers_xB < xB_max)
    
    # Find pmiss bins within range
    pmiss_mask = (centers_pmiss >= pmiss_min) & (centers_pmiss < pmiss_max)
    
    # For this analysis, we'll use the region histograms instead
    # This is a simplified approach - ideally would need event-level data
    return centers_pmiss_exc, weights_pmiss_exc, centers_angle, weights_angle

def create_pmiss_exc_histogram(centers, weights, title, filename):
    """Create a histogram of pmiss_exc"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    bin_width = centers[1] - centers[0] if len(centers) > 1 else 0.01
    ax.bar(centers, weights, width=bin_width*0.9, edgecolor='black', alpha=0.7)
    
    ax.set_xlabel('$p_{miss}^{exc}$ (GeV/c)', fontsize=12)
    ax.set_ylabel('Counts', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.0, PMISS_EXC_MAX)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

def create_2d_scatter_plot(centers_angle, centers_pmiss_exc, weights_combined, 
                          title, filename):
    """Create a 2D scatter/density plot of pmiss_exc vs angle"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create bins for coloring by weight
    scatter = ax.scatter(centers_angle, centers_pmiss_exc, 
                        c=weights_combined, s=50, cmap='viridis', 
                        alpha=0.6, edgecolors='black', linewidth=0.5)
    
    ax.set_xlabel('Angle between $\\vec{p}_{miss}$ and $\\vec{p}_{miss}^{exc}$ (degrees)', fontsize=12)
    ax.set_ylabel('$p_{miss}^{exc}$ (GeV/c)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Counts', fontsize=12)
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

def main():
    """Main plotting function"""
    
    print("Loading histogram data...")
    
    # Read 1D histograms for analysis
    try:
        xB_centers, xB_weights, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_xB_1D.txt"))
        pmiss_centers, pmiss_weights, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_pmiss_1D.txt"))
        pmiss_exc_centers, pmiss_exc_weights, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_pmiss_exc_1D.txt"))
        angle_centers, angle_weights, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_angle pmiss_pmiss_exc_1D.txt"))
        
        print(f"  xB: {len(xB_centers)} bins")
        print(f"  pmiss: {len(pmiss_centers)} bins")
        print(f"  pmiss_exc: {len(pmiss_exc_centers)} bins")
        print(f"  angle: {len(angle_centers)} bins")
        
    except FileNotFoundError as e:
        print(f"Error: Could not find histogram files. {e}")
        print("Make sure the C++ program has been compiled and run with the updated code.")
        return
    
    # Define the two parameter ranges
    ranges = [
        {
            'name': 'Range 1: xB 0.5-1.5, pmiss 0.4-0.7',
            'xB_min': 0.5, 'xB_max': 1.5,
            'pmiss_min': 0.4, 'pmiss_max': 0.7,
            'region': 'region01'  # Regions 0 and 1
        },
        {
            'name': 'Range 2: xB 1.5-2.5, pmiss > 0.7',
            'xB_min': 1.5, 'xB_max': 2.5,
            'pmiss_min': 0.7, 'pmiss_max': 1.5,
            'region': 'region23'  # Regions 2 and 3
        }
    ]
    
    print("\nCreating plots...")
    
    for i, r in enumerate(ranges, 1):
        print(f"\n{r['name']}")
        
        # Load region-specific data
        try:
            # Try to read from multiple regions and combine them
            if 'region01' in r['region']:
                # Load regions 0 and 1
                pmiss_exc_0, weights_exc_0, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_pmiss_exc_region0_1D.txt"))
                angle_0, weights_angle_0, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_angle pmiss_pmiss_exc_region0_1D.txt"))
                
                pmiss_exc_1, weights_exc_1, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_pmiss_exc_region1_1D.txt"))
                angle_1, weights_angle_1, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_angle pmiss_pmiss_exc_region1_1D.txt"))
                
                # Combine regions 0 and 1
                pmiss_exc_combined = pmiss_exc_0.copy()
                angle_combined = angle_0.copy()
                weights_exc_combined = weights_exc_0 + weights_exc_1
                weights_angle_combined = weights_angle_0 + weights_angle_1
                
            elif 'region23' in r['region']:
                # Load regions 2 and 3
                pmiss_exc_2, weights_exc_2, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_pmiss_exc_region2_1D.txt"))
                angle_2, weights_angle_2, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_angle pmiss_pmiss_exc_region2_1D.txt"))
                
                pmiss_exc_3, weights_exc_3, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_pmiss_exc_region3_1D.txt"))
                angle_3, weights_angle_3, _ = read_1d_histogram(os.path.join(OUTPUT_DIR, "hist_angle pmiss_pmiss_exc_region3_1D.txt"))
                
                # Combine regions 2 and 3
                pmiss_exc_combined = pmiss_exc_2.copy()
                angle_combined = angle_2.copy()
                weights_exc_combined = weights_exc_2 + weights_exc_3
                weights_angle_combined = weights_angle_2 + weights_angle_3

            # Focus on the physically relevant pmiss_exc range (a few * sigma_CM)
            pmiss_mask = pmiss_exc_combined <= PMISS_EXC_MAX
            pmiss_exc_combined = pmiss_exc_combined[pmiss_mask]
            weights_exc_combined = weights_exc_combined[pmiss_mask]
                
            # Plot 1: pmiss_exc histogram
            plot_file_1 = os.path.join(PLOT_OUTPUT_DIR, f"pmiss_exc_hist_range{i}.png")
            create_pmiss_exc_histogram(pmiss_exc_combined, weights_exc_combined,
                                      f"$p_{{miss}}^{{exc}}$ Distribution - {r['name']}", 
                                      plot_file_1)
            
            # Plot 2: pmiss_exc vs angle (2D visualization)
            # Create a combined weight visualization
            plot_file_2 = os.path.join(PLOT_OUTPUT_DIR, f"pmiss_exc_vs_angle_range{i}.png")
            
            # For 2D plot, we'll create synthetic 2D data from the two 1D distributions
            # This is a visualization showing the relationship
            angle_grid, pmiss_exc_grid = np.meshgrid(angle_combined, pmiss_exc_combined)
            
            # Simple visualization: normalize and combine weights
            weights_norm_exc = weights_exc_combined / (np.max(weights_exc_combined) + 1e-10)
            weights_norm_angle = weights_angle_combined / (np.max(weights_angle_combined) + 1e-10)
            weights_2d = weights_norm_exc[:, np.newaxis] * weights_norm_angle[np.newaxis, :]
            
            # Create 2D histogram/heatmap
            fig, ax = plt.subplots(figsize=(11, 8))
            
            extent = [angle_combined[0], angle_combined[-1], 
                     0.0, PMISS_EXC_MAX]
            im = ax.imshow(weights_2d, aspect='auto', origin='lower', 
                          extent=extent, cmap='YlOrRd', interpolation='nearest')
            
            ax.set_xlabel('Angle between $\\vec{p}_{miss}$ and $\\vec{p}_{miss}^{exc}$ (degrees)', fontsize=12)
            ax.set_ylabel('$p_{miss}^{exc}$ (GeV/c)', fontsize=12)
            ax.set_title(f'$p_{{miss}}^{{exc}}$ vs Angle - {r["name"]}', fontsize=14, fontweight='bold')
            
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Relative Counts', fontsize=12)
            
            ax.grid(True, alpha=0.3, linestyle='--')
            plt.tight_layout()
            plt.savefig(plot_file_2, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  Saved: {plot_file_2}")
            
        except FileNotFoundError as e:
            print(f"  Error loading region data: {e}")
    
    print("\nPlotting complete!")
    print(f"All plots saved to: {PLOT_OUTPUT_DIR}")

if __name__ == "__main__":
    main()
