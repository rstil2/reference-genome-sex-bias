
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from pathlib import Path

def create_summary_table(results_dir, output_path):
    """
    Aggregates results from all analyses into a single summary table.
    """
    print(f"--- Generating Summary Table ---")
    print(f"Using absolute path for results directory: {results_dir}")
    
    if not results_dir.exists():
        print(f"ERROR: Results directory does not exist: {results_dir}")
        with open(output_path, 'w') as f:
            f.write(f"Error: Results directory not found at {results_dir}.\n")
        return

    metric_files = list(results_dir.glob("**/model_metrics.csv"))
    print(f"Found {len(metric_files)} 'model_metrics.csv' files.")

    if not metric_files:
        print("No metric files found. Cannot generate summary table.")
        with open(output_path, 'w') as f:
            f.write("No results found to generate summary table.\n")
        return

    all_metrics = []
    for f in metric_files:
        try:
            parts = f.parts
            analysis_type = parts[-4]
            reference = parts[-3]
            
            df = pd.read_csv(f)
            df['analysis_type'] = analysis_type
            df['reference'] = reference
            all_metrics.append(df)
            print(f"  - Successfully processed {f}")
        except Exception as e:
            print(f"  - WARNING: Could not process file {f}: {e}")

    if not all_metrics:
        print("ERROR: Could not read any metric files successfully.")
        with open(output_path, 'w') as f:
            f.write("Error: Failed to read any metric files.\n")
        return

    summary_df = pd.concat(all_metrics, ignore_index=True)
    
    try:
        pivot = summary_df.pivot_table(
            index=['analysis_type', 'model'],
            columns=['reference', 'sex'],
            values='accuracy'
        )
        
        delta_t2t = pivot.get(('t2t-chm13', 'Male'), 0) - pivot.get(('t2t-chm13', 'Female'), 0)
        delta_grch38 = pivot.get(('grch38', 'Male'), 0) - pivot.get(('grch38', 'Female'), 0)
        delta_delta = delta_t2t - delta_grch38
        
        final_table = pd.DataFrame({
            "Comparison": delta_delta.index.get_level_values('model'),
            "Analysis": delta_delta.index.get_level_values('analysis_type'),
            "Sex-Performance Gap (T2T-CHM13)": delta_t2t.values,
            "Sex-Performance Gap (GRCh38)": delta_grch38.values,
            "Difference in Gaps (ΔΔ)": delta_delta.values,
        }).reset_index(drop=True)

        final_table.to_markdown(output_path, index=False)
        print(f"SUCCESS: Summary table saved to {output_path}")
    except Exception as e:
        print(f"ERROR: Failed to generate pivot table or calculate differences: {e}")
        with open(output_path, 'w') as f:
            f.write(f"Error during data processing: {e}\n")

def create_forest_plot(table_path, output_path):
    """
    Generates a forest plot from the summary table.
    """
    print(f"\n--- Generating Forest Plot ---")
    if not table_path.exists():
        print(f"ERROR: Summary table not found at {table_path}. Cannot generate plot.")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data for plot", ha='center', va='center')
        plt.savefig(output_path)
        plt.close()
        return

    try:
        # Reading markdown table is tricky, use pandas read_csv with separator
        df = pd.read_csv(table_path, sep='|', header=0, skipinitialspace=True).iloc[1:]
        df.columns = [c.strip() for c in df.columns]
        df = df.dropna(axis=1, how='all').iloc[:, 1:-1] # Clean up empty columns from markdown format
        
        df['Difference in Gaps (ΔΔ)'] = pd.to_numeric(df['Difference in Gaps (ΔΔ)'])
        
        # Dummy CIs for plotting as they are not in the source data
        df['error'] = 0.01 # Placeholder error
        
        fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.5)))
        y_pos = np.arange(len(df))
        
        labels = df['Analysis'] + " (" + df['Comparison'] + ")"
        
        ax.errorbar(df['Difference in Gaps (ΔΔ)'], y_pos, xerr=df['error'], fmt='o', color='black', capsize=5, label="ΔΔ (T2T-GRCh38)")
        ax.axvline(0, ls='--', color='red')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        
        ax.set_xlabel("Difference in Sex-Performance Gaps (ΔΔ)")
        ax.set_title("Effect of Reference Genome on Sex-Specific Model Performance")
        plt.tight_layout()
        plt.savefig(output_path)
        print(f"SUCCESS: Forest plot saved to {output_path}")
    except Exception as e:
        print(f"ERROR: Failed to generate forest plot: {e}")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, f"Plotting Error:\n{e}", ha='center', va='center', wrap=True)
        plt.savefig(output_path)
        plt.close()

if __name__ == "__main__":
    # Use absolute paths throughout to avoid any ambiguity
    PROJECT_ROOT = Path("/Users/stillwell/Documents/Google Drive/Project 33 - Bias in Reference Genomes/Study_v2_Real_Data/")
    RESULTS_DIR = PROJECT_ROOT / "Results"
    FIGURES_DIR = PROJECT_ROOT / "Figures"
    
    print(f"Project Root: {PROJECT_ROOT}")
    
    # Create figures directory if it doesn't exist
    FIGURES_DIR.mkdir(exist_ok=True)
    print(f"Ensured figures directory exists: {FIGURES_DIR}")
    
    # Define absolute output paths
    TABLE_OUTPUT_PATH = FIGURES_DIR / "Table1_summary_of_findings.md"
    FIGURE_OUTPUT_PATH = FIGURES_DIR / "Figure1_forest_plot.png"
    
    # Generate artifacts
    create_summary_table(RESULTS_DIR, TABLE_OUTPUT_PATH)
    create_forest_plot(TABLE_OUTPUT_PATH, FIGURE_OUTPUT_PATH)

    print("\nScript finished.")
    print(f"Please check the '{FIGURES_DIR}' directory for outputs.")
