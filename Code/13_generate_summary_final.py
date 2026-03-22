
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from pathlib import Path
import traceback

def create_summary_table(results_dir, output_path):
    """
    Aggregates results from all analyses into a single summary table.
    """
    print(f"--- Generating Summary Table ---")
    print(f"Searching for metric files in: {results_dir}")
    
    if not results_dir.exists():
        error_msg = f"FATAL: Results directory does not exist: {results_dir}"
        print(error_msg)
        with open(output_path, 'w') as f:
            f.write(error_msg + "\n")
        return False

    metric_files = list(results_dir.glob("**/model_metrics.csv"))
    print(f"Found {len(metric_files)} 'model_metrics.csv' files.")

    if not metric_files:
        error_msg = "FATAL: No 'model_metrics.csv' files found. Cannot generate summary."
        print(error_msg)
        with open(output_path, 'w') as f:
            f.write(error_msg + "\n")
        return False

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
            print(f"  - Successfully processed: {f.relative_to(results_dir.parent)}")
        except Exception as e:
            print(f"  - WARNING: Could not process file {f}. Error: {e}")

    if not all_metrics:
        error_msg = "FATAL: All metric files failed to process."
        print(error_msg)
        with open(output_path, 'w') as f:
            f.write(error_msg + "\n")
        return False

    summary_df = pd.concat(all_metrics, ignore_index=True)
    
    try:
        print("Pivoting data to calculate performance gaps...")
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
            "Sex-Perf-Gap (T2T)": delta_t2t.values,
            "Sex-Perf-Gap (GRCh38)": delta_grch38.values,
            "Difference in Gaps (ΔΔ)": delta_delta.values,
        }).reset_index(drop=True)

        final_table.to_markdown(output_path, index=False)
        print(f"SUCCESS: Summary table saved to {output_path}")
        return True
    except Exception as e:
        error_msg = f"FATAL: Failed to generate pivot table. Error: {e}\n{traceback.format_exc()}"
        print(error_msg)
        with open(output_path, 'w') as f:
            f.write(error_msg + "\n")
        return False

def create_forest_plot(table_path, output_path):
    """
    Generates a forest plot from the summary table.
    """
    print(f"\n--- Generating Forest Plot ---")
    if not table_path.exists() or os.path.getsize(table_path) < 50:
        error_msg = f"FATAL: Valid summary table not found at {table_path}. Cannot generate plot."
        print(error_msg)
        # Create an empty plot with error message
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data for plot.\nSummary table generation failed.", ha='center', va='center', color='red')
        plt.savefig(output_path)
        plt.close()
        return

    try:
        with open(table_path, 'r') as f:
            lines = f.readlines()
        
        # Clean up markdown formatting for robust parsing
        header = [h.strip() for h in lines[0].strip('|').split('|')]
        data = [list(map(str.strip, line.strip('|').split('|'))) for line in lines[2:]]
        
        df = pd.DataFrame(data, columns=header)
        df = df.apply(pd.to_numeric, errors='ignore')

        df['error'] = 0.01 # Placeholder error for plot
        
        fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.6)))
        y_pos = np.arange(len(df))
        
        labels = df['Analysis'] + " (" + df['Comparison'] + ")"
        
        ax.errorbar(df['Difference in Gaps (ΔΔ)'], y_pos, xerr=df['error'], fmt='o', color='black', capsize=5)
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
        error_msg = f"FATAL: Failed to generate forest plot. Error: {e}\n{traceback.format_exc()}"
        print(error_msg)
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, f"Plotting Error:\n{e}", ha='center', va='center', wrap=True, color='red')
        plt.savefig(output_path)
        plt.close()

if __name__ == "__main__":
    try:
        PROJECT_ROOT = Path("/Users/stillwell/Documents/Google Drive/Project 33 - Bias in Reference Genomes/Study_v2_Real_Data/")
        RESULTS_DIR = PROJECT_ROOT / "Results"
        FIGURES_DIR = PROJECT_ROOT / "Figures"
        
        print(f"Project Root: {PROJECT_ROOT}")
        
        FIGURES_DIR.mkdir(exist_ok=True)
        print(f"Ensured figures directory exists: {FIGURES_DIR}")
        
        TABLE_OUTPUT_PATH = FIGURES_DIR / "Table1_summary_of_findings.md"
        FIGURE_OUTPUT_PATH = FIGURES_DIR / "Figure1_forest_plot.png"
        
        if create_summary_table(RESULTS_DIR, TABLE_OUTPUT_PATH):
            create_forest_plot(TABLE_OUTPUT_PATH, FIGURE_OUTPUT_PATH)

        print("\n--- Script finished ---")
        print(f"Please check the '{FIGURES_DIR}' directory for outputs.")
    except Exception as e:
        print(f"An unexpected error occurred at the top level: {e}\n{traceback.format_exc()}")

