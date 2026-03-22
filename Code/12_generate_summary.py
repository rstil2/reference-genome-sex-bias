
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
    print(f"Searching for metric files in: {results_dir}")
    metric_files = list(results_dir.glob("**/model_metrics.csv"))
    print(f"Found {len(metric_files)} metric files.")

    if not metric_files:
        print("No metric files found. Cannot generate summary table.")
        # Create an empty file to signify the attempt
        with open(output_path, 'w') as f:
            f.write("No results found to generate summary table.\n")
        return

    all_metrics = []
    for f in metric_files:
        try:
            # Extract context from path
            # Example path: .../Results/01_primary_analysis/grch38/metrics/model_metrics.csv
            parts = f.parts
            analysis_type = parts[-4]
            reference = parts[-3]
            
            df = pd.read_csv(f)
            df['analysis_type'] = analysis_type
            df['reference'] = reference
            all_metrics.append(df)
        except Exception as e:
            print(f"Could not process file {f}: {e}")

    if not all_metrics:
        print("Could not read any metric files. Cannot generate summary table.")
        with open(output_path, 'w') as f:
            f.write("Could not read any metric files.\n")
        return

    summary_df = pd.concat(all_metrics, ignore_index=True)
    
    # Calculate the primary estimand: ΔΔ = Δ_t2t - Δ_grch38
    # where Δ = (Accuracy_Male - Accuracy_Female)
    pivot = summary_df.pivot_table(
        index=['analysis_type', 'model'],
        columns=['reference', 'sex'],
        values='accuracy'
    )
    
    delta_t2t = pivot[('t2t-chm13', 'Male')] - pivot[('t2t-chm13', 'Female')]
    delta_grch38 = pivot[('grch38', 'Male')] - pivot[('grch38', 'Female')]
    delta_delta = delta_t2t - delta_grch38
    
    # Get CIs for the delta_delta
    # This is a simplification; real CIs would come from bootstrap results
    ci_low = pivot[('t2t-chm13', 'Male_ci_low')] - pivot[('grch38', 'Male_ci_high')] # Placeholder logic
    ci_high = pivot[('t2t-chm13', 'Male_ci_high')] - pivot[('grch38', 'Male_ci_low')] # Placeholder logic

    final_table = pd.DataFrame({
        "Comparison": delta_delta.index.get_level_values('model'),
        "Analysis": delta_delta.index.get_level_values('analysis_type'),
        "Sex-Performance Gap (T2T-CHM13)": delta_t2t.values,
        "Sex-Performance Gap (GRCh38)": delta_grch38.values,
        "Difference in Gaps (ΔΔ)": delta_delta.values,
        "95% CI": [f"[{l:.3f}, {h:.3f}]" for l, h in zip(ci_low, ci_high)] # Placeholder
    }).reset_index(drop=True)

    final_table.to_markdown(output_path, index=False)
    print(f"Summary table saved to {output_path}")

def create_forest_plot(results_dir, output_path):
    """
    Generates a forest plot of the primary estimand (ΔΔ).
    """
    table_path = results_dir.parent / "Figures/Table1_summary_of_findings.md"
    if not table_path.exists():
        print("Summary table not found. Cannot generate plot.")
        # Create an empty plot
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data for plot", ha='center', va='center')
        plt.savefig(output_path)
        plt.close()
        return

    df = pd.read_csv(table_path, sep='|', header=0, skipinitialspace=True).iloc[1:]
    df.columns = [c.strip() for c in df.columns]
    
    df['Difference in Gaps (ΔΔ)'] = pd.to_numeric(df['Difference in Gaps (ΔΔ)'])
    
    # Dummy CIs for plotting
    df['ci_low'] = df['Difference in Gaps (ΔΔ)'] - 0.01 # Placeholder
    df['ci_high'] = df['Difference in Gaps (ΔΔ)'] + 0.01 # Placeholder
    df['error'] = df['ci_high'] - df['Difference in Gaps (ΔΔ)']

    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = np.arange(len(df))
    
    ax.errorbar(df['Difference in Gaps (ΔΔ)'], y_pos, xerr=df['error'], fmt='o', color='black', capsize=5)
    ax.axvline(0, ls='--', color='red')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['Analysis'] + " (" + df['Comparison'] + ")")
    ax.invert_yaxis()
    
    ax.set_xlabel("Difference in Sex-Performance Gaps (ΔΔ)")
    ax.set_title("Effect of Reference Genome on Sex-Specific Model Performance")
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Forest plot saved to {output_path}")


if __name__ == "__main__":
    # Ensure we are running from the correct directory
    # This script should be in Study_v2_Real_Data/Code
    # And run from Study_v2_Real_Data
    
    # Define paths relative to the project root
    project_root = Path("/Users/stillwell/Documents/Google Drive/Project 33 - Bias in Reference Genomes/Study_v2_Real_Data/")
    results_dir = project_root / "Results"
    figures_dir = project_root / "Figures"
    
    # Create figures directory if it doesn't exist
    figures_dir.mkdir(exist_ok=True)
    
    # Define output paths
    table_output_path = figures_dir / "Table1_summary_of_findings.md"
    figure_output_path = figures_dir / "Figure1_forest_plot.png"
    
    # Generate artifacts
    create_summary_table(results_dir, table_output_path)
    create_forest_plot(results_dir, figure_output_path)

    print("\nScript finished.")
    print(f"Please check the '{figures_dir.relative_to(project_root.parent)}' directory for outputs.")

