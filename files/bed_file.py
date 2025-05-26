import os
import sys
import pandas as pd
import pyranges as pr
 
def load_bed_file(bed_path, value_col_name):
    df = pd.read_csv(bed_path, sep="\t", header=None,
                     names=["Chromosome", "Start", "End", "SVTYPE", value_col_name])
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df_reduced = df[["Chromosome", "Start", "End", value_col_name]]
    return pr.PyRanges(df_reduced)
 
def add_gc_and_coverage(tsv_path, gc_bed_path, coverage_bed_path, output_path):
    print(f"Processing {tsv_path} ...")
    # Load original TSV, keep all columns
    df = pd.read_csv(tsv_path, sep="\t")
 
    # Prepare PyRanges intervals (convert to 0-based start)
    df_bed = df.rename(columns={"CHROM_CALLER": "Chromosome",
                               "POS_CALLER": "Start",
                               "END_CALLER": "End"})
    df_bed["Start"] = df_bed["Start"] - 1  # 0-based start for BED
 
    intervals = pr.PyRanges(df_bed[["Chromosome", "Start", "End"]])
 
    gc_ranges = load_bed_file(gc_bed_path, "GC_CONTENT%")
    cov_ranges = load_bed_file(coverage_bed_path, "COVERAGE_MOSDEPTH")
 
    joined_gc = intervals.join(gc_ranges)
    joined_cov = intervals.join(cov_ranges)
 
    gc_df = joined_gc.df[["Chromosome", "Start", "End", "GC_CONTENT%"]]
    cov_df = joined_cov.df[["Chromosome", "Start", "End", "COVERAGE_MOSDEPTH"]]
 
    gc_df["key"] = gc_df["Chromosome"].astype(str) + ":" + gc_df["Start"].astype(str) + "-" + gc_df["End"].astype(str)
    cov_df["key"] = cov_df["Chromosome"].astype(str) + ":" + cov_df["Start"].astype(str) + "-" + cov_df["End"].astype(str)
    df_bed["key"] = df_bed["Chromosome"].astype(str) + ":" + df_bed["Start"].astype(str) + "-" + df_bed["End"].astype(str)
 
    gc_agg = gc_df.groupby("key")["GC_CONTENT%"].mean().reset_index()
    cov_agg = cov_df.groupby("key")["COVERAGE_MOSDEPTH"].mean().reset_index()
 
    df_bed = df_bed.merge(gc_agg, on="key", how="left")
    df_bed = df_bed.merge(cov_agg, on="key", how="left")
 
    df_bed["GC_CONTENT%"] = df_bed["GC_CONTENT%"].fillna("NA")
    df_bed["COVERAGE_MOSDEPTH"] = df_bed["COVERAGE_MOSDEPTH"].fillna("NA")
 
    df["GC_CONTENT%"] = df_bed["GC_CONTENT%"].values
    df["COVERAGE_MOSDEPTH"] = df_bed["COVERAGE_MOSDEPTH"].values
 
    df.to_csv(output_path, sep="\t", index=False)
    print(f"Saved TSV with added GC_CONTENT% and COVERAGE_MOSDEPTH columns to: {output_path}")
 
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python update_to_bed.py <input_tsv> <gc_bed> <coverage_bed> <output_tsv>")
        sys.exit(1)
 
    input_tsv = sys.argv[1]
    gc_bed = sys.argv[2]
    coverage_bed = sys.argv[3]
    output_tsv = sys.argv[4]
 
    add_gc_and_coverage(input_tsv, gc_bed, coverage_bed, output_tsv)
