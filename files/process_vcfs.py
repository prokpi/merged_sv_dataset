import os  # For file and folder operations
import pandas as pd  # For data manipulation and reading TSV/CSV files
import pysam  # For reading and parsing VCF files
import pyranges as pr  # For efficient genomic interval operations
import argparse  # For handling command-line arguments

# Function to load RepeatMasker data and convert it into a PyRanges object
def load_repeat_masker(file_path):
    columns = ["Chromosome", "Start", "End", "repClass", "repFamily"]  # Column names for RepeatMasker file
    repeat_masker = pd.read_csv(file_path, sep="\t", comment="#", header=None,
                                usecols=[5, 6, 7, 11, 12], names=columns)  # Read specific columns
    repeat_masker["Start"] = repeat_masker["Start"].astype(int)  # Ensure "Start" is an integer
    repeat_masker["End"] = repeat_masker["End"].astype(int)  # Ensure "End" is an integer
    return pr.PyRanges(repeat_masker)  # Convert to PyRanges for fast genomic operations

# Function to load segmental duplications data and convert it into a PyRanges object
def load_segmental_duplications(file_path):
    columns = ["Chromosome", "Start", "End", "name", "score"]  # Column names for segmental duplications file
    duplication_data = pd.read_csv(file_path, sep="\t", comment="#", header=None,
                                   usecols=[1, 2, 3, 4, 5], names=columns)  # Read specific columns
    duplication_data["Start"] = duplication_data["Start"].astype(int)  # Ensure "Start" is an integer
    duplication_data["End"] = duplication_data["End"].astype(int)  # Ensure "End" is an integer
    return pr.PyRanges(duplication_data)  # Convert to PyRanges

# Function to select the representative sample from the VCF header
def extract_representative_sample(sample_names):
    return sample_names[0]  # Simply return the first sample in the list

# Function to process a VCF file and annotate structural variants with repeat and duplication overlaps
def process_vcf(vcf_path, repeats, duplications, output_file, output_bed=False):
    with pysam.VariantFile(vcf_path) as vcf, open(output_file, "w") as writer:  # Open VCF and output file
        sample_names = list(vcf.header.samples)  # Get the sample names from the VCF header
        representative_sample = extract_representative_sample(sample_names)  # Select a representative sample
        
        if not output_bed:  # If the output is not in BED format, write the TSV header
            columns = [
                "SAMPLE_NAME", "CHROM_CALLER", "POS_CALLER", "END_CALLER",
                "SVTYPE_CALLER", "SVLEN_CALLER", "MANTA", "DELLY", "SMOOVE",
                "OVERLAPS_REPEATS", "OVERLAPS_SEG_DUP"
            ]
            writer.write("\t".join(columns) + "\n")  # Write the column names to the output file

        vcf_data = []  # List to store variant data
        for record in vcf:  # Iterate through each record in the VCF
            chrom = record.chrom  # Chromosome name
            start = record.pos  # Start position
            end = record.stop or record.info.get("END")  # End position
            svtype = record.info.get("SVTYPE", "NA")  # Structural variant type
            svlen = record.info.get("SVLEN", "NA")  # Structural variant length
            supp_vec = record.info.get("SUPP_VEC", "000")  # Support vector from callers

            # Parse support vector into individual callers
            if len(supp_vec) < 3:
                manta = delly = smoove = '0'  # Default to '0' if vector length is insufficient
            else:
                manta, delly, smoove = supp_vec[1], supp_vec[2], supp_vec[0]

            # Add variant data to the list
            vcf_data.append({
                "Chromosome": chrom, "Start": start, "End": end,
                "SVTYPE_CALLER": svtype, "SVLEN_CALLER": svlen,
                "MANTA": manta, "DELLY": delly, "SMOOVE": smoove
            })

        # Convert the list of variant data to a PyRanges object
        vcf_ranges = pr.PyRanges(pd.DataFrame(vcf_data))

        # Check overlaps with repeats and duplications
        overlap_repeats = vcf_ranges.join(repeats).as_df().drop_duplicates(subset=["Chromosome", "Start", "End"])
        overlap_duplications = vcf_ranges.join(duplications).as_df().drop_duplicates(subset=["Chromosome", "Start", "End"])

        # Convert overlaps into sets for quick lookup
        repeat_set = set(zip(overlap_repeats["Chromosome"], overlap_repeats["Start"], overlap_repeats["End"]))
        duplication_set = set(zip(overlap_duplications["Chromosome"], overlap_duplications["Start"], overlap_duplications["End"]))

        # Write annotated variants to the output file
        for record in vcf_ranges.as_df().to_dict("records"):  # Iterate through each variant
            overlaps_repeats = (record["Chromosome"], record["Start"], record["End"]) in repeat_set  # Check repeat overlap
            overlaps_duplications = (record["Chromosome"], record["Start"], record["End"]) in duplication_set  # Check duplication overlap
            
            start = record["Start"] - 1 if output_bed else record["Start"]  # Convert to 0-based if BED format
            end = record["End"]

            # Create a row for TSV output
            row = [
                representative_sample, record["Chromosome"], str(start), str(end),
                record["SVTYPE_CALLER"], str(record["SVLEN_CALLER"]),
                record["MANTA"], record["DELLY"], record["SMOOVE"]
            ]

            # Add overlap flags
            row.append("True" if overlaps_repeats else ".")
            row.append("True" if overlaps_duplications else ".")

            # Write row to the file
            if output_bed:
                writer.write("\t".join([record["Chromosome"], str(start), str(end), record["SVTYPE_CALLER"]]) + "\n")
            else:
                writer.write("\t".join(row) + "\n")

# Main function to parse arguments and orchestrate file processing
def main():
    parser = argparse.ArgumentParser(description="Efficiently filter structural variants using PyRanges.")  # Set up argument parser
    parser.add_argument("--repeat_masker", required=True, help="Path to RepeatMasker TSV file")
    parser.add_argument("--duplications", required=True, help="Path to Segmental Duplications TSV file")
    parser.add_argument("--vcf_file", help="Path to a single VCF file")
    parser.add_argument("--output_file", help="Output file path")
    parser.add_argument("--vcf_folder", help="Path to folder containing input VCF files")
    parser.add_argument("--output_folder", help="Path to output folder")
    args = parser.parse_args()  # Parse arguments

    # Load RepeatMasker and duplication data
    repeats = load_repeat_masker(args.repeat_masker)
    duplications = load_segmental_duplications(args.duplications)

    # Process a single VCF file
    if args.vcf_file:
        process_vcf(args.vcf_file, repeats, duplications, args.output_file)

    # Process multiple VCF files in a folder
    elif args.vcf_folder and args.output_folder:
        for vcf_file in os.listdir(args.vcf_folder):  # Iterate through files in the folder
            if vcf_file.endswith(".vcf") or vcf_file.endswith(".vcf.gz"):  # Check file extension
                vcf_path = os.path.join(args.vcf_folder, vcf_file)
                output_path = os.path.join(
                    args.output_folder,
                    f"{vcf_file.replace('.vcf', '_processed.tsv').replace('.vcf.gz', '_processed.tsv')}"
                )
                process_vcf(vcf_path, repeats, duplications, output_path)
    else:
        raise ValueError("You must provide either --vcf_file and --output_file, or --vcf_folder and --output_folder")

# Entry point for the script
if __name__ == "__main__":
    main()
