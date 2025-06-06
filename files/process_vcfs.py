import os  
import pandas as pd  d
import pysam  
import pyranges as pr 
import argparse  

def load_repeat_masker(file_path):
    columns = ["Chromosome", "Start", "End", "repClass", "repFamily"] 
    repeat_masker = pd.read_csv(file_path, sep="\t", comment="#", header=None,
                                usecols=[5, 6, 7, 11, 12], names=columns) 
    repeat_masker["Start"] = repeat_masker["Start"].astype(int)  
    repeat_masker["End"] = repeat_masker["End"].astype(int) 
    return pr.PyRanges(repeat_masker)  

def load_segmental_duplications(file_path):
    columns = ["Chromosome", "Start", "End", "name", "score"]  
    duplication_data = pd.read_csv(file_path, sep="\t", comment="#", header=None,
                                   usecols=[1, 2, 3, 4, 5], names=columns)  
    duplication_data["Start"] = duplication_data["Start"].astype(int)  
    duplication_data["End"] = duplication_data["End"].astype(int) 
    return pr.PyRanges(duplication_data) 

def extract_representative_sample(sample_names):
    return sample_names[0]  

def process_vcf(vcf_path, repeats, duplications, output_file, output_bed=False):
    with pysam.VariantFile(vcf_path) as vcf, open(output_file, "w") as writer:  
        sample_names = list(vcf.header.samples) 
        representative_sample = extract_representative_sample(sample_names) 
        
        if not output_bed: 
            columns = [
                "SAMPLE_NAME", "CHROM_CALLER", "POS_CALLER", "END_CALLER",
                "SVTYPE_CALLER", "SVLEN_CALLER", "MANTA", "DELLY", "SMOOVE",
                "OVERLAPS_REPEATS", "OVERLAPS_SEG_DUP"
            ]
            writer.write("\t".join(columns) + "\n")  

        vcf_data = [] 
        for record in vcf: 
            chrom = record.chrom 
            start = record.pos 
            end = record.stop or record.info.get("END") 
            svtype = record.info.get("SVTYPE", "NA")  
            svlen = record.info.get("SVLEN", "NA") 
            supp_vec = record.info.get("SUPP_VEC", "000")  

            if len(supp_vec) < 3:
                manta = delly = smoove = '0'  
            else:
                manta, delly, smoove = supp_vec[1], supp_vec[2], supp_vec[0]

            vcf_data.append({
                "Chromosome": chrom, "Start": start, "End": end,
                "SVTYPE_CALLER": svtype, "SVLEN_CALLER": svlen,
                "MANTA": manta, "DELLY": delly, "SMOOVE": smoove
            })

        vcf_ranges = pr.PyRanges(pd.DataFrame(vcf_data))

        overlap_repeats = vcf_ranges.join(repeats).as_df().drop_duplicates(subset=["Chromosome", "Start", "End"])
        overlap_duplications = vcf_ranges.join(duplications).as_df().drop_duplicates(subset=["Chromosome", "Start", "End"])

        repeat_set = set(zip(overlap_repeats["Chromosome"], overlap_repeats["Start"], overlap_repeats["End"]))
        duplication_set = set(zip(overlap_duplications["Chromosome"], overlap_duplications["Start"], overlap_duplications["End"]))

        for record in vcf_ranges.as_df().to_dict("records"): 
            overlaps_repeats = (record["Chromosome"], record["Start"], record["End"]) in repeat_set  
            overlaps_duplications = (record["Chromosome"], record["Start"], record["End"]) in duplication_set  
            
            start = record["Start"] - 1 if output_bed else record["Start"]  
            end = record["End"]

            row = [
                representative_sample, record["Chromosome"], str(start), str(end),
                record["SVTYPE_CALLER"], str(record["SVLEN_CALLER"]),
                record["MANTA"], record["DELLY"], record["SMOOVE"]
            ]

            row.append("True" if overlaps_repeats else ".")
            row.append("True" if overlaps_duplications else ".")

            if output_bed:
                writer.write("\t".join([record["Chromosome"], str(start), str(end), record["SVTYPE_CALLER"]]) + "\n")
            else:
                writer.write("\t".join(row) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Efficiently filter structural variants using PyRanges.") 
    parser.add_argument("--repeat_masker", required=True, help="Path to RepeatMasker TSV file")
    parser.add_argument("--duplications", required=True, help="Path to Segmental Duplications TSV file")
    parser.add_argument("--vcf_file", help="Path to a single VCF file")
    parser.add_argument("--output_file", help="Output file path")
    parser.add_argument("--vcf_folder", help="Path to folder containing input VCF files")
    parser.add_argument("--output_folder", help="Path to output folder")
    args = parser.parse_args()  

    repeats = load_repeat_masker(args.repeat_masker)
    duplications = load_segmental_duplications(args.duplications)

    if args.vcf_file:
        process_vcf(args.vcf_file, repeats, duplications, args.output_file)

    elif args.vcf_folder and args.output_folder:
        for vcf_file in os.listdir(args.vcf_folder):  
            if vcf_file.endswith(".vcf") or vcf_file.endswith(".vcf.gz"):  
                vcf_path = os.path.join(args.vcf_folder, vcf_file)
                output_path = os.path.join(
                    args.output_folder,
                    f"{vcf_file.replace('.vcf', '_processed.tsv').replace('.vcf.gz', '_processed.tsv')}"
                )
                process_vcf(vcf_path, repeats, duplications, output_path)
    else:
        raise ValueError("You must provide either --vcf_file and --output_file, or --vcf_folder and --output_folder")

if __name__ == "__main__":
    main()
