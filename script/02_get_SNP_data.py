#!/usr/bin/env python3
import os
import sys
import subprocess
import csv

# === Configuration (Adapt if necessary!) ===

# Input TSV file with regions (Requires header: Chromosome, Position_GRCh37)
# Replace with the actual path to your TSV file
REGIONS_FILE_TSV = "data/snp_coordinates_grch37.tsv" # <--- !! UPDATE THIS !!

# Local panel file (provided by user - CHECK NAME!)
PANEL_FILE_LOCAL = "data/integrated_call_samples_v3.20130502.ALL.panel.txt"

# File to save EUR sample IDs
SAMPLES_FILE = "ID_files/sas_samples.txt"

# Target super-population code
POPULATION_CODE = "SAS"

# --- VCF URL is now constructed dynamically below ---
# Base path for the 1000 Genomes VCF files (GRCh37 release 20130502)
# (!!! VERIFY THIS BASE PATH AND FILENAME STRUCTURE IF USING A DIFFERENT RELEASE/SOURCE !!!)
BASE_VCF_FTP_PATH = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
VCF_FILENAME_TEMPLATE = "ALL.chr{CHROMOSOME}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

# Command for bcftools (Use the full path if needed)
BCFTOOLS_CMD = "/opt/anaconda3/envs/bcftools/bin/bcftools" # <--- PUT THE PATH FOUND WITH 'which bcftools' HERE

# Size of the window around the position in the TSV file (e.g., 5000 means Position +/- 5000 bp)
REGION_WINDOW_SIZE = 10000

# Directory to store output VCF files
OUTPUT_DIR = "output_vcfs"

# === Helper Functions ===

def extract_samples(panel_file, population_code, samples_output_file):
    """
    Extracts sample IDs for a specific population from the panel file.
    Returns the count of samples found.
    """
    print(f"\n--- Step 2: Extracting {population_code} sample IDs ---")
    sample_count = 0
    try:
        with open(panel_file, 'r') as infile, open(samples_output_file, 'w') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            header = next(reader) # Skip header
            # Find column indices (more robust than assuming fixed positions)
            try:
                sample_col_idx = header.index('sample') # or the actual header name for sample ID
                pop_col_idx = header.index('super_pop') # or the actual header name for super population
            except ValueError as e:
                print(f"ERROR: Missing expected column in panel file header: {e}")
                print(f"Expected columns like 'sample' and 'super_pop'. Found: {header}")
                sys.exit(1)

            for row in reader:
                 # Check row length before accessing indices
                 if len(row) > max(sample_col_idx, pop_col_idx) and row[pop_col_idx] == population_code:
                    outfile.write(row[sample_col_idx] + '\n')
                    sample_count += 1

        if sample_count > 0:
            print(f"Extraction successful. {sample_count} {population_code} IDs saved to {samples_output_file}")
            return sample_count
        else:
            print(f"ERROR: No samples found for {population_code} in '{panel_file}' (check column {pop_col_idx+1} and population code).")
            sys.exit(1)

    except FileNotFoundError:
        print(f"ERROR: Panel file '{panel_file}' not found!")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to process panel file '{panel_file}': {e}")
        sys.exit(1)


def run_bcftools(bcftools_cmd, vcf_url_for_chr, region, samples_file, output_vcf_file):
    """
    Runs the bcftools view command for a specific region and sample set,
    using the VCF URL appropriate for the chromosome in the region.
    """
    print(f"\n--- Step 3: Extracting region {region} with bcftools ---")
    print(f"  Region target (GRCh37): {region}")
    # Note: vcf_url_for_chr is now specific to the chromosome being processed
    print(f"  VCF source (GRCh37): {vcf_url_for_chr}")
    print(f"  Target samples: {samples_file}")
    print(f"  Output file: {output_vcf_file}")

    command = [
        bcftools_cmd,
        "view",
        vcf_url_for_chr, # Use the dynamically generated URL passed to this function
        "-r", region,
        "-S", samples_file,
        "-O", "z",         # Output compressed VCF
        "-o", output_vcf_file
    ]

    try:
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_vcf_file), exist_ok=True)

        # Run bcftools command
        result = subprocess.run(command, check=False, capture_output=True, text=True) # Use check=False to handle errors manually

        # Check if bcftools failed
        if result.returncode != 0:
            print(f"ERROR: bcftools failed for region {region}.")
            print(f"  Return Code: {result.returncode}")
            print(f"  Stderr: {result.stderr.strip()}") # Strip whitespace from stderr
            print(f"  Stdout: {result.stdout.strip()}") # Strip whitespace from stdout
            print("  Check VCF URL (is it valid for this chromosome? is the server accessible?), region format,")
            print("  existence of the index (.tbi) on the server for this chromosome file,")
            print(f"  and the sample file ('{samples_file}').")
            # Indicate failure for this region
            return False
        else:
            print(f"bcftools extraction completed successfully for region {region}.")
            # bcftools might print non-error info to stderr (e.g., index loading messages)
            if result.stderr:
                 print(f"  bcftools stderr output:\n{result.stderr.strip()}")
            # Indicate success
            return True

    except FileNotFoundError:
        print(f"ERROR: bcftools command '{bcftools_cmd}' not found.")
        print("  Please ensure bcftools is installed and the path in BCFTOOLS_CMD is correct.")
        sys.exit(1) # Exit here as bcftools is fundamental
    except Exception as e:
        print(f"ERROR: An unexpected error occurred while running bcftools for region {region}: {e}")
        # Indicate failure for this region
        return False

# === Main Script Execution ===

def main():
    print("--- Step 1: Checking local panel file ---")
    if not os.path.isfile(PANEL_FILE_LOCAL):
        print(f"ERROR: Local panel file '{PANEL_FILE_LOCAL}' not found!")
        sys.exit(1)
    else:
        print(f"Using local panel file: {PANEL_FILE_LOCAL}")

    # Extract sample IDs *once* before the loop
    extract_samples(PANEL_FILE_LOCAL, POPULATION_CODE, SAMPLES_FILE)

    print(f"\n--- Step 4: Processing regions from TSV file: {REGIONS_FILE_TSV} ---")
    if not os.path.isfile(REGIONS_FILE_TSV):
        print(f"ERROR: Regions TSV file '{REGIONS_FILE_TSV}' not found!")
        sys.exit(1)

    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output VCFs will be saved in: {OUTPUT_DIR}")

    successful_extractions = []
    failed_regions = []

    try:
        with open(REGIONS_FILE_TSV, 'r', newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            header = next(reader) # Read header row

            # Find column indices dynamically
            try:
                chr_col_idx = header.index("Chromosome")
                pos_col_idx = header.index("Position_GRCh37")
                # Optional: Get rsID for naming output files
                rsid_col_idx = header.index("rsID") if "rsID" in header else -1
            except ValueError as e:
                print(f"ERROR: Missing expected column in regions TSV file header: {e}")
                print(f"Expected columns like 'Chromosome' and 'Position_GRCh37'. Found: {header}")
                sys.exit(1)

            # Process each row (region) from the TSV
            for i, row in enumerate(reader):
                try:
                    chromosome = row[chr_col_idx]
                    # Handle potential 'chr' prefix if needed, although 1000G usually uses just numbers/X/Y/MT
                    chromosome_cleaned = chromosome.replace('chr', '')
                    position = int(row[pos_col_idx])
                    rsid = row[rsid_col_idx] if rsid_col_idx != -1 else f"region_{i+1}"

                    # Define the region boundaries
                    start_pos = max(1, position - REGION_WINDOW_SIZE) # Ensure start is not negative
                    end_pos = position + REGION_WINDOW_SIZE
                    # Use the cleaned chromosome for region string if needed by bcftools, but original for filename likely
                    region_str = f"{chromosome_cleaned}:{start_pos}-{end_pos}" # bcftools usually expects e.g. '1' not 'chr1'

                    # --- Dynamically create the VCF URL for the current chromosome ---
                    vcf_filename_on_server = VCF_FILENAME_TEMPLATE.replace("{CHROMOSOME}", chromosome_cleaned)
                    dynamic_vcf_url = BASE_VCF_FTP_PATH + vcf_filename_on_server
                    # --- End of dynamic URL creation ---

                    # Define the output file name (use original chromosome string from TSV for clarity)
                    output_filename = os.path.join(
                        OUTPUT_DIR,
                        f"{rsid}_chr{chromosome}_{start_pos}_{end_pos}_GRCh37_{POPULATION_CODE}_subset.vcf.gz"
                    )

                    # Run bcftools for this region using the dynamically constructed URL
                    if run_bcftools(BCFTOOLS_CMD, dynamic_vcf_url, region_str, SAMPLES_FILE, output_filename):
                         successful_extractions.append(output_filename)
                    else:
                         failed_regions.append(f"{region_str} (Source: {dynamic_vcf_url})") # Add more context on failure

                except (IndexError, ValueError) as e:
                    print(f"WARNING: Skipping invalid row {i+2} in {REGIONS_FILE_TSV}: {row}. Error: {e}")
                    failed_regions.append(f"Row {i+2} (Invalid Data: {row})")
                except Exception as e:
                    print(f"WARNING: An unexpected error occurred processing row {i+2} ({row}) from {REGIONS_FILE_TSV}: {e}")
                    failed_regions.append(f"Row {i+2} (Processing Error: {row})")


    except FileNotFoundError:
        print(f"ERROR: Regions TSV file '{REGIONS_FILE_TSV}' not found during processing!")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read or process regions TSV file '{REGIONS_FILE_TSV}': {e}")
        sys.exit(1)

    print("\n--- Optional Cleanup ---")
    print(f"Sample ID file ({SAMPLES_FILE}) kept.")
    # If you want to remove the sample file after completion, uncomment the next line:
    # try:
    #     os.remove(SAMPLES_FILE)
    #     print(f"Removed temporary sample file: {SAMPLES_FILE}")
    # except OSError as e:
    #     print(f"Warning: Could not remove sample file {SAMPLES_FILE}: {e}")


    print("\n=== SCRIPT FINISHED ===")
    print(f"Attempted to process {i+1} regions from TSV.") # i is the index from enumerate starting at 0
    print(f"Successfully created VCFs ({len(successful_extractions)}):")
    # Optional: print list of successful files if needed
    # for fname in successful_extractions:
    #     print(f"  - {fname}")

    if failed_regions:
        print(f"\nFailed or skipped regions/rows ({len(failed_regions)}):")
        for region_info in failed_regions:
            print(f"  - {region_info}")
    else:
        print("\nAll regions processed successfully (check bcftools output above for details on empty regions).")

if __name__ == "__main__":
    main()