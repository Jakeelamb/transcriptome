#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import glob
from datetime import datetime

def setup_directories(root_dir=None):
    """Set up necessary directories for transcriptome analysis."""
    if root_dir is None:
        root_dir = os.getcwd()
    
    dirs = {
        'project': root_dir,
        'raw_reads': os.path.join(root_dir, 'data', 'raw_reads'),
        'trimmed_reads': os.path.join(root_dir, 'results', 'trimmed_reads'),
        'transcriptome': os.path.join(root_dir, 'data', 'transcriptome'),
        'results': os.path.join(root_dir, 'results'),
        'trimming_results': os.path.join(root_dir, 'results', 'trimming'),
        'merging_results': os.path.join(root_dir, 'results', 'merging'),
        'normalization_results': os.path.join(root_dir, 'results', 'normalization'),
        'assembly_results': os.path.join(root_dir, 'results', 'assembly'),
        'busco_results': os.path.join(root_dir, 'results', 'busco'),
        'logs': os.path.join(root_dir, 'logs'),
        'main_logs': os.path.join(root_dir, 'logs', 'main'),
        'trimming_logs': os.path.join(root_dir, 'logs', 'trimming'),
        'merging_logs': os.path.join(root_dir, 'logs', 'merging'),
        'normalization_logs': os.path.join(root_dir, 'logs', 'normalization'),
        'assembly_logs': os.path.join(root_dir, 'logs', 'assembly'),
        'busco_logs': os.path.join(root_dir, 'logs', 'busco')
    }
    
    # Create directories if they don't exist
    for directory in dirs.values():
        os.makedirs(directory, exist_ok=True)
    
    return dirs

def find_paired_reads(raw_reads_dir):
    """Identify paired read files and report unpaired ones."""
    # Update patterns to match compressed fastq files with various naming conventions
    r1_patterns = ['*_R1.fastq.gz', '*_R1_001.fastq.gz']
    r2_patterns = ['*_R2.fastq.gz', '*_R2_001.fastq.gz']
    
    print(f"Searching for files in: {raw_reads_dir}")
    
    # Collect all files matching any R1 pattern
    r1_files = []
    for pattern in r1_patterns:
        matching_files = glob.glob(os.path.join(raw_reads_dir, pattern))
        print(f"Pattern '{pattern}' matched {len(matching_files)} files")
        r1_files.extend(matching_files)
    
    # Collect all files matching any R2 pattern
    r2_files = []
    for pattern in r2_patterns:
        r2_files.extend(glob.glob(os.path.join(raw_reads_dir, pattern)))
    
    # Extract sample names by removing the R1/R2 part and extensions
    r1_samples = set()
    for f in r1_files:
        basename = os.path.basename(f)
        sample = basename.replace('_R1_001.fastq.gz', '').replace('_R1.fastq.gz', '')
        r1_samples.add(sample)
    
    r2_samples = set()
    for f in r2_files:
        basename = os.path.basename(f)
        sample = basename.replace('_R2_001.fastq.gz', '').replace('_R2.fastq.gz', '')
        r2_samples.add(sample)
    
    paired_samples = r1_samples & r2_samples
    unpaired_r1 = r1_samples - r2_samples
    unpaired_r2 = r2_samples - r1_samples
    
    print(f"Found {len(paired_samples)} paired samples")
    if unpaired_r1:
        print(f"Unpaired R1 files: {', '.join(unpaired_r1)}")
    if unpaired_r2:
        print(f"Unpaired R2 files: {', '.join(unpaired_r2)}")
    
    # Create pairs with correct paths
    pairs = []
    for sample in paired_samples:
        # Try to find the matching files for this sample
        r1_file = None
        r2_file = None
        
        # Check for _R1_001.fastq.gz pattern first
        r1_candidate = os.path.join(raw_reads_dir, f"{sample}_R1_001.fastq.gz")
        if os.path.exists(r1_candidate):
            r1_file = r1_candidate
            r2_file = os.path.join(raw_reads_dir, f"{sample}_R2_001.fastq.gz")
        else:
            # Fall back to _R1.fastq.gz pattern
            r1_file = os.path.join(raw_reads_dir, f"{sample}_R1.fastq.gz")
            r2_file = os.path.join(raw_reads_dir, f"{sample}_R2.fastq.gz")
        
        if os.path.exists(r1_file) and os.path.exists(r2_file):
            pairs.append((sample, r1_file, r2_file))
    
    return pairs

def all_trimmed_files_exist(pairs, trimmed_reads_dir):
    """Check if all trimmed files exist and are non-empty."""
    for sample, _, _ in pairs:
        out_r1 = os.path.join(trimmed_reads_dir, f"{sample}_trimmed_R1.fastq.gz")
        out_r2 = os.path.join(trimmed_reads_dir, f"{sample}_trimmed_R2.fastq.gz")
        if not (os.path.exists(out_r1) and os.path.getsize(out_r1) > 0 and
                os.path.exists(out_r2) and os.path.getsize(out_r2) > 0):
            return False
    return True

def submit_trimming_jobs(pairs, dirs):
    """Submit trimming jobs for each pair using fastp."""
    job_ids = []
    for sample, r1, r2 in pairs:
        sbatch_script = f"""#!/bin/bash
#SBATCH --partition=day-long-cpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=fastp_{sample}
#SBATCH --output={dirs['trimming_logs']}/{sample}_fastp.out
#SBATCH --error={dirs['trimming_logs']}/{sample}_fastp.err

# Source bashrc and activate conda environment
source ~/.bashrc
conda activate transcriptome

fastp \\
    -i "{r1}" \\
    -I "{r2}" \\
    -o "{dirs['trimmed_reads']}/{sample}_trimmed_R1.fastq.gz" \\
    -O "{dirs['trimmed_reads']}/{sample}_trimmed_R2.fastq.gz" \\
    --html "{dirs['trimming_results']}/{sample}_fastp.html" \\
    --json "{dirs['trimming_results']}/{sample}_fastp.json" \\
    --detect_adapter_for_pe \\
    --dedup \\
    -w 4 \\
    -q 20 \\
    --compression 6
"""
        result = subprocess.run(['sbatch'], input=sbatch_script, text=True, capture_output=True)
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            print(f"Submitted trimming job for {sample}: job ID {job_id}")
        else:
            print(f"Failed to submit trimming job for {sample}: {result.stderr}")
            sys.exit(1)
    return job_ids

def submit_normalization_job(dirs, dependency=None):
    """Submit job to normalize reads using BBNorm with repair step."""
    sbatch_script = f"""#!/bin/bash
#SBATCH --partition=week-long-highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --job-name=bbnorm
#SBATCH --output={dirs['normalization_logs']}/bbnorm.out
#SBATCH --error={dirs['normalization_logs']}/bbnorm.err
"""
    if dependency:
        sbatch_script += f"#SBATCH --dependency=afterok:{':'.join(dependency)}\n"
    
    sbatch_script += f"""
# Source bashrc and activate conda environment
source ~/.bashrc
conda activate transcriptome

# Define input and output files
in1="{dirs['normalization_results']}/left.fa"
in2="{dirs['normalization_results']}/right.fa"
out1="{dirs['normalization_results']}/left.norm.fa"
out2="{dirs['normalization_results']}/right.norm.fa"
repaired_left="{dirs['normalization_results']}/repaired_left.fa"
repaired_right="{dirs['normalization_results']}/repaired_right.fa"
singletons="{dirs['normalization_results']}/singletons.fa"
hist_in="{dirs['normalization_results']}/histogram_in.txt"
hist_out="{dirs['normalization_results']}/histogram_out.txt"
peaks="{dirs['normalization_results']}/peaks.txt"
target=100

# Create input FASTA files if they don't exist
if [ ! -f "$in1" ] || [ ! -f "$in2" ] || [ ! -s "$in1" ] || [ ! -s "$in2" ]; then
    echo "Merged FASTA files don't exist or are empty. Creating them now..."
    
    # Remove any existing but empty files
    rm -f "$in1" "$in2"
    
    # Ensure the directory exists
    mkdir -p "{dirs['normalization_results']}"
    
    # Convert and merge R1 files
    echo "Converting and merging R1 files at $(date)"
    for R1_FILE in {dirs['trimmed_reads']}/*_trimmed_R1.fastq.gz; do
        if [ -f "$R1_FILE" ]; then
            BASENAME=$(basename "$R1_FILE" _trimmed_R1.fastq.gz)
            echo "Processing $R1_FILE"
            pigz -dc "$R1_FILE" | seqtk seq -A - >> "$in1"
            if [ $? -ne 0 ]; then
                echo "Error processing $R1_FILE"
                exit 1
            fi
        fi
    done
    
    # Convert and merge R2 files
    echo "Converting and merging R2 files at $(date)"
    for R2_FILE in {dirs['trimmed_reads']}/*_trimmed_R2.fastq.gz; do
        if [ -f "$R2_FILE" ]; then
            BASENAME=$(basename "$R2_FILE" _trimmed_R2.fastq.gz)
            echo "Processing $R2_FILE"
            pigz -dc "$R2_FILE" | seqtk seq -A - >> "$in2"
            if [ $? -ne 0 ]; then
                echo "Error processing $R2_FILE"
                exit 1
            fi
        fi
    done
    
    # Check if files were created successfully
    if [ ! -s "$in1" ] || [ ! -s "$in2" ]; then
        echo "Failed to create merged FASTA files. Check if trimmed reads exist."
        echo "Trimmed reads directory: {dirs['trimmed_reads']}"
        ls -la {dirs['trimmed_reads']}
        exit 1
    fi
    
    echo "Successfully created merged FASTA files:"
    du -h "$in1" "$in2"
fi

# Now check if input files exist and are non-empty
if [ ! -s "$in1" ] || [ ! -s "$in2" ]; then
    echo "Error: Input files $in1 or $in2 do not exist or are empty."
    exit 1
fi

# Step 1: Repair the input reads to ensure pairing
echo "Starting repair.sh at $(date)"
repair.sh in1="$in1" in2="$in2" out1="$repaired_left" out2="$repaired_right" outs="$singletons" \\
          threads=64 -Xmx200g verbose=t overwrite=t
echo "Finished repair.sh at $(date)"

# Verify the number of reads in the repaired files
left_count=$(grep -c '^>' "$repaired_left")
right_count=$(grep -c '^>' "$repaired_right")
singleton_count=$(grep -c '^>' "$singletons")

echo "Repaired left reads: $left_count"
echo "Repaired right reads: $right_count"
echo "Singletons: $singleton_count"

if [ "$left_count" -ne "$right_count" ]; then
    echo "Error: Mismatch in read counts between $repaired_left ($left_count) and $repaired_right ($right_count)"
    exit 1
fi

if [ "$singleton_count" -gt 0 ]; then
    echo "Warning: $singleton_count singletons found. Check $singletons for unpaired reads."
fi

# Step 2: Normalize the repaired reads, capping kmer depth at target value
echo "Starting bbnorm.sh at $(date)"
bbnorm.sh in1="$repaired_left" in2="$repaired_right" out1="$out1" out2="$out2" \\
          target="$target" mindepth=5 maxdepth=100 passes=2 \\
          k=31 prefilter=t prefiltersize=0.35 prehashes=2 prefilterbits=2 buildpasses=1 \\
          bits=16 hashes=3 \\
          threads=64 interleaved=false \\
          ecc=f tossbadreads=f fixspikes=t deterministic=t \\
          hist="$hist_in" histout="$hist_out" peaks="$peaks" \\
          zerobin=t pzc=t histlen=10000 \\
          minq=6 minprob=0.5 \\
          -Xmx240g -eoom -da overwrite=t
echo "Finished bbnorm.sh at $(date)"

# Check output files and verify depth cap
if [ ! -f "$out1" ] || [ ! -f "$out2" ]; then
    echo "Error: Normalization failed. Output files $out1 or $out2 not generated."
    exit 1
fi

echo "Normalization completed successfully at $(date)"
echo "Kmer depth capped at $target. Check $hist_out to confirm depth distribution."

# Get sizes of normalized files for reporting
norm_left_size=$(stat -c %s "$out1")
norm_right_size=$(stat -c %s "$out2")
echo "Normalized left file size: $(numfmt --to=iec-i --suffix=B $norm_left_size)"
echo "Normalized right file size: $(numfmt --to=iec-i --suffix=B $norm_right_size)"
"""
    result = subprocess.run(['sbatch'], input=sbatch_script, text=True, capture_output=True)
    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]
        print(f"Submitted BBNorm normalization job: job ID {job_id}")
        return job_id
    else:
        print(f"Failed to submit BBNorm normalization job: {result.stderr}")
        sys.exit(1)

def normalized_files_exist(normalization_results_dir):
    """Check if normalized FASTA files exist and are non-empty."""
    norm_left = os.path.join(normalization_results_dir, 'left.norm.fa')
    norm_right = os.path.join(normalization_results_dir, 'right.norm.fa')
    return (os.path.exists(norm_left) and os.path.getsize(norm_left) > 0 and
            os.path.exists(norm_right) and os.path.getsize(norm_right) > 0)

def submit_assembly_jobs(dirs, dependency=None):
    """Submit assembly jobs for rnaSPAdes and Trinity with FASTA input."""
    job_ids = {}
    
    # rnaSPAdes
    rnaspades_script = f"""#!/bin/bash
#SBATCH --partition=week-long-highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --job-name=rnaspades
#SBATCH --output={dirs['assembly_logs']}/rnaspades.out
#SBATCH --error={dirs['assembly_logs']}/rnaspades.err
"""
    if dependency:
        rnaspades_script += f"#SBATCH --dependency=afterok:{dependency}\n"
    
    rnaspades_script += f"""
# Source bashrc and activate conda environment
source ~/.bashrc
conda activate transcriptome

rnaspades.py \\
    --rna \\
    -1 "{dirs['normalization_results']}/left.norm.fa" \\
    -2 "{dirs['normalization_results']}/right.norm.fa" \\
    -o "{dirs['assembly_results']}/rnaspades" \\
    -t $SLURM_CPUS_PER_TASK \\
    -m 250 \\
    2>> {dirs['assembly_logs']}/rnaspades.log
"""
    result = subprocess.run(['sbatch'], input=rnaspades_script, text=True, capture_output=True)
    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]
        job_ids['rnaspades'] = job_id
        print(f"Submitted rnaSPAdes job: job ID {job_id}")
    else:
        print(f"Failed to submit rnaSPAdes job: {result.stderr}")
        sys.exit(1)
    
    # Trinity
    trinity_script = f"""#!/bin/bash
#SBATCH --partition=week-long-highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --job-name=trinity
#SBATCH --output={dirs['assembly_logs']}/trinity.out
#SBATCH --error={dirs['assembly_logs']}/trinity.err
"""
    if dependency:
        trinity_script += f"#SBATCH --dependency=afterok:{dependency}\n"
    
    trinity_script += f"""
# Source bashrc and activate conda environment
source ~/.bashrc
conda activate transcriptome

Trinity \\
    --seqType fa \\
    --left "{dirs['normalization_results']}/left.norm.fa" \\
    --right "{dirs['normalization_results']}/right.norm.fa" \\
    --no_normalize \\
    --CPU $SLURM_CPUS_PER_TASK \\
    --max_memory 250G \\
    --output "{dirs['assembly_results']}/trinity" \\
    2>> {dirs['assembly_logs']}/trinity.log
"""
    result = subprocess.run(['sbatch'], input=trinity_script, text=True, capture_output=True)
    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]
        job_ids['trinity'] = job_id
        print(f"Submitted Trinity job: job ID {job_id}")
    else:
        print(f"Failed to submit Trinity job: {result.stderr}")
        sys.exit(1)
    
    return job_ids

def assembly_files_exist(assembly_results_dir, assembler):
    """Check if assembly output files exist and are non-empty."""
    if assembler == 'rnaspades':
        # Check for the main transcripts.fasta file from rnaSPAdes
        transcripts_file = os.path.join(assembly_results_dir, 'rnaspades', 'transcripts.fasta')
        return os.path.exists(transcripts_file) and os.path.getsize(transcripts_file) > 0
    elif assembler == 'trinity':
        # Check for the Trinity.fasta file from Trinity
        trinity_file = os.path.join(assembly_results_dir, 'trinity', 'Trinity.fasta')
        return os.path.exists(trinity_file) and os.path.getsize(trinity_file) > 0
    return False

def busco_files_exist(busco_results_dir, assembler):
    """Check if BUSCO output files exist and are non-empty."""
    # Check for the short_summary file which indicates BUSCO completed
    summary_file = os.path.join(busco_results_dir, assembler, 'short_summary.*.txt')
    summary_files = glob.glob(summary_file)
    return len(summary_files) > 0

def submit_busco_jobs(dirs, assembly_job_ids):
    """Submit BUSCO analysis jobs for each assembler's output."""
    busco_job_ids = {}
    
    for assembler, job_id in assembly_job_ids.items():
        # Define input file path based on assembler
        if assembler == 'rnaspades':
            input_fasta = f"{dirs['assembly_results']}/rnaspades/transcripts.fasta"
        elif assembler == 'trinity':
            input_fasta = f"{dirs['assembly_results']}/trinity/Trinity.fasta"
        else:
            continue
        
        busco_script = f"""#!/bin/bash
#SBATCH --partition=day-long-cpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --job-name=busco_{assembler}
#SBATCH --output={dirs['busco_logs']}/{assembler}_busco.out
#SBATCH --error={dirs['busco_logs']}/{assembler}_busco.err
"""
        if job_id:
            busco_script += f"#SBATCH --dependency=afterok:{job_id}\n"
        
        busco_script += f"""
# Source bashrc and activate conda environment
source ~/.bashrc
conda activate transcriptome

# Create output directory
mkdir -p {dirs['busco_results']}/{assembler}

busco \\
    -i "{input_fasta}" \\
    -o "{assembler}" \\
    -l eukaryota_odb10 \\
    -m transcriptome \\
    -c $SLURM_CPUS_PER_TASK \\
    --out_path "{dirs['busco_results']}" \\
    2>> {dirs['busco_logs']}/{assembler}_busco.log
"""
        result = subprocess.run(['sbatch'], input=busco_script, text=True, capture_output=True)
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            busco_job_ids[assembler] = job_id
            print(f"Submitted BUSCO job for {assembler}: job ID {job_id}")
        else:
            print(f"Failed to submit BUSCO job for {assembler}: {result.stderr}")
            sys.exit(1)
    
    return busco_job_ids

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run transcriptome assembly pipeline')
    parser.add_argument('-i', '--input-dir', help='Path to the directory containing raw reads')
    args = parser.parse_args()
    
    # Setup directories
    dirs = setup_directories()
    
    # If input directory is specified, use it instead of the default raw_reads directory
    if args.input_dir:
        if os.path.isdir(args.input_dir):
            raw_reads_dir = args.input_dir
            print(f"Using raw reads from: {raw_reads_dir}")
        else:
            print(f"Error: Input directory {args.input_dir} does not exist.")
            sys.exit(1)
    else:
        raw_reads_dir = dirs['raw_reads']
        print(f"Using default raw reads directory: {raw_reads_dir}")
    
    # Find read pairs
    pairs = find_paired_reads(raw_reads_dir)
    
    if not pairs:
        print("No read pairs found. Exiting.")
        sys.exit(1)
    
    print(f"Found {len(pairs)} paired read files.")
    
    # Check if all trimmed files exist
    if all_trimmed_files_exist(pairs, dirs['trimmed_reads']):
        print("All trimmed files already exist, skipping trimming step.")
        trimming_job_ids = []
    else:
        # Submit trimming jobs
        print("Submitting trimming jobs...")
        trimming_job_ids = submit_trimming_jobs(pairs, dirs)
    
    # Check if normalized files exist
    if normalized_files_exist(dirs['normalization_results']):
        print("Normalized files already exist, skipping normalization step.")
        normalization_job_id = None
    else:
        # Submit normalization job
        print("Submitting normalization job...")
        normalization_job_id = submit_normalization_job(dirs, trimming_job_ids)
    
    # Check if assembly files exist
    rnaspades_exists = assembly_files_exist(dirs['assembly_results'], 'rnaspades')
    trinity_exists = assembly_files_exist(dirs['assembly_results'], 'trinity')
    
    if rnaspades_exists and trinity_exists:
        print("All assembly files already exist, skipping assembly step.")
        assembly_job_ids = {}
    else:
        # Submit assembly jobs
        print("Submitting assembly jobs...")
        assembly_job_ids = submit_assembly_jobs(dirs, normalization_job_id)
    
    # Check if BUSCO results exist
    rnaspades_busco_exists = busco_files_exist(dirs['busco_results'], 'rnaspades')
    trinity_busco_exists = busco_files_exist(dirs['busco_results'], 'trinity')
    
    if rnaspades_busco_exists and trinity_busco_exists:
        print("All BUSCO results already exist, skipping BUSCO step.")
    else:
        # Submit BUSCO jobs
        print("Submitting BUSCO jobs...")
        busco_job_ids = submit_busco_jobs(dirs, assembly_job_ids)
    
    print("Pipeline submitted successfully. Monitor jobs with 'squeue -u $USER'.")

if __name__ == "__main__":
    main()