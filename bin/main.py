#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import glob

def setup_directories(base_dir):
    """Define and create directory structure if not exists."""
    dirs = {
        'data': os.path.join(base_dir, 'data'),
        'raw_reads': os.path.join(base_dir, 'data', 'raw_reads'),
        'trimmed_reads': os.path.join(base_dir, 'data', 'trimmed_reads'),
        'transcriptome': os.path.join(base_dir, 'data', 'transcriptome'),
        'results': os.path.join(base_dir, 'results'),
        'trimming_results': os.path.join(base_dir, 'results', 'trimming'),
        'merging_results': os.path.join(base_dir, 'results', 'merging'),
        'normalization_results': os.path.join(base_dir, 'results', 'normalization'),
        'assembly_results': os.path.join(base_dir, 'results', 'assembly'),
        'busco_results': os.path.join(base_dir, 'results', 'busco'),
        'logs': os.path.join(base_dir, 'logs'),
        'main_logs': os.path.join(base_dir, 'logs', 'main'),
        'trimming_logs': os.path.join(base_dir, 'logs', 'trimming'),
        'merging_logs': os.path.join(base_dir, 'logs', 'merging'),
        'normalization_logs': os.path.join(base_dir, 'logs', 'normalization'),
        'assembly_logs': os.path.join(base_dir, 'logs', 'assembly'),
        'busco_logs': os.path.join(base_dir, 'logs', 'busco')
    }
    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)
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
    """Submit job to normalize reads using BBNorm instead of Trinity."""
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

# Check if trimmed files exist
if [ -z "$(ls {dirs['trimmed_reads']}/*_trimmed_R1.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No trimmed R1 files found in {dirs['trimmed_reads']}" >&2
    exit 1
fi
if [ -z "$(ls {dirs['trimmed_reads']}/*_trimmed_R2.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No trimmed R2 files found in {dirs['trimmed_reads']}" >&2
    exit 1
fi

# Create directories for intermediate files
mkdir -p {dirs['normalization_results']}/tmp_fasta

# Prepare input files
echo "Preparing R1 FASTA files"
for R1_FILE in {dirs['trimmed_reads']}/*_trimmed_R1.fastq.gz; do
    BASENAME=$(basename "$R1_FILE" _trimmed_R1.fastq.gz)
    echo "Converting $R1_FILE to FASTA"
    pigz -dc "$R1_FILE" | seqtk seq -A - > {dirs['normalization_results']}/tmp_fasta/"$BASENAME"_R1.fa
done

echo "Preparing R2 FASTA files"
for R2_FILE in {dirs['trimmed_reads']}/*_trimmed_R2.fastq.gz; do
    BASENAME=$(basename "$R2_FILE" _trimmed_R2.fastq.gz)
    echo "Converting $R2_FILE to FASTA"
    pigz -dc "$R2_FILE" | seqtk seq -A - > {dirs['normalization_results']}/tmp_fasta/"$BASENAME"_R2.fa
done

# Combine all FASTA files
echo "Combining all R1 files"
cat {dirs['normalization_results']}/tmp_fasta/*_R1.fa > {dirs['normalization_results']}/left.fa
echo "Combining all R2 files"
cat {dirs['normalization_results']}/tmp_fasta/*_R2.fa > {dirs['normalization_results']}/right.fa

# Define input and output files for BBNorm
in1="{dirs['normalization_results']}/left.fa"
in2="{dirs['normalization_results']}/right.fa"
out1="{dirs['normalization_results']}/left.norm.fa"
out2="{dirs['normalization_results']}/right.norm.fa"
target=100

# Run BBNorm for normalization
echo "Starting BBNorm normalization"
bbnorm.sh in1="$in1" in2="$in2" out1="$out1" out2="$out2" \\
          target="$target" mindepth=3 maxdepth=-1 passes=2 \\
          k=25 prefilter=t prefiltersize=0.4 buildpasses=2 bits=32 \\
          threads=64 interleaved=false \\
          ecc=f tossbadreads=f \\
          -Xmx240g \\
          2>> "{dirs['normalization_logs']}/bbnorm.log"

# Clean up temporary files if normalization was successful
if [ -f "$out1" ] && [ -f "$out2" ]; then
    echo "Normalization successful, cleaning up temporary files"
    rm -rf {dirs['normalization_results']}/tmp_fasta
else
    echo "Normalization failed, keeping temporary files for debugging"
fi
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
    # Updated to check .fa extensions for BBNorm output
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
    """Orchestrate the optimized transcriptome analysis pipeline."""
    parser = argparse.ArgumentParser(description='Transcriptome assembly pipeline')
    parser.add_argument('-d', '--debug', action='store_true', help='Skip completed steps in debug mode')
    parser.add_argument('-b', '--base-dir', default=None, 
                        help='Base directory for the project (defaults to parent directory of script)')
    args = parser.parse_args()
    
    if args.base_dir:
        base_dir = os.path.abspath(args.base_dir)
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.dirname(script_dir)
    
    print(f"Using base directory: {base_dir}")
    dirs = setup_directories(base_dir)
    
    print(f"Looking for reads in: {dirs['raw_reads']}")
    
    # Step 1: Identify paired reads
    pairs = find_paired_reads(dirs['raw_reads'])
    if not pairs:
        print("No paired samples found. Exiting.")
        sys.exit(1)
    
    # Step 2: Trimming
    if args.debug and all_trimmed_files_exist(pairs, dirs['trimmed_reads']):
        print("All trimmed files exist, skipping trimming")
        trimming_job_ids = []
    else:
        trimming_job_ids = submit_trimming_jobs(pairs, dirs)
    
    # Step 3: Normalization (Merging integrated here)
    if args.debug and normalized_files_exist(dirs['normalization_results']):
        print("Normalized files exist, skipping normalization")
        normalization_job_id = None
    else:
        normalization_job_id = submit_normalization_job(dirs, trimming_job_ids if trimming_job_ids else None)
    
    # Step 4: Assembly - Fix the logic to submit both assemblies at once
    run_assembly = False
    for assembler in ['rnaspades', 'trinity']:
        if not (args.debug and assembly_files_exist(dirs['assembly_results'], assembler)):
            run_assembly = True
            break
    
    if run_assembly:
        assembly_job_ids = submit_assembly_jobs(dirs, normalization_job_id)
    else:
        print("All assembly outputs exist, skipping assembly")
        assembly_job_ids = {'rnaspades': None, 'trinity': None}
    
    # Step 5: BUSCO
    run_busco = {}
    for assembler in ['rnaspades', 'trinity']:
        if args.debug and busco_files_exist(dirs['busco_results'], assembler):
            print(f"BUSCO results for {assembler} exist, skipping BUSCO for {assembler}")
            run_busco[assembler] = None
        else:
            run_busco[assembler] = assembly_job_ids[assembler]
    
    if any(run_busco.values()) or not args.debug:
        submit_busco_jobs(dirs, run_busco)

if __name__ == "__main__":
    main()