#!/usr/bin/env python3

import os
import glob
import json
import subprocess
import re
from datetime import datetime
import shutil

def format_size(size):
    """Convert size in bytes to a human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return f"{size:.2f} {unit}"
        size /= 1024.0
    return f"{size:.2f} PB"

def get_file_sizes(directory, pattern):
    """Get list of (filename, size) for files matching the pattern in the directory."""
    files = glob.glob(os.path.join(directory, pattern))
    sizes = [(os.path.basename(f), os.path.getsize(f)) for f in files if os.path.exists(f)]
    return sizes

def get_fastp_params(json_file):
    """Extract fastp command from the JSON report."""
    if not os.path.exists(json_file):
        return "N/A"
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data.get('command', 'N/A')

def submit_kmer_job(fasta_file, output_dir, k=25):
    """Submit a Slurm job to compute k-mer statistics using BBTools kmercountexact.sh."""
    if not os.path.exists(fasta_file):
        return {'unique_kmers': 'N/A', 'avg_coverage': 'N/A', 'job_id': 'N/A'}

    # Check if kmercountexact.sh is available
    if not shutil.which('kmercountexact.sh'):
        return {'unique_kmers': 'BBTools not found', 'avg_coverage': 'N/A', 'job_id': 'N/A'}

    base_name = os.path.basename(fasta_file).replace('.fa', '')
    histo_output = os.path.join(output_dir, f"{base_name}_k{k}_histogram.txt")
    log_file = os.path.join(output_dir, f"{base_name}_kmer.log")

    # Slurm script for k-mer counting with BBTools
    sbatch_script = f"""#!/bin/bash
#SBATCH --partition=week-long-highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --job-name=kmer_{base_name}
#SBATCH --output={log_file}
#SBATCH --error={log_file}

source ~/.bashrc
conda activate transcriptome

# Run kmercountexact.sh
kmercountexact.sh \\
    in={fasta_file} \\
    k={k} \\
    threads=64 \\
    khist={histo_output} \\
    2>> {log_file}
"""

    # Submit the job
    result = subprocess.run(['sbatch'], input=sbatch_script, text=True, capture_output=True)
    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]
        print(f"Submitted k-mer job for {base_name}: job ID {job_id}")
        return {'unique_kmers': 'Pending', 'avg_coverage': 'Pending', 'job_id': job_id}
    else:
        print(f"Failed to submit k-mer job for {base_name}: {result.stderr}")
        return {'unique_kmers': 'Failed', 'avg_coverage': 'N/A', 'job_id': 'N/A'}

def parse_kmer_stats(histo_file):
    """Parse BBTools histogram to compute unique k-mers and average coverage."""
    if not os.path.exists(histo_file):
        return {'unique_kmers': 'N/A', 'avg_coverage': 'N/A'}
    
    unique_kmers = 0
    total_kmers = 0
    with open(histo_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                count, freq = int(parts[0]), int(parts[1])
                unique_kmers += 1
                total_kmers += count * freq
    
    avg_coverage = total_kmers / unique_kmers if unique_kmers > 0 else 0
    return {
        'unique_kmers': unique_kmers,
        'avg_coverage': f"{avg_coverage:.2f}"
    }

def main():
    # Define base directory relative to script location
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dirs = {
        'raw_reads': os.path.join(base_dir, 'data', 'raw_reads'),
        'trimmed_reads': os.path.join(base_dir, 'data', 'trimmed_reads'),
        'trimming_results': os.path.join(base_dir, 'results', 'trimming'),
        'normalization_results': os.path.join(base_dir, 'results', 'normalization'),
        'assembly_results': os.path.join(base_dir, 'results', 'assembly'),
        'statistics': os.path.join(base_dir, 'results', 'statistics')
    }

    # Create statistics directory if it doesn't exist
    os.makedirs(dirs['statistics'], exist_ok=True)
    report_file = os.path.join(dirs['statistics'], 'report.txt')

    with open(report_file, 'w') as report:
        # Header
        report.write("Transcriptome Analysis Statistics Report\n")
        report.write("=======================================\n")
        report.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Raw Read Files
        report.write("Raw Read Files:\n")
        report.write("---------------\n")
        raw_files = get_file_sizes(dirs['raw_reads'], '*_R[12].fastq.gz')
        if raw_files:
            for name, size in sorted(raw_files):
                report.write(f"{name}: {format_size(size)}\n")
        else:
            report.write("No raw read files found.\n")
        report.write("\n")

        # Trimmed Read Files
        report.write("Trimmed Read Files:\n")
        report.write("-------------------\n")
        trimmed_files = get_file_sizes(dirs['trimmed_reads'], '*_trimmed_R[12].fastq.gz')
        if trimmed_files:
            for name, size in sorted(trimmed_files):
                report.write(f"{name}: {format_size(size)}\n")
        else:
            report.write("No trimmed read files found.\n")
        report.write("\n")

        # fastp Parameters
        report.write("fastp Parameters:\n")
        report.write("-----------------\n")
        json_files = glob.glob(os.path.join(dirs['trimming_results'], '*_fastp.json'))
        if json_files:
            params = get_fastp_params(json_files[0])
            report.write(f"{params}\n")
        else:
            report.write("No fastp JSON files found.\n")
        report.write("\n")

        # Merged FASTA Files
        report.write("Merged FASTA Files:\n")
        report.write("-------------------\n")
        merged_files = [
            ('left.fa', os.path.join(dirs['normalization_results'], 'left.fa')),
            ('right.fa', os.path.join(dirs['normalization_results'], 'right.fa'))
        ]
        for name, path in merged_files:
            if os.path.exists(path):
                size = os.path.getsize(path)
                report.write(f"{name}: {format_size(size)}\n")
            else:
                report.write(f"{name}: Not found\n")
        report.write("\n")

        # Normalized FASTA Files
        report.write("Normalized FASTA Files:\n")
        report.write("-----------------------\n")
        norm_files = [
            ('left.norm.fa', os.path.join(dirs['normalization_results'], 'left.norm.fa')),
            ('right.norm.fa', os.path.join(dirs['normalization_results'], 'right.norm.fa'))
        ]
        for name, path in norm_files:
            if os.path.exists(path):
                size = os.path.getsize(path)
                report.write(f"{name}: {format_size(size)}\n")
            else:
                report.write(f"{name}: Not found\n")
        report.write("\n")

        # K-mer Statistics
        report.write("K-mer Statistics (k=25):\n")
        report.write("------------------------\n")
        kmer_jobs = {}
        
        report.write("Before normalization:\n")
        for name, path in merged_files:
            report.write(f"  {name}:\n")
            stats = submit_kmer_job(path, dirs['statistics'])
            kmer_jobs[f"{name}_k25_histogram.txt"] = stats['job_id']
            report.write(f"    Unique k-mers: {stats['unique_kmers']}\n")
            report.write(f"    Average k-mer coverage: {stats['avg_coverage']}\n")
            if stats['job_id'] != 'N/A':
                report.write(f"    Job ID: {stats['job_id']} (Check {dirs['statistics']}/{name}_k25_histogram.txt after completion)\n")
        
        report.write("After normalization:\n")
        for name, path in norm_files:
            report.write(f"  {name}:\n")
            stats = submit_kmer_job(path, dirs['statistics'])
            kmer_jobs[f"{name}_k25_histogram.txt"] = stats['job_id']
            report.write(f"    Unique k-mers: {stats['unique_kmers']}\n")
            report.write(f"    Average k-mer coverage: {stats['avg_coverage']}\n")
            if stats['job_id'] != 'N/A':
                report.write(f"    Job ID: {stats['job_id']} (Check {dirs['statistics']}/{name}_k25_histogram.txt after completion)\n")
        report.write("\n")

        # Assembly Files
        report.write("Assembly Files:\n")
        report.write("---------------\n")
        assembly_files = [
            ('rnaspades/transcripts.fasta', os.path.join(dirs['assembly_results'], 'rnaspades', 'transcripts.fasta')),
            ('trinity/Trinity.fasta', os.path.join(dirs['assembly_results'], 'trinity', 'Trinity.fasta'))
        ]
        for name, path in assembly_files:
            if os.path.exists(path):
                size = os.path.getsize(path)
                report.write(f"{name}: {format_size(size)}\n")
            else:
                report.write(f"{name}: Not found\n")
        report.write("\n")

        # Instructions for post-processing
        report.write("Post-Processing Instructions:\n")
        report.write("-----------------------------\n")
        report.write("K-mer jobs have been submitted. After they complete, re-run this script to parse the histograms:\n")
        report.write(f"  python {os.path.abspath(__file__)} --parse-only\n")
        report.write("This will update the report with final k-mer statistics.\n")

    # Optional: Parse existing histograms if --parse-only flag is provided
    if '--parse-only' in os.sys.argv:
        with open(report_file, 'a') as report:
            report.write("\nUpdated K-mer Statistics (k=25):\n")
            report.write("--------------------------------\n")
            report.write("Before normalization:\n")
            for name, path in merged_files:
                histo_file = os.path.join(dirs['statistics'], f"{name}_k25_histogram.txt")
                report.write(f"  {name}:\n")
                stats = parse_kmer_stats(histo_file)
                report.write(f"    Unique k-mers: {stats['unique_kmers']}\n")
                report.write(f"    Average k-mer coverage: {stats['avg_coverage']}\n")
            report.write("After normalization:\n")
            for name, path in norm_files:
                histo_file = os.path.join(dirs['statistics'], f"{name}_k25_histogram.txt")
                report.write(f"  {name}:\n")
                stats = parse_kmer_stats(histo_file)
                report.write(f"    Unique k-mers: {stats['unique_kmers']}\n")
                report.write(f"    Average k-mer coverage: {stats['avg_coverage']}\n")

    print(f"Statistics report generated: {report_file}")
    if '--parse-only' not in os.sys.argv:
        print("K-mer jobs submitted. Monitor with 'squeue -u $USER'. Re-run with '--parse-only' after completion.")

if __name__ == "__main__":
    main()