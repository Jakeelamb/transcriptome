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

def compute_kmer_stats(fasta_file, output_dir, k=25):
    """Compute k-mer statistics using Jellyfish."""
    if not os.path.exists(fasta_file):
        return {'unique_kmers': 'N/A', 'avg_coverage': 'N/A'}

    # Check if Jellyfish is available
    if not shutil.which('jellyfish'):
        return {'unique_kmers': 'Jellyfish not found', 'avg_coverage': 'N/A'}

    # Define output files
    base_name = os.path.basename(fasta_file).replace('.fa', '')
    kmer_output = os.path.join(output_dir, f"{base_name}_k{k}.jf")
    histo_output = os.path.join(output_dir, f"{base_name}_k{k}.histo")

    # Run Jellyfish count
    count_cmd = [
        'jellyfish', 'count', '-m', str(k), '-s', '100M', '-t', '4',
        '-C', '-o', kmer_output, fasta_file
    ]
    subprocess.run(count_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Run Jellyfish histo
    histo_cmd = ['jellyfish', 'histo', kmer_output, '-o', histo_output]
    subprocess.run(histo_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Parse histogram to compute unique k-mers and average coverage
    unique_kmers = 0
    total_kmers = 0
    with open(histo_output, 'r') as f:
        for line in f:
            count, freq = map(int, line.strip().split())
            unique_kmers += 1
            total_kmers += count * freq

    # Clean up temporary files
    os.remove(kmer_output)
    os.remove(histo_output)

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
    kmer_temp_dir = os.path.join(dirs['statistics'], 'kmer_temp')
    os.makedirs(kmer_temp_dir, exist_ok=True)

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
        report.write("Before normalization:\n")
        for name, path in merged_files:
            report.write(f"  {name}:\n")
            stats = compute_kmer_stats(path, kmer_temp_dir)
            report.write(f"    Unique k-mers: {stats['unique_kmers']}\n")
            report.write(f"    Average k-mer coverage: {stats['avg_coverage']}\n")
        report.write("After normalization:\n")
        for name, path in norm_files:
            report.write(f"  {name}:\n")
            stats = compute_kmer_stats(path, kmer_temp_dir)
            report.write(f"    Unique k-mers: {stats['unique_kmers']}\n")
            report.write(f"    Average k-mer coverage: {stats['avg_coverage']}\n")
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

    # Clean up temporary k-mer directory
    shutil.rmtree(kmer_temp_dir)
    print(f"Statistics report generated: {report_file}")

if __name__ == "__main__":
    main()