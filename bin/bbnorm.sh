#!/bin/bash
#SBATCH --partition=week-long-highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --job-name=trinity_norm
#SBATCH --output=%j_bbnorm.out
#SBATCH --error=%j_bbnorm.err

# Activate the environment
source ~/.bashrc
conda activate transcriptome

# Define input and output files
in1=/nfs/home/jlamb/Projects/nf_denovo_transcriptome/results/merging/merged_R1.fastq.gz
in2=/nfs/home/jlamb/Projects/nf_denovo_transcriptome/results/merging/merged_R2.fastq.gz
out1=/nfs/home/jlamb/Projects/transcriptome/results/normalization/left.norm.fa
out2=/nfs/home/jlamb/Projects/transcriptome/results/normalization/right.norm.fa
repaired_left=/nfs/home/jlamb/Projects/transcriptome/results/normalization/repaired_left.fa
repaired_right=/nfs/home/jlamb/Projects/transcriptome/results/normalization/repaired_right.fa
singletons=/nfs/home/jlamb/Projects/transcriptome/results/normalization/singletons.fa
hist_in=/nfs/home/jlamb/Projects/transcriptome/results/normalization/histogram_in.txt
hist_out=/nfs/home/jlamb/Projects/transcriptome/results/normalization/histogram_out.txt
peaks=/nfs/home/jlamb/Projects/transcriptome/results/normalization/peaks.txt
target=100

# Check if input files exist
if [ ! -f "$in1" ] || [ ! -f "$in2" ]; then
    echo "Error: Input files $in1 or $in2 do not exist."
    exit 1
fi

# Step 1: Repair the input reads to ensure pairing
echo "Starting repair.sh at $(date)"
repair.sh in1="$in1" in2="$in2" out1="$repaired_left" out2="$repaired_right" outs="$singletons" \
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

# Step 2: Normalize the repaired reads, capping kmer depth at 100
echo "Starting bbnorm.sh at $(date)"
bbnorm.sh in1="$repaired_left" in2="$repaired_right" out1="$out1" out2="$out2" \
          target="$target" mindepth=5 maxdepth=100 passes=2 \
          k=31 prefilter=t prefiltersize=0.35 prehashes=2 prefilterbits=2 buildpasses=1 \
          bits=16 hashes=3 \
          threads=64 interleaved=false \
          ecc=f tossbadreads=f fixspikes=t deterministic=t \
          hist="$hist_in" histout="$hist_out" peaks="$peaks" \
          zerobin=t pzc=t histlen=10000 \
          minq=6 minprob=0.5 \
	  tossbrokenreads=t \
          -Xmx240g -eoom -da overwrite=t
echo "Finished bbnorm.sh at $(date)"

# Check output files and verify depth cap
if [ ! -f "$out1" ] || [ ! -f "$out2" ]; then
    echo "Error: Normalization failed. Output files $out1 or $out2 not generated."
    exit 1
fi

echo "Normalization completed successfully at $(date)"
echo "Kmer depth capped at 100. Check $hist_out to confirm depth distribution."
