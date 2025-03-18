#!/bin/bash
#SBATCH --partition=week-long-highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --job-name=subsample_norm
#SBATCH --output=%j_subsample_norm.out
#SBATCH --error=%j_subsample_norm.err

# Activate the environment
source ~/.bashrc
conda activate transcriptome

# Define input and output files
original_left=/nfs/home/jlamb/Projects/transcriptome/results/normalization/left.fa
original_right=/nfs/home/jlamb/Projects/transcriptome/results/normalization/right.fa
subsampled_left=/nfs/home/jlamb/Projects/transcriptome/results/normalization/subsample/left.subset.fa
subsampled_right=/nfs/home/jlamb/Projects/transcriptome/results/normalization/subsample/right.subset.fa
norm_out_left=/nfs/home/jlamb/Projects/transcriptome/results/normalization/subsample/left.subset.norm.fa
norm_out_right=/nfs/home/jlamb/Projects/transcriptome/results/normalization/subsample/right.subset.norm.fa
hist_in=/nfs/home/jlamb/Projects/transcriptome/results/normalization/subsample/histogram_in.txt
hist_out=/nfs/home/jlamb/Projects/transcriptome/results/normalization/subsample/histogram_out.txt
peaks=/nfs/home/jlamb/Projects/transcriptome/results/normalization/subsample/peaks.txt
target=100

# Create output directory if it doesn't exist
mkdir -p /nfs/home/jlamb/Projects/transcriptome/results/normalization/subsample

# Check if input files exist
if [ ! -f "$original_left" ] || [ ! -f "$original_right" ]; then
    echo "Error: Input files $original_left or $original_right do not exist."
    exit 1
fi

# Calculate the number of sequences needed for 100G output
# Since we don't know exact sequence lengths, we'll estimate based on file size
# and adjust with a scaling factor to target ~100G per file

# Get total size of original files in bytes
original_left_size=$(stat -c %s "$original_left")
original_right_size=$(stat -c %s "$original_right")

# Calculate the sampling fraction (aiming for 100G per file)
target_size=$((100 * 1024 * 1024 * 1024)) # 100G in bytes
left_fraction=$(echo "scale=6; $target_size / $original_left_size" | bc)
right_fraction=$(echo "scale=6; $target_size / $original_right_size" | bc)

echo "Original left file size: $(numfmt --to=iec-i --suffix=B $original_left_size)"
echo "Original right file size: $(numfmt --to=iec-i --suffix=B $original_right_size)"
echo "Sampling fraction for left: $left_fraction"
echo "Sampling fraction for right: $right_fraction"

# Step 1: Subsample the left and right files
echo "Starting seqtk subsampling at $(date)"
echo "Subsampling left file..."
seqtk sample -s 42 "$original_left" $left_fraction > "$subsampled_left"
echo "Subsampling right file..."
seqtk sample -s 42 "$original_right" $right_fraction > "$subsampled_right"
echo "Finished subsampling at $(date)"

# Check if subsampled files were created successfully
if [ ! -f "$subsampled_left" ] || [ ! -f "$subsampled_right" ]; then
    echo "Error: Subsampling failed. Output files $subsampled_left or $subsampled_right not generated."
    exit 1
fi

# Get sizes of subsampled files
subsampled_left_size=$(stat -c %s "$subsampled_left")
subsampled_right_size=$(stat -c %s "$subsampled_right")
echo "Subsampled left file size: $(numfmt --to=iec-i --suffix=B $subsampled_left_size)"
echo "Subsampled right file size: $(numfmt --to=iec-i --suffix=B $subsampled_right_size)"

# Step 2: Run bbnorm on the subsampled files
echo "Starting bbnorm.sh at $(date)"
bbnorm.sh in1="$subsampled_left" in2="$subsampled_right" out1="$norm_out_left" out2="$norm_out_right" \
          target="$target" mindepth=5 maxdepth=100 passes=2 \
          k=31 prefilter=t prefiltersize=0.35 prehashes=2 prefilterbits=2 buildpasses=1 \
          bits=16 hashes=3 \
          threads=64 interleaved=false \
          ecc=f tossbadreads=f fixspikes=t deterministic=t \
          hist="$hist_in" histout="$hist_out" peaks="$peaks" \
          zerobin=t pzc=t histlen=10000 \
          minq=6 minprob=0.5 \
          -Xmx240g -eoom -da overwrite=t
echo "Finished bbnorm.sh at $(date)"

# Check output files and verify depth cap
if [ ! -f "$norm_out_left" ] || [ ! -f "$norm_out_right" ]; then
    echo "Error: Normalization failed. Output files $norm_out_left or $norm_out_right not generated."
    exit 1
fi

# Check if normalization ran successfully
echo "Normalization completed successfully at $(date)"
echo "Kmer depth capped at 100. Check $hist_out to confirm depth distribution."

# Get sizes of normalized files
norm_left_size=$(stat -c %s "$norm_out_left")
norm_right_size=$(stat -c %s "$norm_out_right")
echo "Normalized left file size: $(numfmt --to=iec-i --suffix=B $norm_left_size)"
echo "Normalized right file size: $(numfmt --to=iec-i --suffix=B $norm_right_size)" 