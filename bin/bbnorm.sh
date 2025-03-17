#SBATCH --partition=week-long-highmem
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --job-name=trinity_norm
#SBATCH --output=%j_bbnorm.out
#SBATCH --error=%j_bbnorm.err

source ~/.bashrc
conda activate transcriptome

in1=/nfs/home/jlamb/Projects/transcriptome/results/normalization/left.fa
in2=/nfs/home/jlamb/Projects/transcriptome/results/normalization/right.fa
out1=/nfs/home/jlamb/Projects/transcriptome/results/normalization/left.norm.fa
out2=/nfs/home/jlamb/Projects/transcriptome/results/normalization/right.norm.fa
target=100

bbnorm.sh in1="$in1" in2="$in2" out1="$out1" out2="$out2" \
          target="$target" mindepth=3 maxdepth=-1 passes=2 \
          k=25 prefilter=t prefiltersize=0.4 buildpasses=2 bits=32 \
          threads=48 interleaved=false \
          ecc=f tossbadreads=f \
          -Xmx240g