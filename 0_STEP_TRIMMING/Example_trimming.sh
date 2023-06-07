#!/bin/bash
#SBATCH --job-name="Trimming_test"
#SBATCH -o TRIM_Test.out
#SBATCH -e TRIM_Test.err
#SBATCH --cpus-per-task=8
#SBATCH --time=5:00:00
#SBATCH --mem=2GB
#SBATCH --exclude=gn[05-09]
#SBATCH --partition=FAST
#SBATCH --mail-user=irondon@cicbiogune.es
#SBATCH --mail-type=FAIL

source /share/apps/anaconda3/cic-env
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n cutadaptenv cutadapt
conda activate cutadaptenv


input1="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs/AC123_1.fastq.gz"
input2="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs/AC123_2.fastq.gz"

output1="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/AC123_1_out.fastq.gz"
output2="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/AC123_2_out.fastq.gz"


cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $output1 -p $output2 $input1 $input2

output11="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/AC123_1.fastq.gz"
output22="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/AC123_2.fastq.gz"

cutadapt -q 10 -m 20 -u -3 -U 3 -o $output11 -p $output22 $output1 $output2
              
