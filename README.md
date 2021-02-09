# nf-rnaseq


```sh
conda create -n rnaseq -y
conda activate rnaseq
conda install -c bioconda -c conda-forge bowtie2 fastp fastqc multiqc biopython==1.77 pandas samtools diamond=2.0.6 bbmap subread bioconda dna_features_viewer
```

# Usage

```sh
nextflow run rnaseq.nf --reads "data/*_R{1,2}.fq.gz" -profile hpc_slurm --mode "fastqc"
nextflow run rnaseq.nf --reads "data/*_R{1,2}.fq.gz" -profile hpc_slurm --mode "clean" -resume
nextflow run rnaseq.nf --reads "data/*_R{1,2}.fq.gz" -profile hpc_slurm --mode "all" -resume
```
