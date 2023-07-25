# STARsolo have already installed and the reference genome was prepared.
# Process the scRNA-seq data using STARsolo (soloFeature) to extract the features.

# STARsolo command:
STAR --genomeDir /path/to/reference_genome/ \
     --readFilesIn /path/to/input_R1.fastq,/path/to/input_R2.fastq \
     --soloType soloFeature \
     --soloBarcodeReadLength 12 \
     --outFileNamePrefix /path/to/soloFeature_output/output_prefix \
     --outSAMtype BAM Unsorted
