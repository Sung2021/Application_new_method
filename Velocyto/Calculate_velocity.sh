# The scRNA-seq data was prepared using soloFeature, use the output of soloFeature as input for Velocyto.

# Convert the output of soloFeature (BAM file) to Velocyto compatible format (LOOM file).
velocyto prepare /path/to/soloFeature_output/output_prefix.Aligned.sortedByCoord.out.bam \
               /path/to/velocyto_output/velocyto_input.loom \
               /path/to/reference_genome/GTF_annotation.gtf

# Calculate RNA velocity using Velocyto.

velocyto run10x /path/to/velocyto_output/velocyto_input.loom /path/to/velocyto_output
