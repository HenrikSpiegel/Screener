# Extract reads from a single chromosome (here bgc)
samtools view -h sorted.aln.bam NZ_CP053893.1.region002 > NZ.CP053893.1.region002.sam

# Convert to fasta (minimap2 input)
samtools fasta  NZ.CP053893.1.region002.sam > NZ_CP053893.1.region006.fa

# get reference genome:
cp /home/projects/dtu_00009/people/henspi/git/AntibioticaScreening/project/data/simulated_data/input_genomes/NZ_CP053893_1.fa .

#map to genome
minimap2 -t 4 -a NZ_CP053893_1.fa NZ_CP053893.1.region006.fa > NZ_CP053893_1.region002.aln.sam

# sort and coverage

samtools sort --threads 4 --write-index -o sorted.NZ_CP053893_1.aln.bam NZ_CP053893_1.region002.aln.sam

samtools mpileup -a sorted.NZ_CP053893_1.aln.bam | awk '{{print $1"\t"$2"\t"$4}}' > genome_coverage.tsv 