assembly_dir="data/simulated_data_300genomes/genera_tables"
full_dump="$assembly_dir/assembly_summary.txt"

if [ -f "$full_dump" ]; then
    echo "$full_dump exists"
else
    echo "Downloading dump from NCBI"
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt --directory-prefix "$assembly_dir"
fi

genus_extract=("Faecalibacterium" "Escherichia coli" "Bifidobacterium")

for genus in "${genus_extract[@]}"; do
    echo $genus
    head -n2 $full_dump | tail -n 1 | cut -c 3- > "$assembly_dir/$genus.tsv"
    grep "$genus" $full_dump >> "$assembly_dir/$genus.tsv"
done
