singularity pull docker://evolbioinfo/goalign:dev0537492
singularity pull docker://evolbioinfo/gotree:devb324e73
singularity pull docker://evolbioinfo/r-extended:v4.3.3

goalign="singularity exec goalign_dev0537492.sif goalign"
gotree="singularity exec gotree_devb324e73.sif gotree"
rscript="singularity exec r-extended_v4.3.3 Rscript"

## Random selection of 800 sequences among dataset of the paper doi:10.5802/crbiol.29
mkdir results_random_800
wget https://github.com/evolbioinfo/phylocovid/raw/CRAS/data/20200425/alignments/duplicates.txt
sed 's/,/\n/g' duplicates.txt > results_random_800/ids_all.txt
shuf -n 800 results_random_800/ids_all.txt > results_random_800/ids_800.txt
rm results_random_800/ids_all.txt duplicates.txt

# We extract & align sequences from GISAID data
nextflow run -c select.config select.nf \
	 --gisaidfasta "data/sequences_fasta_2022_03_09.tar.xz" \
	 --gisaidmeta "data/metadata_tsv_2022_03_09.tar.xz" \
	 --ids "results_random_800/ids_800.txt" \
	 --results "results_random_800"


$goalign rename -e '^.*\|(.*)$' -b '$1' -i results_random_800/nextalign.aligned.fasta  | xz -c - > results_random_800/nextalign.aligned_renamed.fasta.xz
# We clean sequences
$goalign mask -s 0 -l 55 -i results_random_800/nextalign.aligned_renamed.fasta.xz \
    | $goalign mask -s 29803 -l 101 \
    | $goalign mask --pos 186,1058,2093,3036,3129,6989,8021,10322,10740,11073,13407,14785,19683,20147,21136,24033,24377,25562,26143,26460,26680,28076,28825,28853,29699 \
    | xz -c - > results_random_800/nextalign.aligned_renamed_masked.fasta.xz

# We run the workflow
nextflow run main.nf --msa results_random_800/nextalign.aligned_renamed_masked.fasta.xz --collapse 0.1 --results results_random_800

# We compute branch metrics
$gotree stats edges -i results_random_800/reftree_bootsupport.nw > results_random_800/edges_ref.txt
$gotree stats edges -i results_random_800/reftree_bootsupport_nocollapse.nw > results_random_800/edges_ref_nocollapse.txt
$gotree stats edges -i results_random_800/reftree_weightsupport.nw > results_random_800/edges_ref_weight.txt
$gotree stats edges -i results_random_800/reftree_bootsupport_raxml.nw > results_random_800/edges_ref_rax.txt
$gotree compare edges -i results_random_800/reftree_bootsupport.nw -c results_random_800/reftree_weightsupport.nw > results_random_800/comptrees.txt
$gotree compare edges -i results_random_800/reftree_bootsupport_raxml.nw -c results_random_800/reftree_bootsupport.nw > results_random_800/comptreesrax.txt

# We compute homoplasy metrics
rm homoplasies.txt
ALI=$($goalign stats char --per-sites -i results_random_800/nextalign.aligned_renamed_masked.fasta.xz | bin/pars.pl)
PARS=$(grep Pars results_random_800/nextalign.aligned_renamed_masked.fasta.phylip_phyml_stats.txt | awk '{print $3}')
LEN=$($goalign stats length -i results_random_800/nextalign.aligned_renamed_masked.fasta.xz)
ML=$($gotree stats -i results_random_800/nextalign.aligned_renamed_masked.fasta.phylip_phyml_tree.txt | cut -f 6| tail -n 1)
echo $(awk -v ml=$ML -v len=$LEN -v ali=$ALI 'BEGIN{print (ml*len)}')
MLHOMO=$(awk -v ml=$ML -v len=$LEN -v ali=$ALI 'BEGIN{print (ml*len-ali)*100/ali}')
echo -e "$ALI\t$PARS\t$MLHOMO" >> homoplasies.txt

$rscript script_covid.R

# Size of the tree
$gotree stats -i reftree_bootsupport.nw
$goalign stats char --per-sites -i nextalign.aligned_renamed_masked.fasta.xz | bin/pars.pl

# Number of deduplicated sequences
$goalign replace -s 'N' -n '-' -i results_random_800/nextalign.aligned_renamed_masked.fasta.xz | $goalign replace -s '?' -n '-' | $goalign replace -s 'M' -n '-' | $goalign replace -s 'W' -n '-' | $goalign replace -s 'R' -n '-' | $goalign replace -s 'K' -n '-' | $goalign replace -s 'Y' -n '-' > results_random_800/tmp
# We compute a distance matrix by considering only mutations (gap to character = distance 0)
$goalign compute distance --gap-mut 0 -i results_random_800/tmp -p -t 10 -m pdist  > results_random_800/dist
# We keep connected components from this matrix
$rscript dedup_sequences.R results_random_800/dist results_random_800/names_tokeep.txt
