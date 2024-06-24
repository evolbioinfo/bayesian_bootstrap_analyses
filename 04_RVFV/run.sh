singularity pull docker://evolbioinfo/goalign:dev0537492
singularity pull docker://evolbioinfo/gotree:devb324e73
singularity pull docker://evolbioinfo/r-extended:v4.3.3

goalign="singularity exec goalign_dev0537492.sif goalign"
gotree="singularity exec gotree_devb324e73.sif gotree"
rscript="singularity exec r-extended_v4.2.3_2 Rscript"

nextflow run convert.nf --ncbi data/sequences.fasta --gb data/sequence.gb

nextflow run main.nf --alignment "results/full_L_dedup_aligned_clean.fasta" --results "results_L" --collapse 0.1 --collapseref "false"
$gotree stats edges -i results_L/reftree_bootsupport.nw > results_L/edges_ref.txt
$gotree stats edges -i results_L/reftree_weightsupport.nw > results_L/edges_ref_weight.txt
$gotree stats edges -i results_L/reftree_bootsupport_nocollapse.nw > results_L/edges_ref_nocollapse.txt

nextflow run main.nf --alignment "results/full_M_dedup_aligned_clean.fasta" --results "results_M" --collapse 0.1 --collapseref "false"
$gotree stats edges -i results_M/reftree_bootsupport.nw > results_M/edges_ref.txt
$gotree stats edges -i results_M/reftree_weightsupport.nw > results_M/edges_ref_weight.txt
$gotree stats edges -i results_M/reftree_bootsupport_nocollapse.nw > results_M/edges_ref_nocollapse.txt

nextflow run main.nf --alignment "results/full_S_dedup_aligned_clean.fasta" --results "results_S" --collapse 0.1 --collapseref "false"
$gotree stats edges -i results_S/reftree_bootsupport.nw > results_S/edges_ref.txt
$gotree stats edges -i results_S/reftree_weightsupport.nw > results_S/edges_ref_weight.txt
$gotree stats edges -i results_S/reftree_bootsupport_nocollapse.nw > results_S/edges_ref_nocollapse.txt

rm homoplasies.txt
for d in L M S
do
    ALI=$($goalign stats char -p --per-sites -i results_${d}/full_${d}_dedup_aligned_clean.phylip | bin/pars.pl)
    PARS=$(grep Pars results_${d}/full_${d}_dedup_aligned_clean.phylip_phyml_stats.txt | awk '{print $3}')
    LEN=$($goalign stats length -p -i results_${d}/full_${d}_dedup_aligned_clean.phylip)
    ML=$($gotree stats -i results_${d}/full_${d}_dedup_aligned_clean.phylip_phyml_tree.txt | cut -f 6| tail -n 1)
    echo $(awk -v ml=$ML -v len=$LEN -v ali=$ALI 'BEGIN{print (ml*len)}')
    MLHOMO=$(awk -v ml=$ML -v len=$LEN -v ali=$ALI 'BEGIN{print (ml*len-ali)*100/ali}')
    echo -e "$d\t$ALI\t$PARS\t$MLHOMO" >> homoplasies.txt
done

$rscript script_rvfv.R

$goalign stats -i results_S/full_S_dedup_aligned_clean.phylip -p
$goalign stats -i results_M/full_M_dedup_aligned_clean.phylip -p
$goalign stats -i results_L/full_L_dedup_aligned_clean.phylip -p
$gotree stats -i results_S/full_S_dedup_aligned_clean.phylip_phyml_tree.txt
$gotree stats -i results_M/full_M_dedup_aligned_clean.phylip_phyml_tree.txt
$gotree stats -i results_L/full_L_dedup_aligned_clean.phylip_phyml_tree.txt
$goalign stats char --per-sites -i results_S/full_S_dedup_aligned_clean.phylip -p  | bin/pars.pl
$goalign stats char --per-sites -i results_M/full_M_dedup_aligned_clean.phylip -p  | bin/pars.pl
$goalign stats char --per-sites -i results_L/full_L_dedup_aligned_clean.phylip -p  | bin/pars.pl
