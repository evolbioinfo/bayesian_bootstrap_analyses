
singularity pull docker://evolbioinfo/goalign:dev0537492
singularity pull docker://evolbioinfo/gotree:devb324e73
singularity pull docker://evolbioinfo/r-extended:v4.3.3


goalign="singularity exec goalign_dev0537492.sif goalign"
gotree="singularity exec gotree_devb324e73.sif gotree"
rscript="singularity exec r-extended_v4.2.3_2 Rscript"

$TREE = "../03_ebola/results/aln.ids_subsamp.phylip_phyml_tree.txt"
$ALIGN = "../03_ebola/results/aln.ids_subsamp.phylip"

for scale in 1 4 16 64
do
    LEN=$((18996/$scale))
    RES=results_${scale}_$LEN
    nextflow run main.nf --tree $TREE --align $ALIGN --gtr true --results $RES --minlen 0 --maxlen 2.632133e-05 --meanlen 2.632133e-05 --scale $scale --nboot 200 --length $LEN -resume
    $gotree compare edges -i $RES/msa.phylip_phyml_tree.txt -c $RES/tree_realized.nw > $RES/edges.txt
    $gotree compare edges -i $RES/reftree_weightsupport.nw  -c $RES/tree_realized.nw > $RES/edges_ref_weight.txt
    $gotree compare edges -i $RES/reftree_bootsupport.nw  -c $RES/tree_realized.nw > $RES/edges_ref.txt
    $gotree compare edges -i $RES/reftree_bootsupport_nocollapse.nw  -c $RES/tree_realized.nw > $RES/edges_ref_nocollapse.txt
    $gotree compare edges -i $RES/msa.phylip_phyml_tree.txt -c <($gotree collapse length -l 0 -i $RES/tree_realized.nw) > $RES/edges_refcol.txt
    $gotree compare edges -i $RES/reftree_weightsupport.nw  -c <($gotree collapse length -l 0 -i $RES/tree_realized.nw) > $RES/edges_ref_weight_refcol.txt
    $gotree compare edges -i $RES/reftree_bootsupport.nw    -c <($gotree collapse length -l 0 -i $RES/tree_realized.nw) > $RES/edges_ref_refcol.txt
    $gotree compare edges -i $RES/reftree_bootsupport_nocollapse.nw  -c <($gotree collapse length -l 0 -i $RES/tree_realized.nw) > $RES/edges_ref_nocollapse_refcol.txt
done

# Compute homoplasy metrics
rm homoplasies.txt
for d in $(ls | grep results_)
do
    if [[ -f "${d}/msa.phylip_phyml_stats.txt" ]]
    then 
        ALI=$($goalign stats char -p --per-sites -i ${d}/msa.phylip | bin/pars.pl)
        PARS=$(grep Pars ${d}/msa.phylip_phyml_stats.txt| awk '{print $3}')
        LEN=$($goalign stats length -p -i ${d}/msa.phylip)
        ML=$($gotree stats -i ${d}/msa.phylip_phyml_tree.txt | cut -f 6| tail -n 1)
	echo $d
	echo $(awk -v ml=$ML -v len=$LEN -v ali=$ALI 'BEGIN{print (ml*len)}')	
        HOMO=$(awk -v pars=$PARS -v ali=$ALI 'BEGIN{print (pars-ali)*100/ali}')
        MLHOMO=$(awk -v ml=$ML -v len=$LEN -v ali=$ALI 'BEGIN{print (ml*len-ali)*100/ali}')
        echo -e "$d\t$ALI\t$PARS\t$HOMO\t$MLHOMO" >> homoplasies.txt
    fi
done

$rscript script_simulations.R
$rscript script_simulations_nocollapse.R

# Number of branches in common between CorrectedEbola and inferred tree
$gotree compare trees -i results_1_18996/aln.ids_subsamp.phylip_phyml_tree_scale.nw -c <(gotree collapse length -i results_1_18996/reftree_bootsupport.nw -l 0.000005264266161)
# Number of minimal parsimony steps
$goalign stats char --per-sites -i results_1_18996/msa.phylip -p | bin/pars.pl
# Size of the inferred tree (column sumbrlen to multiply by 18996)
$gotree stats -i results_1_18996/reftree_bootsupport.nw
# Size of the X4 inferred tree (column sumbrlen to multiply by 4749)
$gotree stats -i results_4_4749/reftree_bootsupport.nw
# Number of minimal parsimony steps
$goalign stats char --per-sites -i results_4_4749/msa.phylip -p | bin/pars.pl
# Size of the X16 inferred tree (column sumbrlen to multiply by 1187)
$gotree stats -i results_16_1187/reftree_bootsupport.nw
# Number of minimal parsimony steps
$goalign stats char --per-sites -i results_16_1187/msa.phylip -p | bin/pars.pl
# Size of the X64 inferred tree (column sumbrlen to multiply by 296)
$gotree stats -i results_64_296/reftree_bootsupport.nw
# Number of minimal parsimony steps
$goalign stats char --per-sites -i results_64_296/msa.phylip -p | bin/pars.p
