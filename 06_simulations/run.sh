
singularity pull docker://evolbioinfo/goalign:dev0537492
singularity pull docker://evolbioinfo/gotree:devb324e73

goalign="singularity exec goalign_dev0537492.sif goalign"
gotree="singularity exec gotree_devb324e73.sif gotree"

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
