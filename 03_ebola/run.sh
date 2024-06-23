singularity pull docker://evolbioinfo/goalign:dev0537492
singularity pull docker://evolbioinfo/gotree:devb324e73
singularity pull docker://evolbioinfo/r-extended:v4.2.3_2

goalign="singularity exec goalign_dev0537492.sif goalign"
gotree="singularity exec gotree_devb324e73.sif gotree"
rscript="singularity exec r-extended_v4.2.3_2 Rscript"

git clone git@github.com:evolbioinfo/bdei.git
nextflow run main.nf --msa bdei/ebola/data/aln.ids.fa --collapse 0.1 --collapseref false --results results
$gotree stats edges -i results/reftree_bootsupport.nw > results/edges_ref.txt
$gotree stats edges -i results/reftree_bootsupport_nocollapse.nw > results/edges_ref_nocollapse.txt
$gotree stats edges -i results/reftree_weightsupport.nw > results/edges_ref_weight.txt
$gotree stats edges -i results/reftree_bootsupport_raxml.nw > results/edges_ref_rax.txt
$gotree compare edges -i results/reftree_bootsupport.nw -c results/reftree_weightsupport.nw > results/comptrees.txt
$gotree reformat newick --format nexus -i results/align.nx.con.tre | head -n 1 | $gotree compare edges -i results/reftree_bootsupport.nw -c - > results/bayessample.txt
$gotree reformat newick --format nexus -i results/align.nx.con.tre | head -n 1 | $gotree compare edges -i results/reftree_weightsupport.nw -c - > results/bayesweight.txt
$gotree compare edges -i results/reftree_bootsupport_raxml.nw -c results/reftree_bootsupport.nw > results/comptreesrax.txt
$gotree compare edges -i results/reftree_bootsupport_raxml_nocollapse.nw -c results/reftree_bootsupport_nocollapse.nw > results/comptreesrax_nocollapse.txt

rm homoplasies.txt
ALI=$($goalign stats char -p --per-sites -i results/aln.ids_subsamp.phylip | bin/pars.pl)
PARS=$(grep Pars results/aln.ids_subsamp.phylip_phyml_stats.txt| awk '{print $3}')
LEN=$($goalign stats length -p -i results/aln.ids_subsamp.phylip)
ML=$($gotree stats -i results/aln.ids_subsamp.phylip_phyml_tree.txt | cut -f 6| tail -n 1)
HOMO=$(awk -v pars=$PARS -v ali=$ALI 'BEGIN{print (pars-ali)*100/ali}')
MLHOMO=$(awk -v ml=$ML -v len=$LEN -v ali=$ALI 'BEGIN{print (ml*len-ali)*100/ali}')
echo -e "$ALI\t$PARS\t$HOMO\t$MLHOMO" >> homoplasies.txt

$rscript script_ebola.R

# Number of min parsimony steps
$goalign stats char --per-sites -i aln.ids_subsamp.phylip -p | bin/pars.pl
# Size of the inferred tree (column sumbrlen to multiply by 18996)
$gotree stats -i reftree_bootsupport.nw

