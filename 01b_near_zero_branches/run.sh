mkdir results_ebola_near_zero

OUDIR=results_ebola_near_zero
ALIGN="../03_ebola/results/aln.ids_subsamp.phylip"

singularity pull docker://evolbioinfo/goalign:dev0537492
singularity pull docker://evolbioinfo/gotree:devb324e73
singularity pull docker://evolbioinfo/r-extended:v4.3.3

goalign="singularity exec goalign_dev0537492.sif goalign"
gotree="singularity exec gotree_devb324e73.sif gotree"
rscript="singularity exec r-extended_v4.3.3 Rscript"

# We replace all ambigous characters by -
$goalign replace -s 'N' -n '-' -p -i $ALIGN | $goalign replace -s '?' -n '-' -p | $goalign replace -s 'M' -n '-' -p | $goalign replace -s 'W' -n '-' -p | $goalign replace -s 'R' -n '-' -p | $goalign replace -s 'K' -n '-' -p | $goalign replace -s 'Y' -n '-'  -p > $OUTDIR/tmp
# We compute a distance matrix by considering only mutations (gap to character = distance 0)
$goalign compute distance --gap-mut 0 -i $OUTDIR/tmp -p -t 10 -m pdist  > $OUTDIR/dist
# We keep connected components from this matrix
$rscript dedup_sequences.R $ALIGN $OUTDIR/names_tokeep.txt

# We keep the deduplicated sequences by names
$goalign subset -p -i $ALIGN -f $OUTDIR/names_tokeep.txt | $goalign clean sites -c 1 -p > $OUTDIR/aln.ids_subsamp_quasi_duplicates.phylip

# We launch bootstrap computations
nextflow run main.nf --align $OUTDIR/aln.ids_subsamp_quasi_duplicates.phylip --nboot 200 --results $OUTDIR/

# Stats PHYML NO COLLAPSE
$gotree stats edges -i $OUTDIR/reftree_phyml_bootsupport.nw > $OUTDIR/phyml_nocollapse_supports.txt
# Stats PHYML COLLAPSE
$gotree stats edges -i <(gotree collapse length -i $OUTDIR/reftree_phyml_bootsupport_collapse.nw -l 5.264820469621986e-06) > $OUTDIR/phyml_collapse_supports.txt
# Stats RAxML NO COLLAPSE
$gotree stats edges -i $OUTDIR/reftree_raxml_bootsupport.nw > $OUTDIR/raxml_nocollapse_supports.txt
# Stats RAxML COLLAPSE
$gotree stats edges -i <(gotree collapse length -i $OUTDIR/reftree_raxml_bootsupport_collapse.nw -l 5.264820469621986e-06) > $OUTDIR/raxml_collapse_supports.txt
# Stats IQTREE NO COLLAPSE
$gotree stats edges -i $OUTDIR/reftree_original_iqtree.nw > $OUTDIR/iqtree_nocollapse_supports.txt
# Stats IQTREE COLLAPSE
$gotree stats edges -i $OUTDIR/reftree_collapsed_iqtree.nw  > $OUTDIR/iqtree_collapse_supports.txt
# Compare IQTREE & PHYML NO COLLAPSE
$gotree compare edges -i $OUTDIR/reftree_phyml_bootsupport.nw -c $OUTDIR/reftree_original_iqtree.nw  > $OUTDIR/phyml_iqtree_nocollapse_supports.txt
# Compare IQTREE & PHYML COLLAPSE
$gotree compare edges -i <($gotree collapse length -i $OUTDIR/reftree_phyml_bootsupport_collapse.nw -l 5.264820469621986e-06) -c $OUTDIR/reftree_collapsed_iqtree.nw  > $OUTDIR/phyml_iqtree_collapse_supports.txt
# Compare RAxML & PHYML NO COLLAPSE
$gotree compare edges -i $OUTDIR/reftree_phyml_bootsupport.nw -c $OUTDIR/reftree_raxml_bootsupport.nw  > $OUTDIR/phyml_raxml_nocollapse_supports.txt
# Compare RAxML & PHYML COLLAPSE
$gotree compare edges -i <($gotree collapse length -i $OUTDIR/reftree_phyml_bootsupport_collapse.nw -l 5.264820469621986e-06) -c <($gotree collapse length -i $OUTDIR/reftree_raxml_bootsupport_collapse.nw -l 5.264820469621986e-06) > $OUTDIR/phyml_raxml_collapse_supports.txt

# Compare trees
$gotree compare trees -i $OUTDIR/reftree_phyml_bootsupport.nw -c $OURDIR/reftree_original_iqtree.nw
$gotree compare trees -i  <($gotree collapse length -i $OUTDIR/reftree_phyml_bootsupport_collapse.nw -l 5.264820469621986e-06) -c $OUTDIR/reftree_collapsed_iqtree.nw
$gotree compare trees -i $OUTDIR/reftree_phyml_bootsupport.nw -c $OUTDIR/reftree_raxml_bootsupport.nw
$gotree compare trees -i  <($gotree collapse length -i $OUTDIR/reftree_phyml_bootsupport_collapse.nw -l 5.264820469621986e-06) -c <($gotree collapse length -i $OUTDIR/reftree_raxml_bootsupport_collapse.nw -l 5.264820469621986e-06)

# R statistics & plots
$rscript stats.R $OUTDIR
