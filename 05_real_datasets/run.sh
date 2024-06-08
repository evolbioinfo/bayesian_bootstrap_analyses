
singularity pull docker://evolbioinfo/goalign:dev0537492
singularity pull docker://evolbioinfo/gotree:devb324e73

goalign="singularity exec goalign_dev0537492.sif goalign"
gotree="singularity exec gotree_devb324e73.sif gotree"

wget https://cme.h-its.org/exelixis/material/raxml_adaptive_data.tar.gz
tar -xzvf raxml_adaptive_data.tar.gz


# We extract dataset informations
echo "Sample\tDifficulty\tNSeq\tLength\tAlphabet" > datasets_infos.txt
for d in $(tail -n+2 raxml-adaptive-datasets/final_results_empirical.csv | cut -f 1,2 -d ',')
do
    SAMPLE=$(echo $d | cut -f 1 -d ',')
    DIFFICULTY=$(echo $d | cut -f 2 -d ',')
    NSEQ=$($goalign stats nseq -p -i raxml-adaptive-datasets/empirical/$SAMPLE/alignment.phy)
    LEN=$($goalign stats length -p -i raxml-adaptive-datasets/empirical/$SAMPLE/alignment.phy)
    ALPHABET=$($goalign stats alphabet -p -i raxml-adaptive-datasets/empirical/$SAMPLE/alignment.phy)

    echo -e "$SAMPLE\t$DIFFICULTY\t$NSEQ\t$LEN\t$ALPHABET" >> datasets_infos.txt
done

## 300 largest nucleotidic datasets
########################################################
sort -k 3 -n datasets_infos.txt | grep nucleotide | tail -n 300 > largest_nucleotide_datasets_300
# List samples not already analyzed
rm -f analyzed_nuc.txt
rm -f non_analyzed_nuc.txt
for d in $(cat largest_nucleotide_datasets_300 | cut -f 1)
do
    if [[ -f "results_largest_nuc_${d}/alignment_renamed.phy_phyml_stats.txt" ]]
    then
	echo -e "$d\traxml-adaptive-datasets/empirical/$d/alignment.phy" >> analyzed_nuc.txt
    else
	echo -e "$d\traxml-adaptive-datasets/empirical/$d/alignment.phy" >> non_analyzed_nuc.txt
    fi
done

# We run the workflow
nextflow run main.nf --msa non_analyzed_nuc.txt --collapse 0.1 --collapseref false --results results_largest_nuc_ -resume

# We compute homoplasy metrics
rm homoplasies.txt
for d in $(ls | grep results_largest_nuc)
do
    if [[ -f "${d}/alignment_renamed.phy_phyml_stats.txt" ]]
    then 
        ALI=$($goalign stats char -p --per-sites -i ${d}/alignment_renamed.phy | bin/pars.pl)
        PARS=$(grep Pars ${d}/alignment_renamed.phy_phyml_stats.txt| awk '{print $3}')
        LEN=$($goalign stats length -p -i ${d}/alignment_renamed.phy)
        ML=$($gotree stats -i ${d}/alignment_renamed.phy_phyml_tree.txt | cut -f 6| tail -n 1)
        HOMO=$(awk -v pars=$PARS -v ali=$ALI 'BEGIN{print (pars-ali)*100/ali}')
        MLHOMO=$(awk -v ml=$ML -v len=$LEN -v ali=$ALI 'BEGIN{print (ml*len-ali)*100/ali}')
        echo -e "$d\t$ALI\t$PARS\t$HOMO\t$MLHOMO" >> homoplasies.txt
        echo $HOMO > ${d}/homoplasy.txt
        echo $MLHOMO > ${d}/mlhomoplasy.txt
	$gotree stats edges -i ${d}/reftree_bootsupport.nw > ${d}/edges_ref.txt
	$gotree stats edges -i ${d}/reftree_weightsupport.nw > ${d}/edges_ref_weight.txt
	$gotree stats edges -i ${d}/reftree_bootsupport_nocollapse.nw > ${d}/edges_ref_nocollapse.txt
	$goalign stats length -p -i ${d}/boot_0.ph > ${d}/alilen.txt
	$goalign stats nseq -p -i ${d}/boot_0.ph > ${d}/ntaxa.txt
    fi
done


## 100 largest proteic datasets
########################################################
sort -k 3 -n datasets_infos.txt | grep protein| tail -n 100 > largest_proteic_datasets_100
# List samples not already analyzed
rm -f analyzed_prot.txt
rm -f non_analyzed_prot.txt
for d in $(cat largest_proteic_datasets_100 | cut -f 1)
do
    if [[ -f "results_largest_prot_${d}/alignment_renamed.phy_phyml_stats.txt" ]]
    then
	echo -e "$d\traxml-adaptive-datasets/empirical/$d/alignment.phy" >> analyzed_prot.txt
    else
	echo -e "$d\traxml-adaptive-datasets/empirical/$d/alignment.phy" >> non_analyzed_prot.txt
    fi
done
# Run the workflow
nextflow run main.nf --msa non_analyzed_prot.txt --collapse 0.1 --collapseref false --results results_largest_prot_ -resume

# Compute homoplasy metrics
rm homoplasies_prot.txt
for d in $(ls | grep results_largest_prot)
do
    if [[ -f "${d}/alignment_renamed.phy_phyml_stats.txt" ]]
    then 
        ALI=$($goalign stats char -p --per-sites -i ${d}/alignment_renamed.phy | bin/pars_prot.pl)
        PARS=$(grep Pars ${d}/alignment_renamed.phy_phyml_stats.txt| awk '{print $3}')
        LEN=$($goalign stats length -p -i ${d}/alignment_renamed.phy)
        ML=$($gotree stats -i ${d}/alignment_renamed.phy_phyml_tree.txt | cut -f 6| tail -n 1)
        HOMO=$(awk -v pars=$PARS -v ali=$ALI 'BEGIN{print (pars-ali)*100/ali}')
        MLHOMO=$(awk -v ml=$ML -v len=$LEN -v ali=$ALI 'BEGIN{print (ml*len-ali)*100/ali}')
        echo -e "$d\t$ALI\t$PARS\t$HOMO\t$MLHOMO" >> homoplasies_prot.txt
        echo $HOMO > ${d}/homoplasy.txt
        echo $MLHOMO > ${d}/mlhomoplasy.txt
	$gotree stats edges -i ${d}/reftree_bootsupport.nw > ${d}/edges_ref.txt
	$gotree stats edges -i ${d}/reftree_weightsupport.nw > ${d}/edges_ref_weight.txt
	$gotree stats edges -i ${d}/reftree_bootsupport_nocollapse.nw > ${d}/edges_ref_nocollapse.txt
	$goalign stats length -p -i ${d}/boot_0.ph > ${d}/alilen.txt
	$goalign stats nseq -p -i ${d}/boot_0.ph > ${d}/ntaxa.txt
    fi
done

