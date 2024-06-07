#!/bin/bash

INDIR=$1
ITOLKEY=$2

singularity pull docker://evolbioinfo/gotree:devb324e73

gotree="singularity exec gotree_devb324e73.sif gotree"
itol="$gotree upload itol --user-id $ITOLKEY -i"

## True tree
URL=$($gotree brlen scale -i data/5_cherries.nw  -f 1000 | $gotree brlen round -p 2 | $itol -)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/truetree.svg

# Simulated tree
URL=$(itol $INDIR/5_cherries_scale_3_simu.nw)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///g')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/5_cherries_scale_3_simu.svg

### Collapsed
## PhyML Tree
URL=$($gotree brlen scale -i $INDIR/reftree_phyml_bootsupport_collapse.nw -f 1000 | $gotree brlen round -p 2 | $itol -)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///g')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/reftree_phyml_bootsupport_collapse.svg
$gotree download itol -i $ID -c itol_config.txt --format png -o $INDIR/reftree_phyml_bootsupport_collapse.png

URL=$($gotree brlen scale -i $INDIR/reftree_phyml_weightsupport.nw -f 1000 | $gotree brlen round -p 2 | $itol -)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///g')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/reftree_phyml_weightsupport.svg
$gotree download itol -i $ID -c itol_config.txt --format png -o $INDIR/reftree_phyml_weightsupport.png

URL=$($gotree brlen scale -i $INDIR/reftree_phyml_weightsupport_collapse.nw -f 1000 | $gotree brlen round -p 2 | $itol -)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///g')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/reftree_phyml_weightsupport_collapse.svg
$gotree download itol -i $ID -c itol_config.txt --format png -o $INDIR/reftree_phyml_weightsupport_collapse.png

### RAxML Tree
URL=$($gotree brlen scale -i $INDIR/reftree_raxml_bootsupport_collapse.nw -f 1000 | $gotree brlen round -p 2 | $itol -)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///g')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/reftree_raxml_bootsupport_collapse.svg
$gotree download itol -i $ID -c itol_config.txt --format png -o $INDIR/reftree_raxml_bootsupport_collapse.png

### IQTree
URL=$($gotree brlen scale -i $INDIR/reftree_iqtree_bootsupport_collapse.nw -f 1000 | $gotree brlen round -p 2 | $itol -)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///g')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/reftree_iqtree_bootsupport_collapse.svg
$gotree download itol -i $ID -c itol_config.txt --format png -o $INDIR/reftree_iqtree_bootsupport_collapse.png

#### NOT Collapsed
### PhyML Tree
URL=$($gotree brlen scale -i $INDIR/reftree_phyml_bootsupport.nw -f 1000 | $gotree brlen round -p 2 | $itol -)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///g')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/reftree_phyml_bootsupport.svg
$gotree download itol -i $ID -c itol_config.txt --format png -o $INDIR/reftree_phyml_bootsupport.png

### RAxML Tree
URL=$($gotree brlen scale -i $INDIR/reftree_raxml_bootsupport.nw -f 1000 | $gotree brlen round -p 2 | $itol -)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///g')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/reftree_raxml_bootsupport.svg
$gotree download itol -i $ID -c itol_config.txt --format png -o $INDIR/reftree_raxml_bootsupport.png

### IQTree
URL=$($gotree brlen scale -i $INDIR/reftree_iqtree_bootsupport.nw -f 1000 | $gotree brlen round -p 2 | $itol -)
ID=$(echo $URL|sed 's/https:\/\/itol.embl.de\/tree\///g')
$gotree download itol -i $ID -c itol_config.txt --format svg -o $INDIR/reftree_iqtree_bootsupport.svg
$gotree download itol -i $ID -c itol_config.txt --format png -o $INDIR/reftree_iqtree_bootsupport.png
