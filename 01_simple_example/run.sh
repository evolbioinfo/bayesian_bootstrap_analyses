ITOLKEY="XXXXXXXXXXX"

nextflow run main.nf --length 1000 --nboot 1000 --align data/Test4X3.txt --shuffle false --results results -resume
bin/drawtrees.sh results $ITOLKEY
