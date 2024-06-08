nextflow.enable.dsl=2

params.tree = "../03_ebola/results_001/aln.ids_subsamp.phylip_phyml_tree.txt"
params.align = "../03_ebola/results_001/aln.ids_subsamp.phylip"
params.results="results"
params.collapse=0.1
params.collapseref=false
params.length = 18996
params.nboot=200
params.scale = 1.0 // tree branches scaling factor
params.minlen = -1 // set random exponential of mean meanlen to branches having length >= min-len
params.maxlen = -1 // set random exponential of mean meanlen to branches having length <= max-len
params.meanlen = 2.632133e-05 // set random exponential of mean meanlen
params.gtr=false // If true, simulation is done with gtr and all parameters inferred on true data

length = params.length
nboot=params.nboot
scale = params.scale
minlen=params.minlen
maxlen=params.maxlen
meanlen=params.meanlen
tree = file(params.tree)
align=file(params.align)
results=params.results
collapse=params.collapse
collapseref=params.collapseref
gtr=params.gtr

process preprocessTrueTree {
    publishDir "$results/", mode: 'link'

    label 'gotree'

    input:
    path tree
    val scale
    val minlen
    val maxlen
    val meanlen

    output:
    path "${tree.baseName}_scale.nw"

    script:
    """
    gotree brlen setrand -i $tree --min-len $minlen --max-len $maxlen --mean $meanlen -o tmp --seed 123456789
    gotree brlen scale -i tmp -f $scale | gotree reroot midpoint -o ${tree.baseName}_scale.nw
    """
}

process RootSequence {
    publishDir "$results/", mode: 'link'

    label 'goalign'

    input:
    path align
    val length
    
    output:
    path "root.fasta"

    script:
    """
    goalign consensus -i $align -p | goalign reformat fasta -p > consensus.fasta
    # We keep positions with N, ?, and -
    goalign stats char  --per-sites --only "N" -i consensus.fasta | grep "\t1" | cut -f 1 > n_pos
    goalign stats char  --per-sites --only "?" -i consensus.fasta | grep "\t1" | cut -f 1 >> n_pos
    goalign stats char  --per-sites --only "-" -i consensus.fasta | grep "\t1" | cut -f 1 >> n_pos
    # we generate random nucleotides to replace N, ?, and -
    goalign random --seed 123456789 -l \$(wc -l n_pos) -n 1| grep -v ">" |  fold -w1 > seq
    # we write the replacement file
    paste n_pos seq | awk '{print "consensus\t" \$1 "\t" \$2}' > replace.txt
    # we replace all the N,?,- positions with the random nucleotides
    # And we take a subsample of length \$length
    goalign replace -i consensus.fasta --posfile replace.txt | goalign sample sites -l $length --consecutive=false > root.fasta
    """
}

process simulateMSA {
    publishDir "$results/", mode: 'link'

    label 'snag'

    input:
    path tree
    path root
    val gtr

    output:
    path "msa.phylip"
    path "tree_realized.nw"
    path "rates.txt"

    script:
    snagparms="jc -gamma=false -gamma-cat=1"
    if(gtr){snagparms="gtr -parameters 0.09084,0.72307,0.07212,0.02603,0.97950,0.00979,0.31900,0.21423,0.19823,0.26854 -gamma -gamma-cat 4 -alpha 0.244"}
    """
    head -n 1 ${tree} > tree.nw
    snag -model $snagparms -intree ${tree} -root-seq $root -out-align msa.phylip -out-rates rates.txt -out-trees tree_realized.nw -seed 123456789
    """
}


process inferRefTree {
    publishDir "$results/", mode: 'link'

    label 'phyml'

    input:
    path msa

    output:
    path "${msa}_phyml_stats.txt"
    path "${msa}_phyml_tree.txt"

    script:
    """
    phyml -i ${msa} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
    """
}

process inferRefTreeRAxML {
    publishDir "$results/", mode: 'link'

    label 'raxml'

    input:
    path msa

    output:
    path "${msa}.raxml.bestTree"

    script:
    """
    raxml-ng --msa ${msa} --model GTR+G4 --threads 1 --seed 123456
    """
}

process genWeightBoot {
    publishDir "$results/", mode: 'link'

    label 'goalign'

    input:
    path msa
    val nboot

    output:
    path "weights.txt"

    script:
    """
    goalign build weightboot -p -i ${msa} -n ${nboot} > weights.txt
    """
}

process inferWeightBootTrees {
     label 'phyml'

     input:
     path weight
     path msa

     output:
     path "${msa}_phyml_stats.txt"
     path "${msa}_phyml_tree.txt"

     script:
     """
     phyml -i ${msa} --weights=${weight} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
     """
}

process computeWeightSupports {
     label 'gotree'
     
     publishDir "$results/", mode: 'link'
     
     input:
     path reftre
     path boot
     val length
     val collapse
     val collapseref

     output:
     path "reftree_weightsupport.nw"

     script:
     if( collapseref )
     """
     gotree collapse length -i ${boot} -l ${collapse/length} |  gotree compute support fbp -i <(gotree collapse length -i ${reftre} -l ${collapse/length}) -b - -o reftree_weightsupport.nw
     """
     else
     """
     gotree collapse length -i ${boot} -l ${collapse/length} |  gotree compute support fbp -i ${reftre} -b - -o reftree_weightsupport.nw
     """
}

process genSeqBoot {
    publishDir "$results/", mode: 'link'

    label 'goalign'

    input:
    path msa
    val nboot

    output:
    path "boot_*"

    script:
    """
    goalign build seqboot -p -n ${nboot} -i ${msa} -n ${nboot} -o boot_
    """
}

process inferSeqBootTrees {
     label 'phyml'

     input:
     path msa

     output:
     path "${msa}_phyml_stats.txt"
     path "${msa}_phyml_tree.txt"

     script:
     """
     phyml -i ${msa} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
     """
}

process inferSeqBootTreesRAxML {
     label 'raxml'

     input:
     path msa

     output:
     path "${msa}.raxml.bestTree"

     script:
     """
     raxml-ng --msa ${msa} --model GTR+G4 --threads 1 --seed 123456
     """
}

process computeSeqbootSupports {
     publishDir "$results/", mode: 'link'

     label 'gotree'

     input:
     path reftre
     path boot
     val length
     val collapse
     val collapseref

     output:
     path "reftree_bootsupport*.nw"

     script:
     if( collapseref )
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i <(gotree collapse length -i ${reftre} -l ${collapse/length}) -b - -o reftree_bootsupport.nw
     gotree compute support fbp -i <(gotree collapse length -i ${reftre} -l ${collapse/length}) -b ${boot} -o reftree_bootsupport_nocollapse.nw
     """
     else
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i ${reftre} -b - -o reftree_bootsupport.nw
     gotree compute support fbp -i ${reftre} -b ${boot} -o reftree_bootsupport_nocollapse.nw
     """
}

process computeSeqbootSupportsRAxML {
     publishDir "$results/", mode: 'link'

     label 'gotree'

     input:
     path reftre
     path boot
     val length
     val collapse
     val collapseref

     output:
     path "reftree_bootsupport_raxml*.nw"

     script:
     if( collapseref ) 
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i <(gotree collapse length -i ${reftre} -l ${collapse/length}) -b - -o reftree_bootsupport_raxml.nw
     gotree compute support fbp -i <(gotree collapse length -i ${reftre} -l ${collapse/length}) -b ${boot} -o reftree_bootsupport_raxml_nocollapse.nw
     """
     else
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i ${reftre} -b - -o reftree_bootsupport_raxml.nw
     gotree compute support fbp -i ${reftre} -b ${boot} -o reftree_bootsupport_raxml_nocollapse.nw
     """
}


workflow {
	 root = RootSequence(align, length)
	 pretree = preprocessTrueTree(tree,scale, minlen, maxlen, meanlen)
	 outsimu = simulateMSA(pretree, root, gtr)
	 simualn = outsimu[0]
	 simutree= outsimu[1]
	 simurate= outsimu[2]

	 outref = inferRefTree(simualn)
	 refstats = outref[0]
	 reftree  = outref[1]

	 reftreerax = inferRefTreeRAxML(simualn)

	 weights = genWeightBoot(simualn, nboot)

	 outbootweight = inferWeightBootTrees(weights.splitText(by: 1, file:true), simualn)
	 bootweightstats = outbootweight[0]
	 bootweighttrees = outbootweight[1]


	 weightsupport = computeWeightSupports(reftree, bootweighttrees.collectFile(name: "boot.nw"), length, collapse, collapseref)
	 
	 seqboot = genSeqBoot(simualn, nboot)
	 
	 outseqboot = inferSeqBootTrees(seqboot.flatten())
	 seqbootstats = outseqboot[0]
	 seqboottrees = outseqboot[1]
	 
	 seqboottreesrax = inferSeqBootTreesRAxML(seqboot.flatten())
	 
	 seqbootsupport = computeSeqbootSupports(reftree, seqboottrees.collectFile(name: "boot.nw"), length, collapse, collapseref)
	 
	 seqbootsupportrax = computeSeqbootSupportsRAxML(reftreerax, seqboottreesrax.collectFile(name: "boot.nw"), length, collapse, collapseref)
}
