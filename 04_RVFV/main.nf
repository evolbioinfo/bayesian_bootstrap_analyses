params.alignment = "data/alignment.fas"
params.results="results"
params.collapse=0.1
params.collapseref=false

nboot = 200
results=params.results
alignment = file(params.alignment)
collapse=params.collapse
collapseref=params.collapseref

process reformat {
    label 'goalign'

    input:
    file msa

    output:
    file "${msa.baseName}.phylip"

    script:
    """
    goalign rename -i $msa -e "([^\\s]+).*" -b "\\\$1" | goalign reformat phylip | goalign replace -p -s "!" -n "-" | goalign replace -p -s "?" -n "-" | goalign clean sites -c 0.6 --char=GAP -p | goalign clean seqs -c 0.6 --char=GAP -p > ${msa.baseName}.phylip
    """
}

process subSample {
    label 'goalign'

    input:
    file msa
	
    output:
    file "${msa.baseName}_subsamp.phylip"

    script:
    """
    goalign sample seqs -n 800 -i ${msa} -p --seed 123456789 | goalign clean sites -p  --char=GAP -c 0.8 -o ${msa.baseName}_subsamp.phylip
    """
}

process alilen {
    label 'goalign'

    input:
    file msa

    output:
    stdout

    script:
    """
    printf \$(goalign stats length -p -i ${msa})
    """
}

process inferRefTree {
    publishDir "$results/", mode: 'link'

    label 'phyml'

    input:
    file msa

    output:
    file "${msa}_phyml_stats.txt"
    file "${msa}_phyml_tree.txt"

    script:
    """
    phyml -i ${msa} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
    """
}

process inferRefTreeRAxML {
    publishDir "$results/", mode: 'link'

    label 'raxml'

    input:
    file msa

    output:
    file "${msa}.raxml.bestTree"

    script:
    """
    raxml-ng --msa ${msa} --model GTR+G4 --threads 1 --seed 123456
    """
}

process genWeightBoot {
    label 'goalign'

    publishDir "$results/", mode: 'link'

    input:
    file msa
    val nboot

    output:
    file "weights.txt"

    script:
    """
    goalign build weightboot -p -i ${msa} -n ${nboot} > weights.txt
    """
}

process inferWeightBootTrees {
     label 'phyml'

     input:
     file w
     file msa

     output:
     file "${msa}_phyml_stats.txt"
     file "${msa}_phyml_tree.txt"

     script:
     """
     phyml -i ${msa} --weights=${w} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
     """
}

process computeWeightSupports {
    label 'gotree'
    
     publishDir "$results/", mode: 'link'
     
     input:
     file ref
     file boot
     val length
     val collapse
     val collapseref

     output:
     file "reftree_weightsupport.nw"

     script:
     if( collapseref )
     """
     gotree collapse length -i ${boot} -l ${collapse/length} |  gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b - -o reftree_weightsupport.nw
     """
     else
     """
     gotree collapse length -i ${boot} -l ${collapse/length} |  gotree compute support fbp -i ${ref} -b - -o reftree_weightsupport.nw
     """
}

process genSeqBoot {
    label 'goalign'
    
    publishDir "$results/", mode: 'link'

    input:
    file msa
    val nboot

    output:
    file "boot_*"

    script:
    """
    goalign build seqboot -p -n ${nboot} -i ${msa} -n ${nboot} -o boot_
    """
}

process inferSeqBootTrees {
     label 'phyml'

     input:
     file msa

     output:
     file "${msa}_phyml_stats.txt"
     file "${msa}_phyml_tree.txt"

     script:
     """
     phyml -i ${msa} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
     """
}

process inferSeqBootTreesRAxML {
     label 'raxml'

     input:
     file msa

     output:
     //file "${msa}_phyml_stats.txt"
     file "${msa}.raxml.bestTree"

     script:
     """
     raxml-ng --msa ${msa} --model GTR+G4 --threads 1 --seed 123456
     """
}

process computeSeqbootSupports {
     label 'gotree'
     
     publishDir "$results/", mode: 'link'
     
     input:
     file ref
     file boot
     val length
     val collapse
     val collapseref

     output:
     file "reftree_bootsupport*.nw"

     script:
     if( collapseref )
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b - -o reftree_bootsupport.nw
     gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b ${boot} -o reftree_bootsupport_nocollapse.nw
     """
     else
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i ${ref} -b - -o reftree_bootsupport.nw
     gotree compute support fbp -i ${ref} -b ${boot} -o reftree_bootsupport_nocollapse.nw
     """
}

process computeSeqbootSupportsRAxML {
     label 'gotree'
     
     publishDir "$results/", mode: 'link'
     
     input:
     file ref
     file boot
     val length
     val collapse
     val collapseref

     output:
     file "reftree_bootsupport_raxml.nw"

     script:
     if( collapseref ) 
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b - -o reftree_bootsupport_raxml.nw
     """
     else
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i ${ref} -b - -o reftree_bootsupport_raxml.nw
     """
}


workflow {
    alignment = file(params.alignment)
    phy = reformat(alignment)
    len = alilen(phy).map{it ->  Integer.parseInt(it.trim())}
    reftree = inferRefTree(phy)
    weights=genWeightBoot(phy,nboot).splitText(by: 1, file:true)
    weightboot = inferWeightBootTrees(weights,phy)
    weightboottrees = weightboot[1].collectFile(name: "boot.nw")
    weightsupport = computeWeightSupports(reftree[1], weightboottrees,len, collapse, collapseref)
    
    freqboot = genSeqBoot(phy, nboot).flatten()
    seqboot = inferSeqBootTrees(freqboot)
    seqboottrees = seqboot[1].collectFile(name: "boot.nw")
    computeSeqbootSupports(reftree[1], seqboottrees, len, collapse, collapseref)    
}
