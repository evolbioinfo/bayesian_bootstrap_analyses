params.msa = "data/alignment.fasta"
params.results="results"
params.collapse=0.1
params.collapseref=false

nboot = 200
results=params.results
msa = file(params.msa)
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
    goalign reformat phylip -i ${msa} > ${msa.baseName}.phylip
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
    goalign sample seqs -n 800 -i ${msa} -p --seed 123456789 -o ${msa.baseName}_subsamp.phylip
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
    goalign stats length -i ${msa}
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
     gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b $boot -o reftree_bootsupport_nocollapse.nw
     """
     else
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i ${ref} -b - -o reftree_bootsupport.nw
     gotree compute support fbp -i ${ref} -b $boot -o reftree_bootsupport_nocollapse.nw
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
     file "reftree_bootsupport_raxml*.nw"

     script:
     if( collapseref ) 
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b - -o reftree_bootsupport_raxml.nw
     gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b $boot -o reftree_bootsupport_raxml_nocollapse.nw
     """
     else
     """
     gotree collapse length -i ${boot} -l ${collapse/length} | gotree compute support fbp -i ${ref} -b - -o reftree_bootsupport_raxml.nw
     gotree compute support fbp -i ${ref} -b $boot -o reftree_bootsupport_raxml_nocollapse.nw
     """
}


workflow {
    msareformat=reformat(msa)
    msasub=subSample(msareformat)
    length=alilen(msa).map{it ->  Integer.parseInt(it.trim()) }
    phy=inferRefTree(msasub)
    refphy=phy[1]
    refrax=inferRefTreeRAxML(msasub)
    weights=genWeightBoot(msasub,nboot).splitText(by: 1, file:true)
    weighttrees=inferWeightBootTrees(weights, msasub)
    weightphy=weighttress[1].collectFile(name: "boot.nw")
    suppweight=computeWeightSupports(refphy, weightphy, collapse, collapseref)
    bootaligns=genSeqBoot(msasub, nboot).flatten()
    boottress=inferSeqBootTrees(bootaligns)
    bootphy=boottrees[1].collectFile(name:"boot.nw")
    bootrax=inferSeqBootTreesRAxML(bootaligns).collectFile(name: "boot.nw")

    computeSeqbootSupports(refphy, bootphy, length, collapse, collapseref)
    computeSeqbootSupportsRAxML(refrax, bootrax, length, collapse, collapseref)
}
