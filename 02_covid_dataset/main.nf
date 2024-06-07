params.msa = "data/alignment.fasta.xz"
params.results="results"
params.collapse=0.1
params.collapseref=false

results=params.results
nboot = 200
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
    xz -c -d ${msa} | goalign reformat phylip > ${msa.baseName}.phylip
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
    xz -c -d ${msa} | goalign stats length
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
    file msa6
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
    msareformat=reformat(msa)
    len=length(msa).map{it ->  Integer.parseInt(it.trim()) }
    reftreeall=inferRefTree(msareformat)
    reftreephy=reftreeall[1]
    reftreerax=inferRefTreeRAxML(msareformat)
    weights=genWeightBoot(msareformat,nboot).splitText(by: 1, file:true)
    weightbootsall=inferWeightBootTrees(weights,msareformat)
    weightboots=weightbootsallall[1].collectFile(name: "boot.nw")
    computeWeightSupports(reftreephy, weightboots, len, collapse, collapseref)
    seqbootaligns = genSeqBoot(msareformat, nboot).flatten()
    seqboottreesall=inferSeqBootTrees(seqbootaligns)
    seqboottrees=seqboottreesall[1].collectFile(name: "boot.nw")
    seqboottreesrax=inferSeqBootTreesRAxML(seqbootaligns).collectFile(name: "boot.nw")
    computeSeqbootSupports(reftreephy, seqboottrees, len, collapse, collapseref)
    computeSeqbootSupportsRAxML(reftreerax, seqboottreesrax, len, collapse, collapseref)
}
