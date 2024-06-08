params.msa = "data/alignment_hmmer_subset.fasta.xz"
params.results="results"
params.collapse=0.1
params.collapseref=false

nboot = 200
results=params.results
msa = file(params.msa)
collapse=params.collapse
collapseref=params.collapseref


process alphabet {
    label 'goalign'

    tag "$id"
    publishDir "${results}${id}/", mode: 'link'

    input:
    tuple val(id), file(msa)

    output:
    tuple val(id), stdout
    tuple val(id), path("alphabet.txt")

    script:
    """
    goalign stats alphabet -p -i ${msa} > alphabet.txt
    cat alphabet.txt
    """
}

process rename {
    label 'goalign'

    tag "$id"
    publishDir "${results}${id}/", mode: 'link'

    input:
    tuple val(id), file(msa)

    output:
    tuple val(id), path("${msa.baseName}_renamed.phy")

    script:
    """
    goalign trim name -p -a -i $msa -o ${msa.baseName}_renamed.phy
    """
}

process alilen {
    label 'goalign'

    tag "$id"
    
    input:
    tuple val(id), file(msa)

    output:
    tuple val(id), stdout

    script:
    """
    printf \$(goalign stats length -p -i ${msa})
    """
}

process inferRefTree {
    tag "$id"
    publishDir "${results}${id}/", mode: 'link'

    label 'phyml'

    input:
    tuple val(id), val(alphabet), file(msa)

    output:
    tuple val(id), file("${msa}_phyml_stats.txt")
    tuple val(id), file("${msa}_phyml_tree.txt")

    script:
    if( alphabet == 'nucleotide' )
    """
    phyml -i ${msa} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
    """
    else
    """
    phyml -i ${msa} -m LG -c 4 -d aa -a e -o tlr -b 0 --r_seed 123456
    """
}

process genWeightBoot {
    label 'goalign'

    tag "$id"
    
    publishDir "${results}${id}/", mode: 'link'

    input:
    tuple val(id), file(msa)
    val nboot

    output:
    tuple val(id), file("weights.txt")

    script:
    """
    goalign build weightboot -p -i ${msa} -n ${nboot} > weights.txt
    """
}

process inferWeightBootTrees {
     tag "$id"
     label 'phyml'

     input:
     tuple val(id), val(alphabet), file(w), file(msa)

     output:
     tuple val(id), file("${msa}_phyml_stats.txt")
     tuple val(id), file("${msa}_phyml_tree.txt")

     script:
     if( alphabet == 'nucleotide' )
     """
     phyml -i ${msa} --weights=${w} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
     """
     else
     """
     phyml -i ${msa} --weights=${w} -m LG -c 4 -d aa -a e -o tlr -b 0 --r_seed 123456
     """
}

process computeWeightSupports {
     label 'gotree'

     tag "$id"

     publishDir "${results}${id}/", mode: 'link'
     
     input:
     tuple val(id), val(length), file(ref), file('boot')
     val collapse
     val collapseref

     output:
     tuple val(id), file("reftree_weightsupport.nw")

     script:
     if( collapseref )
     """
     gotree collapse length -i <(cat $boot) -l ${collapse/length} |  gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b - -o reftree_weightsupport.nw
     """
     else
     """
     gotree collapse length -i <(cat $boot) -l ${collapse/length} |  gotree compute support fbp -i ${ref} -b - -o reftree_weightsupport.nw
     """
}

process genSeqBoot {
    label 'goalign'

    tag "$id"

    publishDir "${results}${id}/", mode: 'link'

    input:
    tuple val(id), file(msa)
    val nboot

    output:
    tuple val(id), file("boot_*")

    script:
    """
    goalign build seqboot -p -n ${nboot} -i ${msa} -n ${nboot} -o boot_
    """
}

process inferSeqBootTrees {
     tag "$id"
     
     label 'phyml'

     input:
     tuple val(id), val(alphabet), file(msa)

     output:
     tuple val(id), file("${msa}_phyml_stats.txt")
     tuple val(id), file("${msa}_phyml_tree.txt")

     script:
     if( alphabet == 'nucleotide' )
     """
     phyml -i ${msa} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
     """
     else
     """
     phyml -i ${msa} -m LG -c 4 -d aa -a e -o tlr -b 0 --r_seed 123456
     """
}

process computeSeqbootSupports {
     label 'gotree'

     tag "$id"
     publishDir "${results}${id}/", mode: 'link'
     
     input:
     tuple val(id), val(length), file(ref), file('boot')
     val collapse
     val collapseref

     output:
     tuple val(id), file("reftree_bootsupport*.nw")

     script:
     if( collapseref )
     """
     gotree collapse length -i <(cat $boot) -l ${collapse/length} | gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b - -o reftree_bootsupport.nw
     gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b <(cat $boot) -o reftree_bootsupport_nocollapse.nw
     """
     else
     """
     gotree collapse length -i <(cat $boot) -l ${collapse/length} | gotree compute support fbp -i ${ref} -b - -o reftree_bootsupport.nw
     gotree compute support fbp -i ${ref} -b <(cat $boot) -o reftree_bootsupport_nocollapse.nw
     """
}

workflow {
    msa = Channel.fromPath(params.msa).splitCsv(header:false, sep:'\t').map{row -> [row[0], file(row[1])]}
    renamed=rename(msa)
    len = alilen(renamed).map{it ->  [it[0],Integer.parseInt(it[1].trim())]}
    alpharaw = alphabet(renamed)
    alphastr=alpharaw[0].map{it -> [it[0],it[1].trim()]}

    reftree = inferRefTree(alphastr.combine(renamed, by:0))

    weights=genWeightBoot(renamed,nboot).splitText(by: 1, file:true)
    weightboot = inferWeightBootTrees(alphastr.combine(weights, by:0).combine(renamed, by:0))
    weightboottrees = weightboot[1].groupTuple()
    weightsupport = computeWeightSupports(len.combine(reftree[1], by:0).combine(weightboottrees, by:0), collapse, collapseref)

    freqboot = genSeqBoot(renamed, nboot).transpose()
    seqboot = inferSeqBootTrees(alphastr.combine(freqboot, by:0))
    seqboottrees = seqboot[1].groupTuple()
    computeSeqbootSupports(len.combine(reftree[1], by:0).combine(seqboottrees, by:0), collapse, collapseref)
}
