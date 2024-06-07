nextflow.enable.dsl=2

params.results="results"
params.align="data/5_cherries_scale_3.phy"
params.length=1000
params.nboot=200
params.shuffle=false
params.phymlopt=""

results=file(params.results)

process InferRefTreeRAxML {
    publishDir "$results/", mode: 'link'

    label 'raxml'

    input:
    file msa

    output:
    path "reftree_raxml.nw"

    script:
    """
    raxml-ng --msa ${msa} --model JC --threads 1 --seed 123456
    mv ${msa}.raxml.bestTree reftree_raxml.nw
    """
}

process InferRefTreeIQTree {
    publishDir "$results/", mode: 'link'

    label 'iqtree'

    input:
    file msa

    output:
    path "reftree_iqtree.nw"

    script:
    """
    iqtree -s ${msa} -m JC69+FO --seed 123456789
    mv ${msa}.treefile reftree_iqtree.nw
    """
}

process InferRefTreeIQTreeUFBoot {
    publishDir "$results/", mode: 'link'

    label 'iqtree'

    input:
    file msa

    output:
    path "reftree_iqtree_ufboot.nw"

    script:
    """
    iqtree -s ${msa} -m JC69+FO -B 1000 --seed 123456789
    mv ${msa}.treefile reftree_iqtree_ufboot.nw
    """
}

process InferRefTreePhyML {
    publishDir "$results/", mode: 'link'

    label 'phyml'

    input:
    path msa
    val phymlopt
    
    output:
    path "reftree_phyml.nw"
    
    script:
    """
    phyml -i ${msa} -m JC69 -c 1 -a e -f e -d nt -o tlr -b 0 ${phymlopt} --r_seed 123456
    mv ${msa}_phyml_tree.txt reftree_phyml.nw
    phyml -i ${msa} -m JC69 -c 1 -a e -f e -d nt -o tlr -b 1000 ${phymlopt} --r_seed 123456
    mv ${msa}_phyml_tree.txt reftree_phyml_boot.nw
    """
}

process SeqBoot {
    publishDir "$results/seqboot/", mode: 'link'

    label 'goalign'

    input:
    file msa
    val nboot
    val shuffle

    output:
    file "boot_*"

    script:
    shuf=(shuffle ? "-S": "")
    """
    goalign build seqboot -p -n $nboot $shuf -i ${msa} -o boot_ --seed 123456789
    """
}

process WeightBoot {
    publishDir "$results/", mode: 'link'

    label 'goalign'

    input:
    file msa
    val nboot

    output:
    file "weights.txt"

    script:
    """
    goalign build weightboot -p -i ${msa} -n $nboot --seed 123456789 > weights.txt
    """
}

process SeqBootTreesRAxML {
     label 'raxml'

     input:
     file msa

     output:
     file "${msa}_boottree_raxml.nw"

     script:
     """
     raxml-ng --msa ${msa} --model JC --threads 1 --seed 123456
     mv ${msa}.raxml.bestTree ${msa}_boottree_raxml.nw
     """
}

process SeqBootTreesPhyML {
    label 'phyml'

    input:
    path msa
    val phymlopt

    output:
    path "${msa}_boottree_phyml.nw"
    
    script:
    """
    phyml -i ${msa} -m JC69 -c 1 -a e -f e -d nt -o tlr -b 0 ${phymlopt} --r_seed 123456
    mv ${msa}_phyml_tree.txt ${msa}_boottree_phyml.nw
    """
}

process SeqBootTreesIQTree {
    label 'iqtree'

    input:
    path msa
    
    output:
    path "${msa}_boottree_iqtree.nw"
    
    script:
    """
    iqtree -s ${msa} -m JC69+FO --seed 123456789
    mv ${msa}.treefile ${msa}_boottree_iqtree.nw
    """
}

process WeightBootTreesPhyML {
    label 'phyml'

    input:
    path msa
    path weight
    val phymlopt

    output:
    path "${msa}_weighttree_phyml.nw"
    
    script:
    """
    phyml -i ${msa} --weights=${weight} -m JC69 -c 1 -a e -f e -d nt -o tlr -b 0 ${phymlopt} --r_seed 123456
    mv ${msa}_phyml_tree.txt ${msa}_weighttree_phyml.nw
    """
}

process SeqBootSupports {
     publishDir "$results/", mode: 'link'

     label 'gotree'

     input:
     file ref
     file boot
     val length

     output:
     file "reftree_raxml_bootsupport.nw"
     file "reftree_raxml_bootsupport_collapse.nw"

     script:
     """
     gotree compute support fbp -i ${ref} -b ${boot} -o reftree_raxml_bootsupport.nw
     gotree collapse length -i ${boot} -l ${0.1/length} | gotree compute support fbp -i ${ref} -b - -o reftree_raxml_bootsupport_collapse.nw
     """
}

process SeqBootSupportsPhyML {
     publishDir "$results/", mode: 'link'

     label 'gotree'

     input:
     file ref
     file boot
     val length

     output:
     file "reftree_phyml_bootsupport.nw"
     file "reftree_phyml_bootsupport_collapse.nw"

     script:
     """
     gotree compute support fbp -i ${ref} -b ${boot} -o reftree_phyml_bootsupport.nw
     gotree collapse length -i ${boot} -l ${0.1/length} | gotree compute support fbp -i ${ref} -b - -o reftree_phyml_bootsupport_collapse.nw
     """
}

process SeqBootSupportsIQTree {
     publishDir "$results/", mode: 'link'

     label 'gotree'

     input:
     file ref
     file boot
     val length

     output:
     file "reftree_iqtree_bootsupport.nw"
     file "reftree_iqtree_bootsupport_collapse.nw"

     script:
     """
     gotree compute support fbp -i ${ref} -b ${boot} -o reftree_iqtree_bootsupport.nw
     gotree collapse length -i ${boot} -l ${0.1/length} | gotree compute support fbp -i ${ref} -b - -o reftree_iqtree_bootsupport_collapse.nw
     """
}


process WeightBootSupportsPhyML {
     publishDir "$results/", mode: 'link'

     label 'gotree'

     input:
     file ref
     file boot
     val length

     output:
     file "reftree_phyml_weightsupport.nw"
     file "reftree_phyml_weightsupport_collapse.nw"

     script:
     """
     gotree compute support fbp -i ${ref} -b ${boot} -o reftree_phyml_weightsupport.nw
     gotree collapse length -i ${boot} -l ${0.1/length} | gotree compute support fbp -i ${ref} -b - -o reftree_phyml_weightsupport_collapse.nw
     """
}

workflow{
	msa = Channel.fromPath(params.align)
	length=params.length
	nboot=params.nboot
	shuffle=params.shuffle
	phymlopt=params.phymlopt

	bootal = SeqBoot(msa,nboot,shuffle)
	weights= WeightBoot(msa,nboot)
	refrax = InferRefTreeRAxML(msa)
	refphy = InferRefTreePhyML(msa,phymlopt)
	refiqtb= InferRefTreeIQTreeUFBoot(msa)
	refiqt  = InferRefTreeIQTree(msa)

	bootrax = SeqBootTreesRAxML(bootal.flatten())
	bootphy = SeqBootTreesPhyML(bootal.flatten(),phymlopt)
	bootiqt = SeqBootTreesIQTree(bootal.flatten())
	weightphy = WeightBootTreesPhyML(msa.first(), weights.splitText(by: 1, file:true), phymlopt)

	bootraxcoll = bootrax.collectFile(name: "boottrees_raxml.nw")
	bootphycoll = bootphy.collectFile(name: "boottrees_phyml.nw")
	bootiqtcoll = bootiqt.collectFile(name: "boottrees_iqtree.nw")
	weightphycoll = weightphy.collectFile(name: "weighttrees_phyml.nw")

	bootraxcoll.subscribe { it -> it.copyTo("$results/") }
	bootphycoll.subscribe { it -> it.copyTo("$results/") }
	bootiqtcoll.subscribe { it -> it.copyTo("$results/") }
	weightphycoll.subscribe { it -> it.copyTo("$results/") }

	SeqBootSupports(refrax,bootraxcoll,length)
	SeqBootSupportsPhyML(refphy,bootphycoll,length)
	SeqBootSupportsIQTree(refiqt,bootiqtcoll,length)
	WeightBootSupportsPhyML(refphy,weightphycoll,length)
}
