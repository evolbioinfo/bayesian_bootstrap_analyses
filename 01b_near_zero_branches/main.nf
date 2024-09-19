nextflow.enable.dsl=2

params.results="results"
params.align="../09_ebola/results_001/aln.ids_subsamp.phylip"
params.nboot=200

results=file(params.results)

process Unique {
    publishDir "$results/", mode: 'link'

    label 'goalign'

    input:
    file msa

    output:
    path "*_unique.phy"
    stdout

    script:
    """
    goalign dedup -i $msa -p | goalign clean sites --char GAP --cutoff 1.0 -p > ${msa.baseName}_unique.phy
    goalign stats length -p -i ${msa.baseName}_unique.phy > len.txt
    printf \$(cat len.txt)
    """
}

process InferRefTreeRAxML {
    publishDir "$results/", mode: 'link'

    label 'raxml'

    input:
    file msa
    val nboot

    output:
    path "reftree_collapsed_raxml.nw"
    path "reftree_original_raxml.nw"

    script:
    """
    raxml-ng --all --msa ${msa} --model GTR+G4 --threads ${task.cpus} --seed 123456 --bs-trees $nboot --bs-metric fbp
    mv ${msa}.raxml.bestTreeCollapsed reftree_collapsed_raxml.nw
    mv ${msa}.raxml.support reftree_original_raxml.nw
    """
}

process InferRefTreeIQTreeCollapse {
    publishDir "$results/", mode: 'link'

    label 'iqtree'

    input:
    file msa
    val nboot

    output:
    path "reftree_collapsed_iqtree.nw"

    script:
    """
    iqtree -s ${msa} -m GTR+G4 --seed 123456789 -T ${task.cpus} -b $nboot -czb
    mv ${msa}.treefile reftree_collapsed_iqtree.nw
    """
}

process InferRefTreeIQTreeNoCollapse {
    publishDir "$results/", mode: 'link'

    label 'iqtree'

    input:
    file msa
    val nboot

    output:
    path "reftree_original_iqtree.nw"

    script:
    """
    iqtree -s ${msa} -m GTR+G4 --seed 123456789 -T ${task.cpus} -b $nboot
    mv ${msa}.treefile reftree_original_iqtree.nw
    """
}

process InferRefTreePhyML {
    publishDir "$results/", mode: 'link'

    label 'phyml'

    input:
    path msa
    
    output:
    path "reftree_phyml.nw"
    
    script:
    """
    phyml -i ${msa} -m GTR -c 4 -a e -f e -d nt -o tlr -b 0 --r_seed 123456
    mv ${msa}_phyml_tree.txt reftree_phyml.nw
    """
}

process SeqBootPhyML {
    publishDir "$results/seqboot/", mode: 'link'

    label 'goalign'

    input:
    file msa
    val nboot

    output:
    file "bootphy_*"

    script:
    """
    goalign build seqboot -S -p -n $nboot -i ${msa} -o bootphy_ --seed 123456789
    """
}

process SeqBootTreesPhyML {
    label 'phyml'

    input:
    path msa
    
    output:
    path "${msa}_boottree_phyml.nw"
    
    script:
    """
    phyml -i ${msa} -m GTR -c 4 -a e -f e -d nt -o tlr -b 0 --r_seed 123456
    mv ${msa}_phyml_tree.txt ${msa}_boottree_phyml.nw
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
     file boot

     script:
     """
     gotree compute support fbp -i ${ref} -b ${boot} -o reftree_phyml_bootsupport.nw
     gotree collapse length -i ${boot} -l ${0.1/length} | gotree compute support fbp -i ${ref} -b - -o reftree_phyml_bootsupport_collapse.nw
     """
}

process SeqBootRAxML {
    publishDir "$results/seqboot/", mode: 'link'

    label 'goalign'

    input:
    file msa
    val nboot

    output:
    file "bootrax_*"

    script:
    """
    goalign build seqboot -S -p -n $nboot -i ${msa} -o bootrax_ --seed 987654321
    """
}

process SeqBootTreesRAxML {
    label 'raxmlshort'

    input:
    path msa
    
    output:
    path "${msa}_boottree_raxml.nw"
    
    script:
    """
    raxml-ng --msa ${msa} --model GTR+G4 --threads ${task.cpus} --seed 123456
    mv ${msa}.raxml.bestTree ${msa}_boottree_raxml.nw
    """
}

process SeqBootSupportsRAxML {
     publishDir "$results/", mode: 'link'

     label 'gotree'

     input:
     file ref
     file boot
     val length

     output:
     file "reftree_raxml_bootsupport.nw"
     file "reftree_raxml_bootsupport_collapse.nw"
     file boot

     script:
     """
     gotree compute support fbp -i ${ref} -b ${boot} -o reftree_raxml_bootsupport.nw
     gotree collapse length -i ${boot} -l ${0.1/length} | gotree compute support fbp -i ${ref} -b - -o reftree_raxml_bootsupport_collapse.nw
     """
}

workflow{
	msa = Channel.fromPath(params.align)
	nboot=params.nboot
	uniqueout = Unique(msa)
	dedup=uniqueout[0]
	length=uniqueout[1].map{it ->  Integer.parseInt(it.trim()) }

	refrax    = InferRefTreeRAxML(dedup,nboot)
	refiqtc   = InferRefTreeIQTreeCollapse(dedup,nboot)
	refiqtnc  = InferRefTreeIQTreeNoCollapse(dedup,nboot)

	phybootmsa   = SeqBootPhyML(dedup,nboot)
	raxbootmsa   = SeqBootRAxML(dedup,nboot)

	refraxnocoll = refrax[1]
	raxboottrees = SeqBootTreesRAxML(raxbootmsa.flatten())
	SeqBootSupportsRAxML(refraxnocoll, raxboottrees.collectFile(name: 'boot_rax.nw'),length)

	refphy = InferRefTreePhyML(dedup)
	phyboottrees = SeqBootTreesPhyML(phybootmsa.flatten())
	SeqBootSupportsPhyML(refphy, phyboottrees.collectFile(name: 'boot_phy.nw'),length)
}
