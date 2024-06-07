params.gisaidfasta="/pasteur/zeus/projets/p01/CNRVIR_bioinfoSEQ/Phylogenet/Omicron/wf_update/data/sequences_fasta_2022_03_06.tar.xz"
params.gisaidmeta="/pasteur/zeus/projets/p01/CNRVIR_bioinfoSEQ/Phylogenet/Omicron/wf_update/data/metadata_tsv_2022_03_06.tar.xz"
params.ids = "data/ids.txt"
params.results="results"

process selectByEPIISL {
   input:
   file ids
   file gisaidfasta
   file gisaidmeta

   output:
   file "selected.fasta"

   script:
   """
   selectByEPIISL.py -i $gisaidfasta -n $ids -m $gisaidmeta -o selected.fasta -t ${task.cpus}
   """
}

process refdata {

    output:
    path "sars-cov-2"

    script:
    """
    nextclade dataset get --name 'sars-cov-2' --output-dir 'sars-cov-2'
    """
}

process align {
    publishDir "${params.results}/", mode: "copy"

    input:
    path sequences
    path ref

    output:
    path "nextalign.aligned.fasta"

    script:
    """
    nextalign run -r ${ref}/reference.fasta  --output-all .  ${sequences}
    """
}

workflow{
    ids=file(params.ids)
    gisaidfasta=file(params.gisaidfasta)
    gisaidmeta=file(params.gisaidmeta)

    seq=selectByEPIISL(ids,gisaidfasta,gisaidmeta)
    ref=refdata()

    align(seq,ref)
}
