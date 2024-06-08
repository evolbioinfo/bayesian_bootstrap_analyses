process RevCompAndExclude {
	label 'goalign'

	input:
        path ncbi
	path revcompid
	path excludelist

        output:
        path "*_revcomp.fasta"

        script:
        """
	goalign subset -i $ncbi -e -r -f $excludelist --unaligned > tmpncbi
	goalign subset -i tmpncbi -f $revcompid -e --unaligned | goalign revcomp --unaligned > revcomp
	goalign subset -r -i tmpncbi -f $revcompid -e --unaligned > tokeep
	cat revcomp tokeep > ${ncbi.baseName}_revcomp.fasta
        """
}

process ExtractNCBISequencesAndMetadata {
	publishDir 'results', mode: 'copy'
	label 'perl'
	
	input:
	path ncbifile
	path gbfile

	output:
	path "metadata.txt"
	path "segment_L.fasta"
	path "segment_M.fasta"
	path "segment_S.fasta"
	path "metadata_dates.txt"
	
	script:
	"""
	strain.pl $ncbifile $gbfile
	cat metadata.txt | cut -f 1,8| sort -u | awk 'BEGIN{IFS="\\t";OFS="\\t"}{n=split(\$2,a,"-");\$2=a[n];if(n>0){print \$0}}' > metadata_dates.txt
	"""
}

process MultipleSequenceAlignment {
	publishDir 'results', mode: 'copy'
	label 'mafft'

	input:
	path segment

	output:
	path "*aligned.fasta"

	script:
	"""
	mafft --auto --thread ${task.cpus} ${segment} > ${segment.baseName}_aligned.fasta
	"""
}

// Remove sequences that have more than 60% gaps
process CleanAlignment {
	publishDir 'results', mode: 'copy'
	label 'goalign'

	input:
	path msa

	output:
	path "*clean.fasta"

	script:
	"""
	goalign clean seqs --cutoff 0.6 -i $msa -o ${msa.baseName}_clean.fasta
	"""
}

process DeduplicateNames{
	publishDir 'results', mode: 'copy'
	label 'goalign'

	input:
	path align

	output:
	path "*_dedup.fasta"

	script:
	"""
	goalign reformat fasta --unaligned -i $align --ignore-identical 0 -o ${align.baseName}_dedup.fasta
	"""
}

process phylogeny{
	publishDir 'results', mode: 'copy'
	label 'phylo'

	input:
	path align
	
	output:
	path "*treefile"
	
	script:
	"""
	iqtree -nt ${task.cpus} -s ${align} -m GTR+G4+FO  --seed 123456789
	"""
}

workflow {
	 ncbi = Channel.fromPath(params.ncbi)
	 gb = Channel.fromPath(params.gb)
	 revcompid = file("data/revcomp.txt")
	 excludelist = file("data/excludelist.txt")
	 revcomp = RevCompAndExclude(ncbi,revcompid,excludelist)
	 seqmeta = ExtractNCBISequencesAndMetadata(revcomp,gb)

	 sequencesL = seqmeta[1].collectFile(name: "full_L.fasta")
	 sequencesM = seqmeta[2].collectFile(name: "full_M.fasta")
	 sequencesS = seqmeta[3].collectFile(name: "full_S.fasta")

	 msa = CleanAlignment(MultipleSequenceAlignment(DeduplicateNames(sequencesL.mix(sequencesM).mix(sequencesS))))
	 tree = phylogeny(msa)
}
