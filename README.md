<H1>Code for splicing efficiency computation</H1>

The code consists of two phases, the first individual sites are computed and in the second the calculation is done on the level of transcripts and introns.
The first part needs to be exectud separately for each sample (but each sample can have 1+ BAM files if there are technical replicates etc)

<H2>First phase</H2>

`java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar scripts.lincs.util.QuantifySplicingEfficiency compute_sites_bed BED_FILE BAM_FILE(s) NAME FR_SECONDSTRAND`

<ul>
<li>`BED_FILE` is the BED file with the transcript annotations
<li>`BAM_FILE(s)` are .bam file names (comma separated)
<li>`NAME` is the base name for output files
<li>`LIBRARY_TYPE` is one of the following: `UNSTRANDED`, `FF_FIRSTSTRAND` - paired-end and both reads are on the + strand (ENCODE data only) `FR_FIRSTSTRAND` - regular dUTP/TruSeq strand specific, `FR_SECONDSTRAND` - single-end stranded or other protocols where the first read is on the + strand
</ul>

This pase will generate `NAME.5P.sites.txt` and `NAME.3P.sites.txt` files. Each has the following columns:
<ul>
<li>Chromosome
<li>Strand
<li>Position
<li>Number of reads that go through this position without splicing
<li>Number of reads that go through this position with splicing
<li>Currently 0
<li>For the spliced reads - the positions to which the reads splice to (e.g., the 3' splice site positions in the 5P file) and the number of reads correspond to each target position
</ul>

<H2>Second phase</H2>
In this phase multiple samples can be processed together

`java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar:jar/bigWig.jar:jar/log4j-1.2.15.jar scripts.lincs.global.LincRNASplicingAnalysis compute_introns_bed [-minReads MIN_READS] -spliceDir SPLICE_DIRECTORY BED_FILE NAME(s) OUTNAME PHYLO_FILE 2BIT_FILE`

where:
<ul>
<li>`MIN_READS` is the minimal number of reads considered per splice site in the transcript-level computation (default 20)
<li>`SPLICE_DIRECTORY` is the location of the directory with the files `all5ss_GTC_only_with_maxent.txt` and `PSSM3ss`
<li>`BED_FILE` is the BED file with the transcript annotations
<li>`NAME(s)` are the base names that were used for the compute_sites above
<li>`PHYLO_FILE` is the conservation bigWig file (phyloP), e.g., for hg38 <A HREF="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw">hg38.phyloP100way.bw</A> from UCSC
<li>`2BIT_FILE` is the genome sequence file. E.g. <A HREF="https://hgdownload.cse.ucsc.edu/goldenpath/hg38//bigZips/hg38.2bit">hg38.2bit</A>
<li>`OUTNAME` is the name of the output file
</ul>

This will generate `OUTNAME.intronStats.txt` containing the following information (for all the unique introns in the BED file):
<ul>
<li>Intron name
<li>Intron position
<li>Intron length
<li>Host gene name
<li>Average number of exons in the transcripts containing this intron
<li>Sequence of the 5' splice-site
<li>Sequence of the 3' splice-site
<li>Is the 5' splice-site canonical?
<li>Is the 3' splice-site canonical?
<li>Senepathy score of the 5' splice-site
<li>Delta-G score of the 5' splice-site
<li>Max entropy score of the 5' splice-site
<li>PSSM Score of the 3' splice-site
<li>Conservation in the 5’ region
<li>Conservation in the 3’ region
<li>Average conservation
<li>Index of the starting exon (average over all transcripts containing this intron)
<li>Index of the end exon (average over all transcripts containing this intron)
<li>Relative index of the starting exon
<li>Relative index of the ending exon
<li>Does this intron overlap some other intron?
<li>Does this intron overlap some other exon?
<li>Is the 5’ site derived from a transposable element?
<li>Is the 3’ site derived from a transposable element?
<li>Per sample: (just for introns with at least MIN_READs on each side)
  <ul>
    <li>Number of reads covering the 5P of the intron (5P)
    <li>Number of reads covering the 3P of the intron (3P)
    <li>Number of reads covering the spliced reads for the intron 
    <li>Splicing efficiency (0.5+Spliced)/(0.5+Spliced+Unspliced_5P+Unspliced_3P)
    <li>Min-based splicing efficiency (0.5+Spliced)/(0.5+Spliced+2*min(Unspliced_5P,mUnspliced_3P))
    <li>Pratio - Unspliced_5P/Unspliced_3P
   </ul>
<li>Total number of spliced reads accross all the samples  
</ul>

and an `OUTNAME.transcriptStats.txt` file with:
<ul>
<li> Transcript ID
<li> Number of exons
<li> For each sample:
  <ul>
    <li>Number of introns with information (above MIN_READS for both 5P and 3P)</li>
    <li>Spliced reads (combining all introns, even those with few reads)</li>
    <li>Unspliced reads (combining all introns, even those with few reads)</li>
    <li>Total reads (combining all introns, even those with few reads)</li>
    <li>Efficiency of the intron with the lowest efficiency (only among introns with at least MIN_READS for both 5P and 3P)</li>
  </ul>
 </ul>
 
 <h3>Sample commands</H3>
 <p>
 <ul>
 <li>
`java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar scripts.lincs.util.QuantifySplicingEfficiency compute_sites_bed data/MANE.GRCh38.v1.0.refseq_genomic.bed data/Cyto.bam Cyto FR_SECONDSTRAND`</li>
 <li>`java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar scripts.lincs.util.QuantifySplicingEfficiency compute_sites_bed data/MANE.GRCh38.v1.0.refseq_genomic.bed data/Nuc.bam Nuc FR_SECONDSTRAND`</li>
 <li>`java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar:jar/bigWig.jar:jar/log4j-1.2.15.jar scripts.lincs.global.LincRNASplicingAnalysis compute_introns_bed -spliceDir splice/ data/MANE.GRCh38.v1.0.refseq_genomic.bed Nuc,Cyto full data/hg38.phyloP100way.bw data/hg38.2bit`</li>
 </p>
    

