<H1>Code for splicing efficiency computation</H1>

The code consists of two phases, the first individual sites are computed and in the second the calculation is done on the level of transcripts and introns.
The first part needs to be exectud separately for each sample (but each sample can have 1+ BAM files if there are technical replicates etc)

<H2>First phase</H2>

`compute_sites_bed BED_FILE BAM_FILE(s) NAME FR_SECONDSTRAND`

<ul>
<li>`BED_FILE` is the BED file with the transcript annotations
<li>`SPECIES` is hg19 or mm9
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

`compute_introns_bed BED_FILE NAME(s) PHYLO_FILE 2BIT_FILE OUTNAME`

where:
<ul>
<li>`BED_FILE` is the BED file with the transcript annotations
<li>`NAME(s)` are the base names that were used for the compute_sites above
<li>`PHYLO_FILE` is the conservation bigWig file (phyloP), e.g., for hg38 <A HREF="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw">hg38.phyloP100way.bw</A> from UCSC
<li>`2BIT_FILE` is the genome sequence file. E.g. <A HREF="https://hgdownload.cse.ucsc.edu/goldenpath/hg38//bigZips/hg38.2bit">hg38.2bit</A>
<li>`OUTNAME` is the name of the output file
</ul>

This will generate `OUTNAME.intronStats.txt` containing the following information:
<ul>
<li>Intron name
<li>Intron position
<li>Intron length
<li>Host gene name
<li>Host manual annotation
<li>Conservation in the 5’ region
<li>Conservation in the 3’ region
<li>Average conservation
<li>Index of the starting exon (average over all transcripts containing this intron)
<li>Index of the end exon (average over all transcripts containing this intron)
<li>Relative index of the starting exon
<li>Relative index of the ending exon
<li>Does this intron overlap some other exon?
<li>Is the 5’ site derived from a transposable element?
<li>Is the 3’ site derived from a transposable element?
</ul>

and an `OUTNAME.transcriptStats.txt` file with the name of the transcript and the total numbers of spliced and unspliced reads respectively for each sample name. 

