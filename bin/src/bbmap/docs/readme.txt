BBMap readme by Brian Bushnell
Last updated February 18, 2014.
Please contact me at bbushnell@lbl.gov if you have any questions or encounter any errors.
BBMap is free to use for noncommercial purposes, and investigators are free to publish results derived from the program, as long as the source code is not published or modified.

This is the official release of BBMAP, version 31.x


NOTE:
Don't recompile unless you run into version problems (such as trying to run with java 1.6).  But if you must, then -
To recompile, run this:
javac -J-Xmx128m align2/*.java jgi/*.java driver/*.java fileIO/*.java dna/*.java org/apache/tools/bzip2/*.java
(*** NOTE - due to email file size limits, I may have omitted the .java files if you received this via email, in which case you can't recompile. ***)


Basic Syntax:

(using shellscript, on Genepool, which autodetects RAM to set -Xmx parameter)
To index:
bbmap.sh build=1 ref=<reference.fa>
To map:
bbmap.sh build=1 in=<reads.fq> out=<mapped.sam>

(without shellscript)
To index:
java -ea -Xmx31g -cp <PATH> align2.BBMap build=1 ref=<reference.fa>
To map:
java -ea -Xmx31g -cp <PATH> align2.BBMap build=1 in=<reads.fq> out=<mapped.sam>

...where "<PATH>" should indicate the path to the directory containing all the source code directories; e.g. "/global/projectb/sandbox/gaag/bbtools/current"

Please note, the reference is only needed for building the index the first time; subsequently, just specify the build number which corresponds to that reference.
So for example the first time you map to e.coli you might specify "ref=ecoli_reference.fa build=3"; after that, just specify "build=3".
The index files would then be stored in ./ref/genome/3/ and ./ref/index/3/
Also, the -Xmx parameter should specify approximately 85% of the physical memory of the target machine; so, 21G for a 24GB node.  The process needs approximately 8 bytes per reference base (plus a several hundred MB overhead).


Advanced Syntax:


Indexing Parameters (required when building the index):
path=<.>        	Base directory to store index files.  Default is the local directory.  The index will always be placed in a subdirectory "ref".
ref=<ref.fasta> 	Use this file to build the index.  Needs to be specified only once; subsequently, the build number should be used.
build=<1>		Write the index to this location (build=1 would be stored in /ref/genome/1/ and /ref/index/1/).  Can be any integer.  This parameter defaults to 1, but using additional numbers allows multiple references to be indexed in the same directory.
k=<13>          	Use length 13 kmers for indexing.  Suggested values are 9-15, with lower typically being slower and more accurate.  13 is usually optimal.  14 is better for RNA-SEQ and very large references >4GB; 12 is better for PacBio and cross-species mapping.
midpad=<300>		Put this many "N" in between scaffolds when making the index.  300 is fine for metagenomes with millions of contigs; for a finished genome like human with 25 scaffolds, this should be set to 100000+ to prevent cross-scaffold mapping.
startpad=<8000> 	Put this many "N" at the beginning of a "chrom" file when making index.  It's best if this is longer than your longest expected read.
stoppad=<8000>		Put this many "N" at the end of a "chrom" file when making index.  It's best if this is longer than your longest expected read.
minscaf=<1>		Do not include scaffolds shorter than this when generating index.  Useful for assemblies with millions of fairly worthless unscaffolded contigs under 100bp.  There's no reason to make this shorter than the kmer length.


Input Parameters:
path=<.>		Base directory to read index files.
build=<1>		Use the index at this location (same as when indexing).
in=<reads.fq>		Use this as the input file for reads.  Also accepts fasta.  "in=sequential length=200" will break a genome into 200bp pieces and map them to itself.  "in=stdin" will accept piped input.  The format of piped input can be specified with e.g. "in=stdin.fq.gz" or "in=stdin.fa"; default is uncompressed fastq.
in2=<reads2.fq> 	Run mapping paired, with reads2 in the file "reads2.fq"
			NOTE:  As a shorthand, "in=reads#.fq" is equivalent to "in=reads1.fq in2=reads2.fq"
interleaved=<auto>	Or "int". Set to "true" to run mapping paired, forcing the reads to be considered interleaved from a single input file.  By default the reader will try to determine whether a file is interleaved based on the read names; so if you don't want this, set interleaved=false.
qin=<auto>       	Set to 33 or 64 to specify input quality value ASCII offset.
fastareadlen=<500>	If fasta is used for input, breaks the fasta file up into reads of about this length.  Useful if you want to map one reference against another, since BBMap currently has internal buffers limited to 500bp.  I can change this easily if desired.
fastaminread=<1>	Ignore fasta reads shorter than this.  Useful if, say, you set fastareadlen=500, and get a length 518 read; this will be broken into a 500bp read and an 18bp read.  But it's not usually worth mapping the 18bp read, which will often be ambiguous.
fakequality=<-1>	Set to a positive number 1-50 to generate fake quality strings for fasta input reads.  Less than one turns this function off.
blacklist=<a.fa,b.fa>	Set a list of comma-delimited fasta files.  Any read mapped to a scaffold name in these files will be considered "blacklisted" and can be handled differently by using the "outm", "outb", and "outputblacklisted" flags.  The blacklist fasta files should also be merged with other fasta files to make a single combined fasta file; this combined file should be specified with the "ref=" flag when indexing.
touppercase=<f>		Set true to convert lowercase read bases to upper case.  This is required if any reads have lowercase letters (which real reads should never have).


Sampling Parameters:
reads=<-1>		Process at most N reads, then stop.  Useful for benchmarking.  A negative number will use all reads.
samplerate=<1.0>	Set to a fraction of 1 if you want to randomly sample reads.  For example, samplerate=0.25 would randomly use a quarter of the reads and ignore the rest.  Useful for huge datasets where all you want to know is the % mapped.
sampleseed=<1>		Set to the RNG seed for random sampling.  If this is set to a negative number, a random seed is used; for positive numbers, the number itself is the seed.  Since the default is 1, this is deterministic unless you explicitly change it to a negative number.	
idmodulo=<1>		Set to a higher number if you want to map only every Nth read (for sampling huge datasets).


Mapping Parameters:
fast=<f>		The fast flag is a macro.  It will set many other paramters so that BBMap will run much faster, at slightly reduced sensitivity for most applications.  Not recommended for RNAseq, cross-species alignment, or other situations where long deletions or low identity matches are expected.
minratio=<0.56>		Alignment sensitivity as a fraction of a read's max possible mapping score.  Lower is slower and more sensitive but gives more false positives.  Ranges from 0 (very bad alignment) to 1 (perfect alignment only).  Default varies between BBMap versions. 
minidentity=<>		Or "minid".  Use this flag to set minratio more easily.  If you set minid=0.9, for example, minratio will be set to a value that will be approximately equivalent to 90% identity alignments.
minapproxhits=<1>	Controls minimum number of seed hits to examine a site.  Higher is less accurate but faster (on large genomes).  2 is maybe 2.5x as fast and 3 is maybe 5x as fast on a genome with of gigabases.  Does not speed up genomes under 100MB or so very much.
padding=<4>		Sets extra padding for slow-aligning.  Higher numbers are more accurate for indels near the tips of reads, but slower.
tipsearch=<100>		Controls how far to look for possible deletions near tips of reads by brute force.  tipsearch=0 disables this function.  Higher is more accurate.
maxindel=<16000>	Sets the maximum size of indels allowed during the quick mapping phase.  Set higher (~100,000) for RNA-SEQ and lower (~20) for large assemblies with mostly very short contigs.  Lower is faster.
strictmaxindel=<f>	Set to true to disallow mappings with indels longer than maxindel.  Alternately, for an integer X, 'strictmaxindel=X' is equivalent to the pair of flags 'strictmaxindel=t maxindel=X'.
pairlen=<32000>  	Maximum distance between mates allowed for pairing.
requirecorrectstrand=<t>	Or "rcs".  Requires correct strand orientation when pairing reads.  Please set this to false for long mate pair libraries!
samestrandpairs=<f>	Or "ssp".  Defines correct strand orientation when pairing reads.  Default is false, meaning opposite strands, as in Illumina fragment libraries.  "ssp=true" mode is not fully tested.
killbadpairs=<f>	Or "kbp".  When true, if a read pair is mapped with an inappropriate insert size or orientation, the read with the lower mapping quality is marked unmapped.
rcompmate=<f>		***TODO*** Set to true if you wish the mate of paired reads to be reverse-complemented prior to mapping (to allow better pairing of same-strand pair libraries).
kfilter=<-1>		If set to a positive number X, all potential mapping locatiosn that do not have X contiguous perfect matches with the read will be ignored.  So, reads that map with "kfilter=51" are assured to have at least 51 contiguous bases that match the reference.  Useful for mapping to assemblies generated by a De Bruijn graph assembly that used a kmer length of X, so that you know which reads were actually used in the assembly.
threads=<?>		Or "t".  Set number of threads.  Default is # of logical cores.  The total number of active threads will be higher than this, because input and output are in seperate threads.
perfectmode=<f>		Only accept perfect mappings.  Everything goes much faster.  
semiperfectmode=<f>	Only accept perfect or "semiperfect" mappings.  Semiperfect means there are no mismatches of defined bases, but up to half of the reference is 'N' (to allow mapping to the edge of a contig).
rescue=<t>		Controls whether paired may be rescued by searching near the mapping location of a mate.  Increases accuracy, with usually a minor speed penalty.
expectedsites=<1>	For BBMapPacBioSkimmer only, sets the expected number of correct mapping sites in the target reference.  Useful if you are mapping reads to other reads with some known coverage depth.
msa=<>			Advanced option, not recommended.  Set classname of MSA to use.
bandwidth=0		Or "bw".  When above zero, restricts alignment band to this width.  Runs faster, but with reduced accuracy for reads with many or long indels.
bandwidthratio=0	Or "bwr".  When above zero, restricts alignment band to this fraction of a read's length.  Runs faster, but with reduced accuracy for reads with many or long indels.
usequality=<t>		Or "uq".  Set to false to ignore quality values when mapping.  This will allow very low quality reads to be attempted to be mapped rather than discarded.
keepbadkeys=<f>		Or "kbk".  With kbk=false (default), read keys (kmers) have their probability of being incorrect evaluated from quality scores, and keys with a 94%+ chance of being wrong are discarded.  This increases both speed and accuracy.


Output Parameters:
out=<outfile.sam>	Write output to this file.  If out=null, output is suppressed.  If you want to output paired reads to paired files, use a "#" symbol, like out=mapped#.sam.  Then reads1 will go to mapped1.sam and reads2 will go to mapped2.sam. (NOTE: split output currently diabled for .sam format, but allowed for native .txt format).  To print to standard out, use "out=stdout"
outm=<>			Write only mapped reads to this file (excluding blacklisted reads, if any).
outu=<>			Write only unmapped reads to this file.
outb=<>			Write only blacklisted reads to this file.  If a pair has one end mapped to a non-blacklisted scaffold, it will NOT go to this file. (see: blacklist)
out2=<>			If you set out2, outu2, outm2, or outb2, the second read in each pair will go to this file.  Not currently allowed for SAM format, but OK for others (such as fasta, fastq, bread).
overwrite=<f>		Or "ow".  Overwrite output file if it exists, instead of aborting.
ambiguous=<best>	Or "ambig". Sets how to handle ambiguous reads.  "first" or "best" uses the first encountered best site (fastest).  "all" returns all best sites.  "random" selects a random site from all of the best sites (does not yet work with paired-ends).  "toss" discards all sites and considers the read unmapped (same as discardambiguous=true).  Note that for all options (aside from toss) ambiguous reads in SAM format will have the extra field "XT:A:R" while unambiguous reads will have "XT:A:U".
ambiguous2=<best>	Or "ambig2". Only for splitter mode.  Ambiguous2 strictly refers to any read that maps to more than one reference set, regardless of whether it has multiple mappings within a reference set.  This may be set to "best" (aka "first"), in which case the read will be written only to the first reference to which it has a best mapping; "all", in which case a read will be written to outputs for all references to which it maps; "toss", in which case it will be considered unmapped; or "split", in which case it will be written to a special output file with the prefix "AMBIGUOUS_" (one per reference).
outputunmapped=<t>	Outputs unmapped reads to primary output stream (otherwise they are dropped).
outputblacklisted=<t>	Outputs blacklisted reads to primary output stream (otherwise they are dropped).
cigar=<t>		Generate cigar strings (for bread format, this means match strings).  cigar=false is faster.  "cigar=" is synonymous with "match=".  This must be enabled if match/insertion/deletion/substitution statistics are desired, but the program will run faster with cigar strings disabled.
keepnames=<f>		Retain original names of paired reads, rather than ensuring both reads have the same name when written in sam format by renaming read2 to the same as read1.  If this is set to true then the output may not be sam compliant.
mdtag=<f>		Generate MD tags for SAM files.  Requires that cigar=true.  I do not recommend generating MD tags for RNASEQ or other data where long deletions are expected because they will be incredibly long.
xstag=<f>		Generate XS (strand) tags for Cufflinks.  This should be used with a stranded RNA-seq protocol.
xmtag=<t>		Generate XM tag.  Indicates number of best alignments.
intronlen=<999999999>	Set to a lower number like 10 to change 'D' to 'N' in cigar strings for deletions of at least that length.  This is used by Cufflinks; 'N' implies an intron while 'D' implies a deletion, but they are otherwise identical.
stoptag=<t>		Allows generation of custom SAM tag YS:i:<read stop location>
idtag=<t>		Allows generation of custom SAM tag YI:f:<percent identity>
ordered=<f>		Set to true if you want reads to be output in the same order they were input.  This takes more memory, and can be slower, due to buffering in multithreaded execution.  Not needed for singlethreaded execution.
ziplevel=<2>		Sets output compression level, from 1 (fast) to 9 (slow).  I/O is multithreaded, and thus faster when writing paired reads to two files rather than one interleaved file.
nodisk=<f>		"true" will not write the index to disk, and may load slightly faster.   Prevents collisions between multiple bbmap instances writing indexes to the same location at the same time.
usegzip=<t>		If gzip is installed, output file compression is done with a gzip subprocess instead of with Java's native deflate method.  Can be faster when set to true.  The output file must end in a compressed file extension for this to have effect. *Temporarily disabled due to bug in UGE.
usegunzip=<t>		If gzip is installed, input file decompression is done with a gzip subprocess instead of with Java's native inflate method.  Can be faster when set to true.  *Temporarily disabled due to bug in UGE.
samversion=<1.3>	SAM specification version. Set to 1.3 for cigar strings with 'M' or 1.4 for cigar strings with '=' and 'X'.  Default is currently 1.3 because samtools 0.1.18 and earlier is incompatible with sam format version 1.4.
bamscript=<filename>	(bs for short) Writes a shell script to <filename> with the command line to translate the sam output of BBMap into a sorted bam file, assuming you have samtools in your path.
maxsites=<5>		Sets maximum alignments to print per read, if secondary alignments are allowed.  Currently secondary alignments may lack cigar strings.
secondary=<f>		Print secondary alignments.
quickmatch=<f>		Generate cigar strings during the initial alignment (before the best site is known).  Currently, this must be enabled to generate cigar strings for secondary alignments.  It increases overall speed but may in some very rare cases yield inferior alignments due to less padding.
local=<f>          	Output local alignments instead of global alignments.  The mapping will still be based on the best global alignment, but the mapping score, cigar string, and mapping coordinate will reflect a local alignment (using the same affine matrix as the global alignment).
sortscaffolds=<f>	Sort scaffolds alphabetically in SAM headers to allow easier comparisons with Tophat (in cuffdif, etc).  Default is in same order as source fasta.
trimreaddescriptions=<f>	(trd) Truncate read names at the first whitespace, assuming that the remaineder is a comment or description.


Statistics Parameters:
showprogress=<f>	Set to true to print out a '.' once per million reads processed.  You can also change the interval with e.g. showprogress=20000.
qhist=<filename>	Output a per-base average quality histogram to <filename>.
mhist=<filename>	Output a per-base match histogram to <filename>.  Requires cigar strings to be enabled.  The columns give fraction of bases at each position having each match string operation: match, substitution, deletion, insertion, N, or other.
ihist=<filename>	Output a per-read-pair insert size histogram to <filename>.
scafstats=<filename>	Track mapping statistics per scaffold, and output to <filename>.
refstats=<filename>	For BBSplitter, enable or disable tracking of read mapping statistics on a per-reference-set basis, and output to <filename>.
verbosestats=<0>	From 0-3; higher numbers will print more information about internal program counters.


Trimming Parameters:
qtrim=<f>		Options are false, left, right, or both.  Allows quality-trimming of read ends before mapping.
			false: Disable trimming.
			left (l): Trim left (leading) end only.
			right (r): Trim right (trailing) end only.  This is the end with lower quality many platforms.
			both (lr): Trim both ends.
trimq=<5>		Set the quality cutoff.  Bases will be trimmed until there are 2 consecutive bases with quality GREATER than this value; default is 5.  If the read is from fasta and has no quality socres, Ns will be trimmed instead, as long as this is set to at least 1.
untrim=<f>		Untrim the read after mapping, restoring the trimmed bases.  The mapping position will be corrected (if necessary) and the restored bases will be classified as soft-clipped in the cigar string.

Java Parameters:
-Xmx       		If running from the shellscript, include it with the rest of the arguments and it will be passed to Java to set memory usage, overriding the shellscript's automatic memory detection.  -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max allowed is typically 85% of physical memory.



Splitting Parameters:
The splitter is invoked by calling bbsplit.sh (or align2.BBSplitter) instead of bbmap.sh, for the indexing phase.  It allows combining multiple references and outputting reads to different files depending on which one they mapped to best.  The order in which references are specified is important in cases of ambiguous mappings; when a read has 2 identically-scoring mapping locations from different references, it will be mapped to the first reference.
All parameters are the same as BBMap with the exception of the ones listed below.  You can still use "outu=" to capture unmapped reads.
ref_<name>=<fasta files>	Defines a named set of organisms with a single fasta file or list.  For example, ref_a=foo.fa,bar.fa defines the references for named set "a"; any read that maps to foo.fasta or bar.fasta will be considered a member of set a.
out_<name>=<output file>	Sets the output file name for reads mapping to set <name>.  out_a=stuff.sam would capture all the reads mapping to ref_a.
basename=<example%.sam>		This shorthand for mass-specifying all output files, where the % symbol is a wildcard for the set name.  For example, "ref_a=a.fa ref_b=b.fa basename=mapped_%.sam" would expand to "ref_a=a.fa ref_b=b.fa out_a=mapped_a.sam out_b=mapped_b.sam"
ref=<fasta files>		When run through the splitter, this is shorthand for writing a bunch of ref_<name> entries.  "ref=a.fa,b.fa" would expand to "ref_a=a.fa ref_b=b.fa".


Formats and Extensions
.gz,.gzip,.zip,.bz2	These file extensions are allowed on input and output files and will force reading/writing compressed data.
.fa,.fasta,.txt,.fq,.fastq	These file extensions are allowed on input and output files.  Having one is REQUIRED.  So, reads.fq and reads.fq.zip are valid, but reads.zip is NOT valid.  Note that outputting in fasta or fastq will not retain mapping locations.
.sam			This is only allowed on output files.
.bam			This is allowed on output files if samtools is installed.  Beware of memory usage; samtools will run in a subprocess, and it can consume over 1kb per scaffold of the reference genome.


Different versions:
BBMap			Fastest version.  Finds single best mapping location.
BBMapAcc		Slower and more accurate.  Finds single best mapping location.  May currently be broken.
BBMapPacBio		Optimized for PacBio's error profile (more indels, fewer substitutions).  Finds single best mapping location.  PacBio reads should be in fasta format.
BBMapPacBioSkimmer	Designed to find ALL mapping locations with alignment score above a certain threshold; also optimized for Pac Bio reads.
BBSplit			Uses BBMap or BBMapPacBio to map to multiple references simultaneously, and output the reads to the file corresponding to the best-matching reference.  Designed to split metagenomes or contaminated datasets prior to assembly.



Notes.
File types are autodetected by parsing the filename.  So you can name files, say, out.fq.gz or out.fastq.gz or reads1.fasta.bz2 or data.sam and it will work as long as the extensions are correct.


Change Log:

v31.
TODO:  Change pipethreads to redirects (where possible), and hash pipethreads by process, not by filename.
TODO:  Improve scoring function by using gembal distribution and/or accounting for read length.
TextStreamWriter was improperly testing for output format 'other'.  Noted by Brian Foster.
Fixed bug for read stream 2 in RTextOutputStream3.  Found by Brian Foster.
Fixed bug in MateReadsMT creating an unwanted read stream 2.  Found by Brian Foster.
TrimRead.testOptimal() mode added, and made default when quality trimming is performed; old mode can be used with 'otf=f' flag.
Fixed a couple cases where output file format was set to "ordered" even though the process was singlethreaded; this had caused an out-of-memory crash noted by Bill A.
Changed shellscripts of MapPacBio classes to remove "interleaved=false" term.
Reduced Shared.READ_BUFFER_LENGTH from 500 to 200 and Shared.READ_BUFFER_MAX_DATA from 1m to 500k, to reduce ram usage of buffers.
Noticed small bug in trimming; somehow a read had a 'T' with quality 0, which triggered assertion error.  I disabled the assertion but I'm not sure how it happened.
Fixed bug in which pigz was not used to decompress fasta files.
All program message information now defaults to stderr.

v30.
Disabled compression/decompression subprocesses when total system threads allowed is less than 3.
Fixed assertion error in calcCorrectness in which SiteScores are not necessarily sorted if AMBIGUOUS_RANDOM=true.  Noted by Brian Foster.
Fixed bug in toLocalAlignment with respect to considering XY as insertions, not subs.  
TODO: XY should be standardized as substitutions.
Added scarf input support. Requested by Alex Copeland.
TODO: Allow sam input with interleaved flag.
TODO: Make pigz a module dependency or script load.
Fixed bug with nodisk mode dropping the name of the first scaffold of every 500MB chunk after the first.  Noted by Vasanth Singan.
Overhaul of I/O channel creation.  Sequence files are now initialized with a FileFormat object which contains information about the format, permission to overwrite, etc.
Increased limit of number of index threads in Windows in nodisk mode (since disk fragmentation is no longer relevant).
Renamed Read.list to sites; added Read.topSite() and Read.numSites(); replaced many instances of things like "r.sites!=null && !r.sites.isEmpty()"
Refactored to put Read and all read-streaming I/O classes in 'stream' package.
Moved kmer hashing and indexing classes to kmer package.
Moved Variation, subclasses, and related classes to var package.
Moved FastaToChrom and ChromToFasta to dna package.
Moved pacbio error correction classes to pacbio package.
Removed stack, stats, primes, and other packages; prefixed all unused pacakges with z_.
TODO: Sites failing Data.isSingleScaffold() test should be clipped, not discarded.
RandomReads3 no longer adds /1 and /2 to paired fastq read names by default (can be enabled with 'addpairnum' flag).
Added "inserttag" flag; adds the insert size to sam output.
Fixed insert size histogram anomaly.  There was a blip at insert==(read1.length+read2.length) because the algorithm used to calculate insert size was different for reads that overlap and reads that don't overlap.
Skimmer now defaults to cigar=true.
Added maxindel1 and maxindel2 (or maxindelsum) flags.
Removed OUTER_DIST_MULT2 because it caused assertion errors when different from OUTER_DIST_MULT; changed OUTER_DIST_MULT from 15 to 14.
Added shellscript for skimmer, bbmapskimmer.sh
TODO:  Document above changes to parameters.



v29.
New version since major refactoring.
Added FRACTION_GENOME_TO_EXCLUDE flag (fgte).  Setting this lower increases sensitivity at expense of speed.  Range is 0-1 and default is around 0.03.
Added setFractionGenometoExclude() to Skimmer index.
LMP librares were not being paired correctly.  Now "rcs=f" may be used to ignore orientation when pairing.  Noted by Kurt LaButti.
Allocating memory to alignment score matrices caused uncaught out-of-memory error on low-memory machines, resulting in a hang.  This is now caught and results in an exit.  Noted by Alicia Clum.
GPINT machines are now detected and restricted to 4 threads max.  This helps prevent out-of-memory errors with PacBio mode.
Fixed sam output bug in which an unmapped read would get pnext of 0 rather than 1 when its mate mapped off the beginning of a scaffold.  Noted by Rob Egan.
Added memory test prior to allocating mapping threads.  Thread count will be reduced if there is not enough memory.  This is to address the issue noted by James Han, in which the PacBio versions would crash after running out of memory on low-memory nodes.
TODO:  Detect and prevent low-memory crashes while loading the index by aborting.
Fixed assertion error caused by strictmaxindel mode (noted by James Han).
Added flag "trd" (trimreaddescriptions) which truncates read names at the first whitespace.
Added "usequality/uq" flag to turn on/off usage of quality information when mapping.  Requested by Rob Egan.
Added "keepbadkeys/kbk" flag to prevent discarding of keys due to low quality.  Requested by Rob Egan.
Fixed crash with very long reads and very small kmers due to exceeding length of various kmer array buffers.
Avg Initial Sites and etc no longer printed for read 2 data.
TODO:  Support for selecting long-mate-pair orientation has been requested by Alex C.
Fixed possible bug in read trimming when the entire read was below the quality threshold.
Fixed trim mode bug: "trim=both" was only trimming the right side.  "qtrim" is also now an alias for "trim".
Fixed bug in ConcurrentGenericReadInputStream causing an incorrect assertion error for input in paired files and read sampling.  Found by Alex Copeland.
TODO:  Found some bz2 files that cannot be read.  Decompressing them and recompressing them fixes the problem, but it's unclear why other programs can read them.  Found by Alex Copeland.
Added insert size histogram: ihist=<file>
Added "machineout" flag for machine-readable output stats.
TODO: reads_B1_100000x150bp_0S_0I_0D_0U_0N_interleaved.fq.gz (ecoli) has 0% rescued for read1 and 0.7% rescued for read 2.  After swapping r1 and r2, .664% of r2 is rescued and .001% of r1 is rescued.  Why are they not symmetric?
Added 'slow' flag to bbmap for increased accuracy.  Still in progress.
Added MultiStateAligner11ts to MSA minIdToMinRatio().
Changed the way files are tested for permission to write (moved to Tools).
Fixed various places in which version string was parsed as an integer.
Added test for "help" and "version" flags.
Fixed bug in testing for file existence; noted by Bryce Foster.
Fixed issue with scaffold names not being trimmed on whitespace boundaries when 'trd=t'.  Noted by Rob Egan.
Added pigz (parallel gzip) support, at suggestion of Rob Egan.
Improved support for subprocesses and pipethreads; they are now automatically killed when not needed, even if the I/O stream is not finished.  This allows gunzip/unpigz when a file is being partially read.
Added shellscript test for the hostname 'gpint'; in that case, memory will be capped at 4G per process.
Changed the way cris/ros are shut down.  All must now go through ReadWrite.closeStreams()
TODO:  Force rtis and tsw to go through that too.
TODO:  Add "Job.fname" field.
Made output threads kill processes also.
Modified TrimRead to require minlength parameter.
Fixed a bug with gathering statistics in BBMapPacBioSkimmer (found by Matt Scholz).
Fixed a bug in which reads with match string containing X/Y were not eligible to be semiperfect (Found by Brian Foster).
Fixed a bug related to improving the prior fix; I had inverted an == operator (Found by Brian Foster).
Added SiteScore.fixXY(), a fast method to fix reads that go out-of-bounds during alignment.  Unfinished; score needs to be altered as a result.
Added "pairsonly" or "po" flag.  Enabling it will treat unpaired reads as unmapped, so they will be sent to 'outu' instead of 'outm'.  Suggested by James Han and Alex Copeland.
Added shellscript support for java -Xmx flag (Suggested by James Han).
Changed behavior: with 'quickmatch' enabled secondary sites will now get cigar strings (mostly, not all of them).
"fast" flag now enables quickmatch (50% speedup in e.coli with low-identity reads).  Very minor effect on accuracy.
Fixed bug with overflowing gref due GREFLIMIT2_CUSHION padding.  Found by Alicia Clum.
Fixed bug in which writing the index would use pigz rather than native gzip, allowing reads from scaffolds.txt.gz before the (buffered) writing finished.  Rare race condition.  Found by Brian Foster.
Fixed stdout.fa.gz writing uncompressed via ReadStreamWriter.
Added "allowSubprocess" flag to all constructors of TextFile and TextStreamWriter, and made TextFile 'tryAllExtensions' flag the last param.
allowSubprocess currently defaults to true for ByteFiles and ReadInput/Output Streams.
TODO: TextFile and TextStreamWriter (and maybe others?) may ignore ReadWrite.killProcess().
TODO: RTextOutputStream3 - make allowSubprocess a parameter
TODO: Assert that first symbol of reference fasta is '>' to help detect corrupt fastas.
Improved TextStreamWriter, TextFile, and all ReadStream classes usage of ReadWrite's InputStream/OutputStream creation/destruction methods.
All InputStream and OutputStream creation/destruction now has an allowSubprocesses flag.
Added verbose output to all ReadWrite methods.
Fixed bug in which realigned SiteScores were not given a new perfect/semiperfect status.  Noted by Brian Foster and Will Andreopoulos.


v28.
New version because the new I/O system seems to be stable now.
Re-enabled bam input/output (via samtools subprocess).  Lowered shellscript memory from 85% to 84% to provide space for samtools.
Added "-l" to "#!/bin/bash" at top.  This may make it less likely for the environment to be messed up.  Thanks to Alex Boyd for the tip.
Addressed potential bug in start/stop index padding calculation for scaffolds that began or ended with non-ACGT bases.
Made superclass for Index.
Made superclass for BBMap.
Removed around 5000 lines of code as a result of dereplication into superclasses.
Added MultiStateAligner11ts, which uses arrays for affine transform instead of if blocks.  Changing insertions gave a ~5% speedup; subs gave an immeasurably small speedup.
Found bug in calculation of insert penalties during mapping.  Fixing this bug increases speed but decreases accuracy, so it was modified toward a compromise.


v27.
Added command line to sam file header.
Added "msa=" flag.  You can specify which msa to use by entering the classname.
Added initial banded mode.  Specify "bandwidth=X" or "bandwidthratio=X" accelerate alignment.
Cleaned up argument parsing a bit.
Improved nodisk mode; now does not use the disk at all for indexing.  BBSplitter still uses the disk.
Added "fast" flag, which changes some paramters to make mapping go faster, with slightly lower sensitivity.
Improved error handling; corrupt input files should be more likely to crash with an error message and less likely to hang.  Noted by Alex Copeland.
Improved SAM input, particularly coordinates and cigar-string parsing; this should now be correct but requires an indexed reference.  Of course this information is irrelevant for mapping so this parsing is turned off by default for bbmap.
Increased maximum read speed with ByteFile2, by using 2 threads per file.  May be useful in input-speed limited scenarios, as when reading compressed input on a node with many cores.  Also accelerates sam input.
TODO:  Make ByteFile auto-select a subtype based on compression and number of cores.
TODO:  Consider moving THREADS to Shared.
Updated match/cigar flag syntax.
Updated shellscript documentation.
Changed ByteFile2 from array lists to arrays; should reduce overhead.
TODO:  Increase speed of sam input.
TODO:  Increase speed of output, for all formats.
TODO:  Finish ReadStreamWriter.addStringList(), which allows formatting to be done in the host.
In progress:  Moving all MapThread fields to abstract class.
MapThread now passes reverse-complemented bases to functions to prevent replication of this array.
Fixed very rare bug when a non-semiperfect site becomes semiperfect after realignment, but subsequently is no longer highest-ranked.
strictmaxindel can now be assigned a number (e.g. stricmaxindel=5).
If a fasta read is broken into pieces, now all pieces will recieve the _# suffix in their name.  Previously, the first piece was exempt.
TODO:  Consider changing SamLine.rname to a String and seq, qual to byte[].
Changed SamLine.seq, qual to byte[].  Now stored in original read order and only reversed for minus strand during I/O.
Added sortscaffolds flag (requested by Vasanth Singan).
Fixed XS tag bug; in some cases read 2 was getting opposite flag (noted by Vasanth Singan).
Fixed bug when reading sam files without qualities (noted by Brian Foster).
Fixed bug where absent cigar strings were printed as "null" instead of "*" as a result of recent changes to sam I/O (noted by Vasanth Singan).
Found error when a read goes off the beginning of a block.  Ref padding seems to be absent, because Ns were replaced by random sequence.  Cause is unknown; cannot replicate.
Fixed Block.getHitList(int, int).
Changed calcAffineScore() to require base array for information when throwing exceptions.
Changed generated bamscript to unload samtools module before loading samtools/0.1.19.
sam file idflag and stopflag are both now faster, particularly for perfect mappings.  But both default to off because they are still slow nonetheless.
Fixed bug in BBIndex in which a site was considered perfect because all bases matched the reference, but some of the bases were N.  Canonically, reads with Ns can never be perfect even if the ref has Ns in the same locations.
Fixed above bug again because it was not fully fixed:  CHECKSITES was allowing a read to be classified as perfect even if it contained an N.
Increased sam read speed by ~2x; 30MB/s to 66MB/s
Increased sam write speed from ~18MB/s to ~32MB/s on my 4-core computer (during mapping), with mapping at peak 42MB/s with out=null.  Standalone (no mapping) sam output seems to run at 51MB/s but it's hard to tell.
Increased fasta write from 118MB/s to 140 MB/s
Increased fastq write from 70MB/s to 100MB/s
Increased fastq read from 120MB/s (I think) to 296MB/s (663 megabytes/sec!) with 2 threads or 166MB/s with 1 thread
Some of these speed increases come from writing byte[] into char[] buffer held in a ThreadLocal, instead of turning them into Strings or appending them byte-by-byte.
All of these speed optimizations caused a few I/O bugs that temporarily affected some users between Oct 1 and Oct 4, 2013.  Sorry!
Flipped XS tag from + to - or vice versa.  I seem to have misinterpreted the Cufflinks documentation (noted by Vasanth Singan).
Fixed bug in which (as a result of speed optimizations) reads outside scaffold boundaries, in sam 1.3 format, were not getting clipped (Noted by Brian Foster).
Changed default behavior of all shellscripts to run with -Xmx4g if maximum memory cannot be detected (typically, because ulimit=infinity).  Was 31.  Unfortunately things will break either way.
Fixed off-by-1 error in sam TLEN calculation; also simplified it to give sign based on leftmost POS and always give a plus and minus even when POS is equal.
Added sam NH tag (when ambig=all).
Disabled sam XM tag because the bowtie documentation and output do not make any sense.
Changed sam MD and NM tags to account for 'N' symbol in cigar strings.
Made sam SM tag score compatible with mapping score.
Fixed bug in SamLine when cigar=f (null pointer when parsing match string). (Found by Vasanth Singan)
Fixed bug in BBMapThread* when local=true and ambiguous=toss (null pointer to read.list). (Found by Alexander Spunde)
Changed synthetic read naming and parsing (parsecustom flag) to use " /1" and " /2" at the end of paired read names. (Requested by Kurt LaButti)
Increased fastq write to 200MB/s (590 megabytes/s)
Increased fasta write to 212MB/s (624 megabytes/s measured by fastq input)
Increased sam write to 167MB/s (492 megabytes/s measured by fastq input)
Increased bread write to 196MB/s (579 megabytes/s measured by fastq input)
bf2 (multithreaded input) is now enabled by default on systems with >4 cores, or in ReformatReads always.
Fixed RTextOutputStream3.finishedSuccessfully() returning false when output was in 2 files.
Changed output streams to unbuffered.  No notable speed increase.
Fixed bug in ByteFile2 in which reads would be recycled when end of file was hit (found by Brian Foster, Bryce Foster, and Kecia Duffy).


v26.
Fixed crash from consecutive newlines in ByteFile.
Made SiteScore clonable/copyable.
Removed @RG line from headers.  It implies that reads should be annotated with addition fields based on the RG line information.
Changed sam flags (at advice of Joel Martin).  Now single-ended reads will never have flags 0x2, 0x40, or 0x80 set.
Added correct insert size average to output stats, in place of old inner distance and mapping length.
Fixed crash when detecting length of SamLines with no cigar string.  (Found by Shayna Stein)
Added flag "keepnames" which keeps the read names unchanged when writing in sam format.  Normally, a trailing "/1", "/2", " 1", or " 2" are stripped off, and if read 2's name differs from read 1's name, read 1's name is used for both.  This is to remain spec-compliant with the sam format.  However, in some cases (such as grading synthetic reads tagged with the correct mapping location) it is useful to retain the original name of each read.
Added local alignment option, "local".  Translates global alignments into a local alignments using the same affine transform (and soft-clips ends).
Changed killbadpairs default to false.  Now by default improperly paired reads are allowed.
Merged TranslateColorspaceRead versions into a single class.
Added interleaved input and output for bread format.  May be useful for error correction pipeline.
TODO:  Mode where reads are mapped to multiple scaffolds, but are mapped at most one time per scaffold.  I.e., remove all but top site per scaffold (and forbid self-mapping).
Fixed yet another instance of negative coordinates appearing in an unmapped read, which the new version of samtools can't handle.
Fixed bug in counting ambiguous reads; was improperly including in statistics reads that were ambiguous but had a score lower than minratio.
Fixed rare crash found related to realignment of reads with ambiguous mappings (found by Rob Egan).
Unified many of the differences between the MapThread variants, and added a new self-checking function (checkTopSite) to ensure a Read is self-consistent.
Added some bitflag fetch functions to SamLine and fixed 'pairedOnSameChrom()' which was not handling the '=' symbol.
TODO: Make GENERATE_BASE_SCORES_FROM_QUALITY a parameter, default false in BBMapPacBio and true elsewhere. (I verified this should work fine)
TODO: Make GENERATE_KEY_SCORES_FROM_QUALITY a parameter, default true (probably even in BBMapPacBio). (I verified this should work fine)
Updated LongM (merged with LongM from Dedupe).
Fixed bug in SamLine in which clipped leading indels were not considered, causing potential negative coordinates.  (Found by Brian Foster)
TODO: Match strings like NNNNNNDDDDDNNNNNmmmmmmmmmmmmmmmmm...mmmmmmm should never exist in the first place.  Why did that happen?
Added "strictmaxindel" flag (default: strictmaxindel=f).  Attempts to kill mappings in which there is a single indel event longer than the "maxindel" setting.  Requested by James Han.
TODO: Ensure strictmaxindel works in all situations, including rescued paired ends and recursively regenerated padded match strings.
TODO: Redo msa to be strictly subtractive.  Start with score=100*bases, then use e.g. 0 for match, -1 for del, -370 for sub, -100 for N, etc.  No need for negative values.
Changed TIMEBITS in MultiStateAligner9PacBio from 10 to 9 to address a score underflow assertion error found by Alicia Clum.  The underflow occuerd around length 5240; new limit should be around 10480.
TODO: Alicia found an error of exceeding gref bounds.
Fixed race condition in TextStreamWriter.
Improved functionality of splitter.  Now you can index once and map subsequently using "basename" without specifying "ref=" every single time.
"Reads Used" in output now dispays the number of reads used.  Before, for paired reads, it would display the number of pairs (half as many).
Added bases used to reads used at Kurt's request.
Improved bam script generation.  Now correctly sets samtools memory based on detected memory, and warns user that crashes may be memory-related.
Fixed an obsolete assertion in SamLine found by Alicia.
Added XS tag option ("xstag=t") for Cufflinks; the need for this was noted by requested by Vasanth Singan.
Added 'N' cigar operation for deletions longer than X bases (intronlen=X).  Also needed by Cufflinks.
Secondary alignments now get "*" for bases and qualities, as recommended by the SAM spec.  This saves space, but may cause problems when converting sam into other formats.
Fixed bug that caused interleaved=true to override in2.  Now if you set in and in2, interleaved input will be disabled.  (noted by Andrew Tritt).
Fixed some low-level bugs in I/O streams.  When shutting down streams I was waiting until !Thread.isAlive() rather than Thread.getState()==Thread.State.TERMINATED, which caused a race condition (since a thread is not alive before it starts execution).
Added debugging file with random name written to /ref/ directory.  This should help debugging if somewhere deep in a pipeline multiple processes try to index at the same location simultaneously.  Suggested by Bryce Foster.
Fixed log file generation causing a crash if the /ref/ directory did not exist, found by Vasanth Singan.  Also logging is now disabled by default but enabled if you set "log=t".
Input sequence data will now translate '.' and '-' to 'N' automatically, as some fasta databases appear to use '.' instead of 'N'.  (Thanks to Kecia Duffy and James Han)
Added capability to convert lowercase reads to upper case (crash on lowercase noted by Vasanth Singan).


v25.
Increased BBMapPacBio max read length to 6000, and BBMapPacBioSkimmer to 4000.
Fixed bugs in padding calculations during match string generation.
Improved some assertion error output.
Added flag "maxsites" for max alignments to print.
Added match field to sitescore.
Made untrim() affect sitescores as well.
Decreased read array buffer from 500 to 20 in MapPacBio.
TODO:  sam secondary alignments.  (PARTIALLY DONE)
TODO:  sam secondary match strings.
TODO:  stitcher for super long reads.
TODO:  msa superclass, "MSA".
TODO:  wrapper for split reference mapping and merging.
Improved fillAndScoreLimited to return additional information.
Added flag "secondary" to print secondary alignments.  Does not yet ensure that all secondary alignments will get cigar strings, but most do.
Added flag "quickmatch" to generate match strings for SiteScores during slow align.  Speeds up the overall process somewhat (at least on my PC; have not tested it on cluster).
Improved pruning during slow align by dynamically increasing msa limit.
Addressed a bug in which reads sometimes have additional sites aligned to the same coordinates as the primary site.  The bug can still occur (typically during match generation or as a result of padding), but is detected and corrected during runtime.
Tracked down and fixed a bug relating to negative coordinates in sam output for unmapped reads paired with reads mapped off the beginning of a scaffold, with help from Rob Egan.
Disabled frowny-face warning message which had caused some confusion.
TODO:  Add verification of match strings on site scores.
Made superclass for MSA.  This will allow merging of redundant code over the various BBMap versions.
Fixed a crash-hang out-of-memory error caused by initialization order.  Now crashes cleanly and terminates.  Found by James Han.
Fixed bug in output related to detecting cigar string length under sam 1.4 specification (found by Rob Egan).
Added flag "killbadpairs"/"kbp".
Added flag "fakequality" for fasta.
Permanently fixed bugs related to unexpected short match strings caused by error messages.
Increased speed of dynamic program phase when dealing with lots of Ns.
TODO:  In-line generation of short match string when printing a read, rather than mutating the read. (mutation is now temporary)
Added flag, "stoptag".  Allows generation of SAM tag YS:i:<read stop location>
Added flag, "idtag".  Allows generation of SAM tag YI:f:<percent identity>

v24.
Fixed bug that slightly reduced accuracy for reads with exactly 1 mismatch.  They were always skipping slow align, sometimes preventing ambiguous reads from being detected.
Increased speed of MakeRocCurve (for automatic grading of sam files from synthetic reads).  Had used 1 pass per quality level; now it uses only 1 pass total.
Increased accuracy of processing reads and contigs with ambiguous bases (in mapping phase).
Adjusted clearzones to use gradient functions and asymptotes rather than step functions.  Reduces false positives and increases true positives, especially near the old step cutoffs.
Fixed trimSitesBelowCutoff assertion that failed for paired reads.
Added single scaffold toggle to RandomReads.  Default 'singlescaffold=true'; forces reads to come from a single scaffold).  This can cause non-termination if no scaffolds are long enough, and may bias against shorter scaffolds.
Added min scaffold overlap to RandomReads.  Default 'overlap=1'; forces reads to overlap a scaffold at least this much.  This can cause non-termination if no scaffolds are long enough, and may bias against shorter scaffolds.
Fixed setPerfect().  Previously, reads with 'N' overlapping 'N' in the reference could be considered perfect matches, but no reads containing 'N' should ever be considered a perfect mapping to anything.
Formalized definition of semiperfect to require read having no ambiguous bases, and fixed "isSemiperfect()" function accordingly.
Shortened and clarified executable names.
Fixed soft-clipped read start position calculation (mainly relevant to grading).
Prevented reads from being double-counted when grading, when a program gives multiple primary alignments for a read.
Fixed a bug in splitter initialization.
Added "ambiguous2".  Reads that map to multiple references can now be written to distinct files (prefixed by "AMBIGUOUS_") or thrown away, independantly of whether they are ambiguous in the normal sense (which includes ambiguous within a single reference).
Added statistics tracking per reference and per scaffold.  Enable with "scafstats=<file>" or "refstats=<file>".
"ambiguous" may now be shortened to "ambig" on the command line.
"true" and "false" may now be shortened to t, 1, or f, 0.  If omitted entirely, "true" is assumed; e.g. "overwrite" is equivalent to "overwrite=true".
Added stderr as a vaild output destination specified from the command line.
BBSplitter now has a flag, "mapmode"; can be set to normal, accurate, pacbio, or pacbioskimmer.
Fixed issue where stuff was being written to stdout instead of stderr and ended up in SAM files (found by Brian Foster).
TODO: Add secondary alignments.
TODO: Unlimited length reads.
TODO: Protein mapping.
TODO: Soft clipping in both bbmap and GradeSamFile.  Should universally adjust coords by soft-clip amount when reported in SAM format.
Fixed assertion error concerning reads containing Ns marked as perfect, when aligned to reference Ns (found by Rob Egan).
Fixed potential null-pointer error in "showprogress" flag.

v23.
Created BBSplitter wrapper for BBMap that allows merging any number references together and splitting the output into different streams.
Added support for ambiguous=random with paired reads (before it was limited to unpaired).
TODO: Iterative anchored alignment for very long reads, with a full master gref.
TODO: untrim=c/m/s/n/r
TODO: mode=vfast/veryfast: k=14 minratio=0.8 minhits=2 maxindel=20
TODO: mode=fast: k=13 minratio=0.7 minhits=2 maxindel=200
TODO: mode=normal: k=13 minratio=0.56 minhits=1 maxindel=16000
TODO: mode=slow/accurate: BBMapi
TODO: mode=pacbio: BBMapPacBio k=12
TODO: mode=perfect
TODO: mode=semiperfect
TODO: mode=rnaseq
TODO: Put untrim in caclStatistics section
TODO: Test with MEGAN.
Finished new random read generator.  Much faster, and solves coordinate problem with multiple indels.
Improved error message on read parsing failures.
TODO: Insert size histogram
TODO: "outp=", output for reads that mapped paired
TODO: "outs=", output for reads that mapped singly
Corrected assertion in "isSingleScaffold()"
Fixed a rare bug preventing recursive realignment when ambiguous=random (found by Brian Foster)
Added samversion/samv flag.  Set to 1.3 for cigar strings with 'M' or 1.4 for cigar strings with '=' and 'X'.  Default is 1.3.
Added enforcement of thread limit when indexing.
Added internal autodetection of gpint machines.  Set default threadcount for gpints at 2.
Improved ability to map with maxindel=0
Added XM:i:<N> optional SAM flag because some programs seem to demand it.  Like all extra flags, this is omitted if the read is not mapped.  Otherwise, it is set to 1 for unambiguously mapped reads, and 2 or more for ambiguously mapped reads.  The number can range as high as the total number of equal-scoring sites, but this is not guaranteed unless the "ambiguous=random" flag is used.
Fixed bug in autodetection of paired ends, found by Rob Egan.



v22.
Added match histogram support.
Added quality histogram support.
Added interleaving support to random read generator.
Added ability to disable pair rescue ("rescue=false" flag), which can speed things up in some cases.
Disabled dynamic-programming slow alignment phase when no indels are allowed.
Accelerated rescue in perfect and semiperfect mode.
Vastly accelerated paired mapping against references with a very low expected mapping rate.
Fixed crash in rescue caused by reads without quality strings (e.g. paired fasta files). (found by Brian Foster)


v21.
If reference specified is same as already-processed reference, the old index will not be deleted.
Added BBMap memory usage estimator to assembly statistics tool:  java -Xmx120m jgi.AssemblyStats2 <fasta file> k=<kmer size for BBMap>
Added support for multiple output read streams: all reads (set by out=), mapped reads (set by outm=), and unmapped reads (set by outu=).  They can be in different formats and any combination can be used at once.  You can set pair output to secondary files with out2, outm2, and outu2.
Changed definition of "out=".  You can no longer specify split output streams implicitly by using a "#" in the filename; it must be explicit.  the "#" wildcard is still allowed for input streams.
Fixed a bug with sam input not working. (found by Brian Foster)
Added additional interleaved autodetection pattern for reads named "xxxxx 1:xxxx" and "xxxxx 2:xxxx"
Fixed a bug with soft-clipped deletions causing an incorrect cigar length. (found by Brian Foster)
Fixed a bug with parsing of negative numbers in byte arrays.
TODO: Found a new situation in which poly-N reads preferentially map to poly-N reference (probably tip search?)
Fixed a bug in which paired reads occasionally are incorrectly considered non-semiperfect. (found by Brian Foster)
Added more assertion tests for perfection/imperfection status.
Added blacklist support.  This allows selection of output stream based on the name of the scaffold to which a read maps.
Created Blacklist class, allowing creation of blacklists and whitelists.
Added outb (aka outblacklist) and outb2 streams, to output reads that mapped to blacklisted scaffolds. 
Added flag "outputblacklisted=<true/false>" which contols whether blacklisted reads are printed to the "out=" stream.  Default is true.
Added support for streaming references. e.g. "cat ref1.fa ref2.fa | java BBMap ref=stdin.fa"
Updated and reorganized this readme.
Removed a dependency on Java 7 libraries (so that the code runs in Java 6).
Added per-read error rate histogram.  Enable with qhist=<filename>
TODO: generate standard deviation.
Added per-base-position M/S/D/I/N rate tracking.  Enable with mhist=<filename>
Added quality trimming.  Reads may be trimmed prior to mapping, and optionally untrimmed after mapping, so that no data is lost.  Trimmed bases are reported as soft-clipped in this case.
Trimming will extend until at least 2 consecutive bases have a quality greater than trimq (default 5).
Added flags:  trim=<left/right/both/false>, trimq=<5>, untrim=<true/false>
TODO:  Correct insert size in realtime for trim length.
TODO:  Consider adding a TrimRead pointer to reads, rather than using obj.
TODO:  Consider extending match string as 'M' rather than 'C' as long as clipped bases match.
Found and made safe some instances where reads could be trimmed to less than kmer length.
Found and fixed instance where rescue was attempted for length-zero reads.
Fixed an instance where perfect reads were not marked perfect (while making match string).


v20.1 (not differentiated from v20 since the differences are minor)
Fixed a minor, longstanding bug that prevented minus-strand alignment of rads that only had a single valid key (due to low complexity or low quality).
Increased accuracy of perfectmode and semiperfectmode, by allowing mapping of reads with only one valid key, without loss of speed.  They still don't quite match normal mode since they use fewer keys.
Added detection of and error messages for reads that are too long to map. 
Improved shell script usage information.


v20.
Made all MapThreads subclasses of MapThread, eliminating duplicate code.
Any exception thrown by a MapThread will now be detected, allowing the process to complete normally without hanging.
Exceptions (e.g. OutOfMemory) when loading reference genome are now detected, typically causing a crash exit instead of a hang.
Exceptions (e.g. OutOfMemory) when generating index are now detected, causing a crash exit instead of a hang.
Exceptions in output stream (RTextOutputStream) subthreads are now detected, throwing an exception.
Added support for soft clipping.  All reads that go off the ends of scaffolds will be soft-clipped when output to SAM format. (The necessity of this was noted by Rob Egan, as negative scaffold indices can cause software such as samtools to crash)


v19.
Added support for leading FASTA comments (denoted by semicolon).
Fixed potential problem in FASTA read input stream with very long reads.
Recognizes additional FASTA file extensions: .seq, .fna, .ffn, .frn, .fsa, .fas
Disabled gzip subprocesses to circumvent a bug in UGE: Forking can cause a program to be terminated.  Gzip is still supported.
Slightly reduced memory allocation in shellscript.
Ported "Analyze Index" improvement over to all versions (except v5).
Added flags: fastaminread, showprogress
Fixed problem noted by Rob Egan in which paired-end reads containing mostly 'N' could be rescued by aligning to the poly-N section off the end of a contig.
Fixed: Synthetic read headers were being improperly parsed by new FASTQ input stream.
Made a new, faster, more correct version of "isSemiperfect".
Added "semiperfect" test for reads changed during findDeletions.
Identified locations in "scoreNoIndels" where call 'N' == ref 'N' is considered a match.  Does not seem to cause problems.
Noted that SAM flag 0x40 and 0x80 definitions differ from my usage.


v18.
Fastq read input speed doubled.
Fasta read input speed increased 50%.
Increased speed of "Analyze Index" by a factor of 3+ (just for BBMap so far; have not yet ported change over to other versions).
Fixed an array out-of-bounds bug found by Alicia Clum.
Added bam output option (relies on Samtools being installed).
Allows gzip subprocesses, which can sometimes improve gzipping and gunzipping speed over Java's implementation (will be used automatically if gzip is installed).  This can be disabled with with the flags "usegzip=false" and "usegunzip=false".
Started a 32-bit mode which allows 4GB per block instead of 2GB, for a slight memory savings (not finished yet).
Added nondeterministic random read sampling option.
Added flags: minscaf, startpad, stoppad, samplerate, sampleseed, kfilter, usegzip, usegunzip


v17.
Changed the way error rate statistics are displayed.  All now use match string length as denominator.
Identified error in random read generator regarding multiple insertions.  It will be hard to fix but does not matter much.
Found out-of-bounds error when filling gref.  Fixed (but maybe not everywhere...).
Added random mapping for ambiguous reads.
Changed index from 2d array to single array (saves a lot of memory).
Increased speed by ~10%.
Improved index generation and loading speed (typically more than doubled).
Changed chrom format to gzipped.
Added "nodisk" flag; index is not written to disk.
Fixed a rare out-of-bounds error.
Increased speed of perfect read mapping.
Fixed rare human PAR bug.


v16.  Changes since last version:
Supports unlimited number of unscaffolded contigs.
Supports piping in and out.  Set "out=stdout.sam" and "in=stdin.fq" to pipe in a fastq file and pipe out a sam file (other extensions are also supported).
Ambiguously named files (without proper extensions) will be autodetected as fasta or fastq (though I suggest not relying on that).
Added additional flags (described in parameters section): minapproxhits, padding, tipsearch, maxindel.
minapproxhits has a huge impact on speed.  Going from 1 to 2 will typically at least double the speed (on a large genome) at some cost to accuracy.


v15.  Changes since last version:
Contig names are retained for output.
SAM header @SQ tags fixed.
SAM header @PG tag added.
An out-of-bounds error was fixed.
An error related to short match strings was found and possibly handled.
All versions now give full statistics related to %matches, %substitutions, %deletions, and %insertions (unless match string generation is disabled).
Increased speed and accuracy for tiny (<20MB) genomes.
Added dynamic detection of scaffold sizes to better partition index, reducing memory in some cases.
Added command-line specification of kmer length.
Added more command line flags and described them in this readme.
Allowed overwriting of existing indices, for ease of use (only when overwrite=true).  For efficiency you should still only specify "ref=" the first time you map to a particular reference, and just specify the build number subsequently.
