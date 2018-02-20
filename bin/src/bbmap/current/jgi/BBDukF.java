package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLongArray;

import kmer.AbstractKmerTable;
import kmer.HashArray;
import kmer.HashForest;
import kmer.KmerTable;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.RTextOutputStream3;
import stream.Read;

import align2.ListNum;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;
import dna.AminoAcid;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextStreamWriter;

/**
 * Separates or trims reads based on matching kmers in a reference.
 * Supports arbitrary K and inexact matches. 
 * Supercedes BBDuk by replacing Java's HashMap with HashArray and HashForest.
 * @author Brian Bushnell
 * @date Aug 30, 2013
 *
 */
public class BBDukF {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		if(Tools.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		
		//Create a new BBDuk instance
		BBDukF bbd=new BBDukF(args);
		
		///And run it
		bbd.process();
	}
	
	/**
	 * Display usage information.
	 */
	private static void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("\njava -ea -Xmx20g -cp <path> jgi.BBDuk in=<input file> out=<output file> ref=<contaminant files>");
		outstream.println("\nOptional flags:");
		outstream.println("in=<file>          \tThe 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.");
		outstream.println("in2=<file>         \tUse this if 2nd read of pairs are in a different file.");
		outstream.println("out=<file>         \t(outnonmatch) The 'out=' flag is needed if the output file is not the second parameter.  'out=stdout' will pipe to standard out.");
		outstream.println("out2=<file>        \t(outnonmatch2) Use this to write 2nd read of pairs to a different file.");
		outstream.println("outmatch=<file>    \t(outm or outb) Write 'bad' reads here (containing contaminant kmers).");
		outstream.println("outmatch2=<file>   \t(outm2 or outb2) Use this to write 2nd read of pairs to a different file.");
		outstream.println("stats=<file>       \tWrite statistics about which contaminants were detected.");
		outstream.println("duk=<file>         \tWrite duk-like output.");
		outstream.println("");
		outstream.println("threads=auto       \t(t) Set number of threads to use; default is number of logical processors.");
		outstream.println("overwrite=t        \t(ow) Set to false to force the program to abort rather than overwrite an existing file.");
		outstream.println("showspeed=t        \t(ss) Set to 'f' to suppress display of processing speed.");
		outstream.println("interleaved=auto   \t(int) If true, forces fastq input to be paired and interleaved.");
		outstream.println("k=31               \tKmer length used for finding contaminants.  Contaminants shorter than k will not be found.");
		outstream.println("maskmiddle=t       \t(mm) Treat the middle base of a kmer as a wildcard.");
		outstream.println("maxbadkmers=0      \t(mbk) Reads with more than this many contaminant kmers will be discarded.");
		outstream.println("minavgquality=0    \t(maq) Reads with average quality (before trimming) below this will be discarded.");
		outstream.println("touppercase=f      \t(tuc) Change all letters in reads and reference to upper-case.");
		outstream.println("ktrim=f            \tTrim reads to remove bases matching reference kmers. ");
		outstream.println("                   \tValues: f (don't trim), r (trim right end), l (trim left end), n (convert to N instead of trimming).");
		outstream.println("useshortkmers=f    \t(usk) Look for shorter kmers at read tips (only for k-trimming).");
		outstream.println("mink=4             \tMinimum length of short kmers.  Setting this automatically sets useshortkmers=t.");
		outstream.println("qtrim=f            \tTrim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers. ");
		outstream.println("                   \tValues: t (trim both ends), f (neither end), r (right end only), l (left end only).");
		outstream.println("minq=4             \tTrim quality threshold.");
		outstream.println("minlength=2        \t(ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.");
		outstream.println("ziplevel=2         \t(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.");
		outstream.println("fastawrap=100      \tLength of lines in fasta output");
		outstream.println("qin=auto           \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
		outstream.println("qout=auto          \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
		outstream.println("rcomp=t            \tLook for reverse-complements of kmers also.");
		outstream.println("forest=t           \tUse HashForest data structure");
		outstream.println("table=f            \tUse KmerTable data structure");
		outstream.println("array=f            \tUse HashArray data structure");
	}
	
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public BBDukF(String[] args){
		for(String s : args){if(s.contains("standardout") || s.contains("stdout")){outstream=System.err;}}
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		FastaReadInputStream.SPLIT_READS=false;
		ByteFile.FORCE_MODE_BF2=Shared.THREADS>2;
		
		/* Initialize local variables with defaults */
		boolean setOut=false, setOutb=false, qtrimRight_=false, qtrimLeft_=false;
		boolean ktrimRight_=false, ktrimLeft_=false, ktrimN_=false, ktrimExclusive_=false;
		boolean addTrimmedToBad_=true;
		boolean rcomp_=true;
		boolean forbidNs_=false;
		boolean useForest_=false, useTable_=false, useArray_=true;
		int k_=28, kbig_=-1;
		int mink_=-1;
		int ways_=-1; //Currently disabled
		int maxBadKmers_=0;
		long skipreads_=0;
		byte qin=-1, qout=-1;
		byte TRIM_SYMBOL_='N';
		
		byte trimq_=4;
		byte minAvgQuality_=0;
		int minReadLength_=20;
		float minLenFraction_=0f;
		int maxNs_=-1;
		boolean removePairsIfEitherBad_=true;
		
		scaffoldNames.add(""); //Necessary so that the first real scaffold gets an id of 1, not zero
		
		{
			boolean b=false;
			assert(b=true);
			EA=b;
		}
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("out") || a.equals("out1") || a.equals("outu") || a.equals("outu1") || a.equals("outnonmatch") || 
					a.equals("outnonmatch1") || a.equals("outunnmatch") || a.equals("outunmatch1") || a.equals("outunnmatched") || a.equals("outunmatched1")){
				out1=b;
				setOut=true;
			}else if(a.equals("out2") || a.equals("outu2") || a.equals("outnonmatch2") || a.equals("outunmatch2") || 
					a.equals("outnonmatched2") || a.equals("outunmatched2")){
				out2=b;
			}else if(a.equals("outb") || a.equals("outm") || a.equals("outb1") || a.equals("outm1") || a.equals("outbad") || 
					a.equals("outbad1") || a.equals("outmatch") || a.equals("outmatch1")){
				outb1=b;
				setOut=true;
			}else if(a.equals("outb2") || a.equals("outm2") || a.equals("outbad2") || a.equals("outmatch2")){
				outb2=b;
			}else if(a.equals("stats")){
				outstats=b;
			}else if(a.equals("duk") || a.equals("outduk")){
				outduk=b;
			}else if(a.equals("rqc")){
				outrqc=b;
			}else if(a.equals("ref")){
				ref=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
			}else if(a.equals("literal")){
				literal=(b==null) ? null : b.split(",");
//				assert(false) : b+", "+Arrays.toString(literal);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("forest")){
				useForest_=Tools.parseBoolean(b);
				if(useForest_){useTable_=useArray_=false;}
			}else if(a.equals("table")){
				useTable_=Tools.parseBoolean(b);
				if(useTable_){useForest_=useArray_=false;}
			}else if(a.equals("array")){
				useArray_=Tools.parseBoolean(b);
				if(useArray_){useTable_=useForest_=false;}
			}else if(a.equals("ways")){
				ways_=Integer.parseInt(b);
			}else if(a.equals("bf1")){
				ByteFile.FORCE_MODE_BF1=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF2=!ByteFile.FORCE_MODE_BF1;
			}else if(a.equals("bf2")){
				ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF1=!ByteFile.FORCE_MODE_BF2;
			}else if(a.equals("usegzip") || a.equals("gzip")){
				ReadWrite.USE_GZIP=Tools.parseBoolean(b);
			}else if(a.equals("usepigz") || a.equals("pigz")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					int zt=Integer.parseInt(b);
					if(zt<1){ReadWrite.USE_PIGZ=false;}
					else{
						ReadWrite.USE_PIGZ=true;
						if(zt>1){
							ReadWrite.MAX_ZIP_THREADS=zt;
							ReadWrite.ZIP_THREAD_DIVISOR=1;
						}
					}
				}else{ReadWrite.USE_PIGZ=Tools.parseBoolean(b);}
			}else if(a.equals("usegunzip") || a.equals("gunzip")){
				ReadWrite.USE_GUNZIP=Tools.parseBoolean(b);
			}else if(a.equals("useunpigz") || a.equals("unpigz")){
				ReadWrite.USE_UNPIGZ=Tools.parseBoolean(b);
			}else if(a.equals("interleaved") || a.equals("int")){
				if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
				else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}else if(a.equals("tuc") || a.equals("touppercase")){
				Read.TO_UPPER_CASE=Tools.parseBoolean(b);
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.equals("k")){
				assert(b!=null) : "\nThe k key needs an integer value greater than 0, such as k=28\n";
				k_=Integer.parseInt(b);
				if(k_>31){
					kbig_=k_;
					k_=31;
				}else{
					kbig_=-1;
				}
				assert(k_>0 && k_<32) : "k must be at least 1; default is 28.";
			}else if(a.equals("mink") || a.equals("kmin")){
				mink_=Integer.parseInt(b);
				assert(mink_<0 || (mink_>0 && mink_<32)) : "kmin must be between 1 and 31; default is 4, negative numbers disable it.";
			}else if(a.equals("useshortkmers") || a.equals("shortkmers") || a.equals("usk")){
				useShortKmers=Tools.parseBoolean(b);
			}else if(a.equals("dist") || a.equals("distance") || a.equals("hdist") || a.equals("hammingdistance")){
				hammingDistance=Integer.parseInt(b);
				assert(hammingDistance>=0 && hammingDistance<4) : "hamming distance must be between 0 and 3; default is 0.";
			}else if(a.equals("edits") || a.equals("edist") || a.equals("editdistance")){
				editDistance=Integer.parseInt(b);
				assert(editDistance>=0 && editDistance<3) : "edit distance must be between 0 and 2; default is 0.";
			}else if(a.equals("maxskip") || a.equals("mxs")){
				maxSkip=Integer.parseInt(b);
			}else if(a.equals("minskip") || a.equals("mns")){
				minSkip=Integer.parseInt(b);
			}else if(a.equals("skipreads")){
				skipreads_=Long.parseLong(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.THREADS : Integer.parseInt(b));
			}else if(a.equals("maxbadkmers") || a.equals("mbk")){
				maxBadKmers_=Integer.parseInt(b);
			}else if(a.equals("minavgquality") || a.equals("maq")){
				minAvgQuality_=(byte)Integer.parseInt(b);
			}else if(a.equals("maxns")){
				maxNs_=Integer.parseInt(b);
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
				assert(WAYS>1) : "WAYS=1 is for debug mode.";
//				verbose=Tools.parseBoolean(b); //123
				if(verbose){outstream=System.err;} //For some reason System.out does not print in verbose mode.
			}else if(a.equals("mm") || a.equals("maskmiddle")){
				maskMiddle=Tools.parseBoolean(b);
			}else if(a.equals("rcomp")){
				rcomp_=Tools.parseBoolean(b);
			}else if(a.equals("forbidns") || a.equals("forbidn") || a.equals("fn")){
				forbidNs_=Tools.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("removeifeitherbad") || a.equals("rieb")){
				removePairsIfEitherBad_=Tools.parseBoolean(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
			}else if(a.equals("fastaminlen") || a.equals("fastaminlength")){
				FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(a.equals("ktrim")){
				if(b==null){b="";}
				if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){ktrimLeft_=true;ktrimRight_=false;ktrimN_=false;}
				else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){ktrimLeft_=false;ktrimRight_=true;ktrimN_=false;}
				else if(b.equalsIgnoreCase("n")){ktrimLeft_=false;ktrimRight_=false;ktrimN_=true;}
				else if(b.length()==1 && !b.equalsIgnoreCase("t") && !b.equalsIgnoreCase("f")){
					ktrimLeft_=false;ktrimRight_=false;ktrimN_=true;
					TRIM_SYMBOL_=(byte)b.charAt(0);
				}else{
					boolean x=Tools.parseBoolean(b);
					assert(!x) : "\nInvalid setting for ktrim - values must be f (false), l (left), r (right), or n.";
					ktrimRight_=ktrimLeft_=ktrimN_=x;
				}
			}else if(a.equals("ktrimright")){
				ktrimRight_=Tools.parseBoolean(b);
				ktrimLeft_=ktrimN_=!(ktrimRight_);
			}else if(a.equals("ktrimleft")){
				ktrimLeft_=Tools.parseBoolean(b);
				ktrimRight_=ktrimN_=!(ktrimLeft_);
			}else if(a.equals("ktrimn")){
				ktrimN_=Tools.parseBoolean(b);
				ktrimLeft_=ktrimRight_=!(ktrimN_);
			}else if(a.equals("ktrimexclusive")){
				ktrimExclusive_=Tools.parseBoolean(b);
			}else if(a.equals("trim") || a.equals("qtrim")){
				if(b==null){qtrimRight_=qtrimLeft_=true;}
				else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){qtrimLeft_=true;qtrimRight_=false;}
				else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){qtrimLeft_=false;qtrimRight_=true;}
				else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){qtrimLeft_=qtrimRight_=true;}
				else if(Character.isDigit(b.charAt(0))){
					if(!qtrimLeft_ && !qtrimRight_){qtrimLeft_=qtrimRight_=true;}
					trimq_=Byte.parseByte(b);
				}else{qtrimRight_=qtrimLeft_=Tools.parseBoolean(b);}
			}else if(a.equals("optitrim") || a.equals("otf") || a.equals("otm")){
				if(b!=null && (b.charAt(0)=='.' || Character.isDigit(b.charAt(0)))){
					TrimRead.optimalMode=true;
					TrimRead.optimalBias=Float.parseFloat(b);
					assert(TrimRead.optimalBias>=0 && TrimRead.optimalBias<1);
				}else{
					TrimRead.optimalMode=Tools.parseBoolean(b);
				}
			}else if(a.equals("trimright") || a.equals("qtrimright")){
				qtrimRight_=Tools.parseBoolean(b);
			}else if(a.equals("trimleft") || a.equals("qtrimleft")){
				qtrimLeft_=Tools.parseBoolean(b);
			}else if(a.equals("trimq") || a.equals("trimquality")){
				trimq_=Byte.parseByte(b);
			}else if(a.equals("q102matrix") || a.equals("q102m")){
				CalcTrueQuality.q102matrix=b;
			}else if(a.equals("qbpmatrix") || a.equals("bqpm")){
				CalcTrueQuality.qbpmatrix=b;
			}else if(a.equals("adjustquality") || a.equals("adjq")){
				TrimRead.ADJUST_QUALITY=Tools.parseBoolean(b);
			}else if(a.equals("otm") || a.equals("outputtrimmedtomatch")){
				addTrimmedToBad_=Tools.parseBoolean(b);
			}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
				minReadLength_=Integer.parseInt(b);
			}else if(a.equals("mlf") || a.equals("minlenfrac") || a.equals("minlenfraction") || a.equals("minlengthfraction")){
				minLenFraction_=Float.parseFloat(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
			}else if(a.equals("ascii") || a.equals("quality") || a.equals("qual")){
				byte x;
				if(b.equalsIgnoreCase("sanger")){x=33;}
				else if(b.equalsIgnoreCase("illumina")){x=64;}
				else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;}
				else{x=(byte)Integer.parseInt(b);}
				qin=qout=x;
			}else if(a.equals("asciiin") || a.equals("qualityin") || a.equals("qualin") || a.equals("qin")){
				byte x;
				if(b.equalsIgnoreCase("sanger")){x=33;}
				else if(b.equalsIgnoreCase("illumina")){x=64;}
				else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=true;}
				else{x=(byte)Integer.parseInt(b);}
				qin=x;
			}else if(a.equals("asciiout") || a.equals("qualityout") || a.equals("qualout") || a.equals("qout")){
				byte x;
				if(b.equalsIgnoreCase("sanger")){x=33;}
				else if(b.equalsIgnoreCase("illumina")){x=64;}
				else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY_OUT=true;}
				else{x=(byte)Integer.parseInt(b);}
				qout=x;
			}else if(a.equals("qauto")){
				FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;
			}else if(i==0 && in1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				in1=args[i];
			}else if(i==1 && out1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				out1=args[i];
				setOut=true;
			}else if(i==2 && ref==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				ref=(new File(args[i]).exists() ? new String[] {args[i]} : args[i].split(","));
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
		/* Set final variables; post-process and validate argument combinations */
		
		useForest=useForest_;
		useTable=useTable_;
		useArray=useArray_;
		hammingDistance=Tools.max(editDistance, hammingDistance);
		minSkip=Tools.max(1, Tools.min(minSkip, maxSkip));
		maxSkip=Tools.max(minSkip, maxSkip);
		addTrimmedToBad=addTrimmedToBad_;
		rcomp=rcomp_;
		forbidNs=(forbidNs_ || hammingDistance<1);
		trimSymbol=TRIM_SYMBOL_;
		skipreads=skipreads_;
		trimq=trimq_;
		minAvgQuality=minAvgQuality_;
		minReadLength=minReadLength_;
		minLenFraction=minLenFraction_;
		removePairsIfEitherBad=removePairsIfEitherBad_;
		maxNs=maxNs_;
//		if(maxReads>0){ReadWrite.USE_GUNZIP=ReadWrite.USE_UNPIGZ=false;}
		
//		if(ways_==-1){
//			while(Shared.THREADS<WAYS){WAYS-=2;}
//		}else{
//			WAYS=ways_;
//		}
//		WAYS=Tools.max(WAYS, 1);
		
		if(ktrimLeft_ || ktrimRight_ || ktrimN_){
			if(kbig_>k_){
				System.err.println("***********************   WARNING   ***********************"); 
				System.err.println("WARNING: When kmer-trimming, the maximum value of K is "+k_+".");
				System.err.println("K has been reduced from kbig_ to "+k_+".");
				System.err.println("***********************************************************"); 
				kbig_=k_;
			}
		}
		
		k=k_;
		k2=k-1;
		kbig=kbig_;
		if(kbig>k){
			minSkip=maxSkip=0;
			if(maskMiddle){
				System.err.println("maskMiddle was disabled because kbig>k");
				maskMiddle=false;
			}
		}
		mink=Tools.min((mink_<1 ? 4 : mink_), k);
		maxBadKmers=maxBadKmers_;
		if(mink_>0 && mink_<k){useShortKmers=true;}
		if(useShortKmers){
			if(maskMiddle){
				System.err.println("maskMiddle was disabled because useShortKmers=true");
				maskMiddle=false;
			}
		}
//		kmask=useShortKmers ? 1L<<(2*k) : 0;
//		assert(useShortKmers==kmask>0);
		
		ktrimRight=ktrimRight_;
		ktrimLeft=ktrimLeft_;
		ktrimN=ktrimN_;
		ktrimExclusive=ktrimExclusive_;
		
		assert(!useShortKmers || ktrimRight || ktrimLeft || ktrimN) : "\nSetting mink or useShortKmers also requires setting a ktrim mode, such as 'r', 'l', or 'n'\n";
		
		qtrimRight=qtrimRight_;
		qtrimLeft=qtrimLeft_;
		
		middleMask=maskMiddle ? ~(3L<<(2*(k/2))) : -1L;
		
		hitCounts=(outduk==null ? null : new AtomicLongArray(1001));
		
		
		/* Adjust I/O settings and filenames */
		
		if(qin!=-1 && qout!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY=false;
		}else if(qin!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.DETECT_QUALITY=false;
		}else if(qout!=-1){
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY_OUT=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(in1!=null && in1.contains("#") && !new File(in1).exists()){
			int pound=in1.lastIndexOf('#');
			String a=in1.substring(0, pound);
			String b=in1.substring(pound+1);
			in1=a+1+b;
			in2=a+2+b;
		}
		if(in2!=null && (in2.contains("=") || in2.equalsIgnoreCase("null"))){in2=null;}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(out1!=null && out1.contains("#")){
			int pound=out1.lastIndexOf('#');
			String a=out1.substring(0, pound);
			String b=out1.substring(pound+1);
			out1=a+1+b;
			out2=a+2+b;
		}
		
		if(outb1!=null && outb1.contains("#")){
			int pound=outb1.lastIndexOf('#');
			String a=outb1.substring(0, pound);
			String b=outb1.substring(pound+1);
			outb1=a+1+b;
			outb2=a+2+b;
		}

		if(!setOut){
			out1="stdout.fq";
			outstream=System.err;
			out2=null;
		}else if("stdout".equalsIgnoreCase(out1) || "standarddout".equalsIgnoreCase(out1)){
			out1="stdout.fq";
			outstream=System.err;
			out2=null;
		}
		if(out1!=null && !Tools.canWrite(out1, overwrite)){throw new RuntimeException("Output file "+out1+" already exists, and overwrite="+overwrite);}

		assert(!in1.equalsIgnoreCase(out1));
		assert(!in1.equalsIgnoreCase(outb1));
		assert(!in1.equalsIgnoreCase(in2));
		assert(out1==null || !out1.equalsIgnoreCase(out2));
		assert(out1==null || !out1.equalsIgnoreCase(outb1));
		assert(outb1==null || !outb1.equalsIgnoreCase(outb2));
		assert(THREADS>0);

		assert(in1==null || in1.toLowerCase().startsWith("stdin") || in1.toLowerCase().startsWith("standardin") || new File(in1).exists()) : "Can't find "+in1;
		assert(in2==null || in2.toLowerCase().startsWith("stdin") || in2.toLowerCase().startsWith("standardin") || new File(in2).exists()) : "Can't find "+in2;
		assert((ref!=null || literal!=null) || qtrimLeft || qtrimRight || minAvgQuality>0 || maxNs>=0) : "No reference files specified, no trimming set, no min avg quality - nothing to do.";
		if(ref!=null){
			for(String s0 : ref){
				assert(s0!=null) : "Specified a null reference.";
				String s=s0.toLowerCase();
				assert(s==null || s.startsWith("stdin") || s.startsWith("standardin") || new File(s0).exists()) : "Can't find "+s0;
			}
		}
		
		//Initialize tables
		for(int i=0; i<keySets.length; i++){
			if(useForest){
				keySets[i]=new HashForest(initialSize, true);
			}else if(useTable){
				keySets[i]=new KmerTable(initialSize, true);
			}else if(useArray){
				keySets[i]=new HashArray(initialSize, true);
			}else{
				throw new RuntimeException("Must use forest, table, or array data structure.");
			}
		}
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void process(){
		
		/* Check for output file collisions */
		Tools.testOutputFiles(overwrite, false, out1, out2, outb1, outb2, outstats, outduk, outrqc);
		
		/* Start overall timer */
		Timer t=new Timer();
		t.start();
		
//		boolean dq0=FASTQ.DETECT_QUALITY;
//		boolean ti0=FASTQ.TEST_INTERLEAVED;
//		int rbl0=Shared.READ_BUFFER_LENGTH;
//		FASTQ.DETECT_QUALITY=false;
//		FASTQ.TEST_INTERLEAVED=false;
//		Shared.READ_BUFFER_LENGTH=16;
		
		process2(t.time1);
		
//		FASTQ.DETECT_QUALITY=dq0;
//		FASTQ.TEST_INTERLEAVED=ti0;
//		Shared.READ_BUFFER_LENGTH=rbl0;
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		
		
		if(showSpeed){
			double rpnano=readsIn/(double)(t.elapsed);
			double bpnano=basesIn/(double)(t.elapsed);

			//Format with k or m suffixes
			String rpstring=(readsIn<100000 ? ""+readsIn : readsIn<100000000 ? (readsIn/1000)+"k" : (readsIn/1000000)+"m");
			String bpstring=(basesIn<100000 ? ""+basesIn : basesIn<100000000 ? (basesIn/1000)+"k" : (basesIn/1000000)+"m");

			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}

			outstream.println("Time:   \t\t\t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException("BBDuk terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	public void process2(long startTime){
		
		/* Start phase timer */
		Timer t=new Timer();
		t.start();

		if(DISPLAY_PROGRESS){
			outstream.println("Initial:");
			printMemory();
			outstream.println();
		}
		
		/* Fill tables with reference kmers */
		fillSet_MT();
		
//		if(THREADS>=4){
//			fillSet_MT();
//		}else{
//			fillSet_ST(); //Still in BBDukF_old
//		}
		
		/* Do kmer matching of input reads */
		findKmers(t);
		
		/* Write statistics to files */
		writeStats(t);
		writeDuk(System.nanoTime()-startTime);
		writeRqc();
		
		outstream.println("Input:                  \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		
		if(ref!=null || literal!=null){
			outstream.println("Contaminants:           \t"+readsKFiltered+" reads ("+String.format("%.2f",readsKFiltered*100.0/readsIn)+"%) \t"+
					basesKFiltered+" bases ("+String.format("%.2f",basesKFiltered*100.0/basesIn)+"%)");
			outstream.flush();
		}
		if(qtrimLeft || qtrimRight){
			outstream.println("QTrimmed:               \t"+readsQTrimmed+" reads ("+String.format("%.2f",readsQTrimmed*100.0/readsIn)+"%) \t"+
					basesQTrimmed+" bases ("+String.format("%.2f",basesQTrimmed*100.0/basesIn)+"%)");
		}
		if(ktrimLeft || ktrimRight || ktrimN){
			outstream.println("KTrimmed:               \t"+readsKTrimmed+" reads ("+String.format("%.2f",readsKTrimmed*100.0/readsIn)+"%) \t"+
					basesKTrimmed+" bases ("+String.format("%.2f",basesKTrimmed*100.0/basesIn)+"%)");
		}
		if(minAvgQuality>0 || maxNs>=0){
			outstream.println("Low quality discards:   \t"+readsQFiltered+" reads ("+String.format("%.2f",readsQFiltered*100.0/readsIn)+"%) \t"+
					basesQFiltered+" bases ("+String.format("%.2f",basesQFiltered*100.0/basesIn)+"%)");
		}
		outstream.println("Result:                 \t"+readsOut+" reads ("+String.format("%.2f",readsOut*100.0/readsIn)+"%) \t"+
				basesOut+" bases ("+String.format("%.2f",basesOut*100.0/basesIn)+"%)");
	}
	
	/**
	 * Write statistics about how many reads matched each reference scaffold.
	 * @param t Phase timer
	 */
	private void writeStats(Timer t){
		if(outstats==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outstats, overwrite, false, false);
		tsw.start();
		
		/* Create StringNum list of scaffold names and hitcounts */
		ArrayList<StringNum> list=new ArrayList<StringNum>();
		for(int i=1; i<scaffoldNames.size(); i++){
			final int num=scaffoldCounts.get(i);
			if(num>0){
				final String s=scaffoldNames.get(i);
				final StringNum sn=new StringNum(s, num);
				list.add(sn);
			}
		}
		Collections.sort(list);
		final double rmult=100.0/(readsIn>0 ? readsIn : 1);
		for(int i=0; i<list.size(); i++){
			StringNum sn=list.get(i);
			double pct=sn.n*rmult;
			tsw.print(sn.toString()+String.format("\t%.3f%%\n",pct));
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Write processing statistics in DUK's format.
	 * @param time Elapsed time, nanoseconds
	 */
	private void writeDuk(long time){
		if(outduk==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outduk, overwrite, false, false);
		tsw.start();
		tsw.println(dukString(time));
		tsw.poisonAndWait();
	}
	
	/**
	 * Write RQCFilter stats.
	 * @param time Elapsed time, nanoseconds
	 */
	private void writeRqc(){
		if(outrqc==null){return;}
		addToRqcMap();
		if(outrqc.endsWith("hashmap")){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outrqc, overwrite, false, false);
		tsw.start();
		tsw.println(rqcString());
		tsw.poisonAndWait();
	}
	
	public static String rqcString(){
		if(RQC_MAP==null){return null;}
		StringBuilder sb=new StringBuilder();
		
		String[] keys=new String[] {"inputReads", "inputBases", "qtrimmedReads", "qtrimmedBases", "qfilteredReads", "qfilteredBases",
				"ktrimmedReads", "ktrimmedBases", "kfilteredReads", "kfilteredBases", "outputReads", "outputBases"};
		
		for(String key : keys){
			String value=RQC_MAP.get(key);
			if(value!=null){
				sb.append(key+"="+value+"\n");
			}
		}
		
		return sb.toString();
	}
	
	private void addToRqcMap(){
		putRqc("inputReads", readsIn, false);
		putRqc("inputBases", basesIn, false);
		if(qtrimLeft || qtrimRight){
			putRqc("qtrimmedReads", readsQTrimmed, false);
			putRqc("qtrimmedBases", basesQTrimmed, false);
		}
		putRqc("qfilteredReads", readsQFiltered, false);
		putRqc("qfilteredBases", basesQFiltered, false);
		
		if(ktrimLeft || ktrimRight || ktrimN){
			putRqc("ktrimmedReads", readsKTrimmed, true);
			putRqc("ktrimmedBases", basesKTrimmed, true);
		}else{
			putRqc("kfilteredReads", readsKFiltered, false);
			putRqc("kfilteredBases", basesKFiltered, false);
		}
		putRqc("outputReads", readsOut, true);
		putRqc("outputBases", basesOut, true);
	}
	
	private static void putRqc(String key, long value, boolean evict){putRqc(key, value+"", evict);}
	
	private static void putRqc(String key, String value, boolean evict){
		if(RQC_MAP==null){RQC_MAP=new HashMap<String,String>();}
		if(evict || !RQC_MAP.containsKey(key)){RQC_MAP.put(key, value);}
	}
	
	/**
	 * Helper method; formats statistics to be duk-compatible
	 * @param time Elapsed time, nanoseconds
	 * @return duk output string
	 */
	private String dukString(long time){
		StringBuilder sb=new StringBuilder();
		sb.append("##INPUT PARAMETERS##\n");
		sb.append("#Reference file:	"+(ref==null || ref.length<1 ? null : ref.length==1 ? ref[0] : Arrays.toString(ref))+"\n");
		sb.append("#Query file:	"+in1+(in2==null ? "" : ","+in2)+"\n");
		sb.append("#Not matched reads file:	"+out1+(out2==null ? "" : ","+out2)+"\n");
		sb.append("#Matched reads file:	"+outb1+(outb2==null ? "" : ","+outb2)+"\n");
		sb.append("#Output file (duk):	"+outduk+"\n");
		sb.append("#Output file (stats):	"+outstats+"\n");
		sb.append("#Mer size:	"+k+"\n");
		long size=0;
		for(AbstractKmerTable x : keySets){size+=x.size();}
		sb.append("#Avg step size:	"+String.format("%.1f", refKmers/(double)(Tools.max(1, size)))+"\n");
		sb.append("#Cut off:	"+maxBadKmers+"\n");
		sb.append("#Mask middle:	"+maskMiddle+"\n");
		sb.append("#Quality trim:	"+((qtrimLeft || qtrimRight) ? trimq : "false")+"\n");
		sb.append("\n");
		
		sb.append("##REFERECE STAT##\n");
		sb.append("#Total Reads:	"+refReads+"\n");
		sb.append("#Total Bases:	"+refBases+"\n");
		sb.append("#Total kmers:	"+refKmers+"\n");
		sb.append("#Total stored kmers:	"+size+"\n");
		sb.append("\n");

		sb.append("## ELAPSED TIME##\n");
		sb.append("# Time:	"+String.format("%.2f", time/1000000000.0)+" seconds\n");
		sb.append("\n");

		sb.append("##QUERY FILE STAT##\n");
		sb.append("# Total number of reads:	"+readsIn+"\n");
		sb.append("# Total number of matched reads:	"+readsKFiltered+"\n");
		sb.append("# Match ratio:	"+String.format("%.6f", readsKFiltered*1.0/readsIn)+"\n");
		sb.append("\n");

		sb.append("##P-VALUE##\n");
		sb.append("#Avg number of Kmer for each read:	"+((basesIn/(Tools.max(readsIn, 1)))-k)+"\n");
//		sb.append("# P value for the given threshold 1 is 4.05231e-14\n"); //duk prints a P value; not sure what it means
		sb.append("\n");

		sb.append("## Histogram of kmer occurance for reads with at least one occurance ##\n");
		sb.append("#NumOcc\tNumReads\tPercentage\n");
		
		long sum=0;
		for(int i=0; i<hitCounts.length(); i++){sum+=hitCounts.get(i);}
		double mult=100.0/(sum<1 ? 1 : sum);
		for(int i=0; i<hitCounts.length(); i++){
			long x=hitCounts.get(i);
			if(x>0){
				sb.append(i).append('\t').append(x).append('\t').append(String.format("%.4f",(x*mult))).append('\n');
			}
		}
		
		return sb.toString();
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	

	/**
	 * Fills tables with kmers from references, using multiple LoadThread.
	 * @return Number of kmers stored.
	 */
	private long fillSet_MT(){
		Timer t=new Timer();
		t.start();
		if((ref==null || ref.length<1) && (literal==null || literal.length<1)){return 0;}
		long added=0;
		
		/* Create load threads */
		LoadThread[] loaders=new LoadThread[WAYS];
		for(int i=0; i<loaders.length; i++){
			loaders[i]=new LoadThread(i);
			loaders[i].start();
		}

		/* For each reference file... */
		if(ref!=null){
			for(String refname : ref){

				/* Start an input stream */
				FileFormat ff=FileFormat.testInput(refname, FileFormat.FASTA, null, false, true);
				ConcurrentReadStreamInterface cris=ConcurrentGenericReadInputStream.getReadInputStream(-1L, false, false, ff, null);
				Thread cristhread=new Thread(cris);
				cristhread.start();
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);

				/* Iterate through read lists from the input stream */
				while(reads!=null && reads.size()>0){
					{
						/* Assign a unique ID number to each scaffold */
						ArrayList<Read> reads2=new ArrayList<Read>(reads);
						for(Read r : reads2){
							final Integer id=scaffoldNames.size();
							scaffoldNames.add(r.id==null ? id.toString() : r.id);
							r.obj=id;
							if(r.mate!=null){r.mate.obj=id;}
						}

						/* Send a pointer to the read list to each LoadThread */
						for(LoadThread lt : loaders){
							boolean b=true;
							while(b){
								try {
									lt.queue.put(reads2);
									b=false;
								} catch (InterruptedException e) {
									//TODO:  This will hang due to still-running threads.
									throw new RuntimeException(e);
								}
							}
						}
					}

					/* Dispose of the old list and fetch a new one */
					cris.returnList(ln, ln.list.isEmpty());
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				/* Cleanup */
				cris.returnList(ln, ln.list.isEmpty());
				errorState|=ReadWrite.closeStream(cris);
			}
		}

		/* If there are literal sequences to use as references */
		if(literal!=null){
			ArrayList<Read> list=new ArrayList<Read>(literal.length);
			if(verbose){System.err.println("Adding literals "+Arrays.toString(literal));}

			/* Assign a unique ID number to each literal sequence */
			for(int i=0; i<literal.length; i++){
				final Integer id=scaffoldNames.size();
				Read r=new Read(literal[i].getBytes(), null, id);
				scaffoldNames.add(id.toString());
				r.obj=id;
				list.add(r);
			}

			/* Send a pointer to the read list to each LoadThread */
			for(LoadThread lt : loaders){
				boolean b=true;
				while(b){
					try {
						lt.queue.put(list);
						b=false;
					} catch (InterruptedException e) {
						//TODO:  This will hang due to still-running threads.
						throw new RuntimeException(e);
					}
				}
			}
		}
		
		/* Signal loaders to terminate */
		for(LoadThread lt : loaders){
			boolean b=true;
			while(b){
				try {
					lt.queue.put(POISON);
					b=false;
				} catch (InterruptedException e) {
					//TODO:  This will hang due to still-running threads.
					throw new RuntimeException(e);
				}
			}
		}
		
		/* Wait for loaders to die, and gather statistics */
		for(LoadThread lt : loaders){
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					lt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			added+=lt.added;
			refKmers+=lt.refKmersT;
			refBases+=lt.refBasesT;
			refReads+=lt.refReadsT;
		}
		//Correct statistics for number of threads, since each thread processes all reference data
		refKmers/=WAYS;
		refBases/=WAYS;
		refReads/=WAYS;
		
		scaffoldCounts=new AtomicIntegerArray(scaffoldNames.size());

		t.stop();
		if(DISPLAY_PROGRESS){
			outstream.println("Added "+added+" kmers; time: \t"+t);
			printMemory();
			outstream.println();
		}
		
		if(verbose){
			TextStreamWriter tsw=new TextStreamWriter("stdout", false, false, false, FileFormat.TEXT);
			tsw.start();
			for(AbstractKmerTable table : keySets){
				table.dumpKmersAsText(tsw, k);
			}
			tsw.poisonAndWait();
		}
		
		return added;
	}
	
	/**
	 * Match reads against reference kmers, using multiple ProcessThread.
	 * @param t
	 */
	private void findKmers(Timer t){
		
		/* Create read input stream */
		final ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, false, ff1, ff2);
			Thread cristhread=new Thread(cris);
			cristhread.start();
		}
		
		/* Create read output streams */
		RTextOutputStream3 ros=null, rosb=null;
		if(out1!=null){
			final int buff=(!ORDERED ? 12 : Tools.max(32, 2*Shared.THREADS));
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, ORDERED);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, ORDERED);
			ros=new RTextOutputStream3(ff1, ff2, null, null, buff, null, true);
			ros.start();
		}
		if(outb1!=null){
			final int buff=(!ORDERED ? 12 : Tools.max(32, 2*Shared.THREADS));
			FileFormat ff1=FileFormat.testOutput(outb1, FileFormat.FASTQ, null, true, overwrite, ORDERED);
			FileFormat ff2=FileFormat.testOutput(outb2, FileFormat.FASTQ, null, true, overwrite, ORDERED);
			rosb=new RTextOutputStream3(ff1, ff2, null, null, buff, null, true);
			rosb.start();
		}
		if(ros!=null || rosb!=null){
			t.stop();
			outstream.println("Started output stream:\t"+t);
			t.start();
		}
		
		/* Optionally skip the first reads, since initial reads may have lower quality */
		if(skipreads>0){
			long skipped=0;

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(skipped<skipreads && reads!=null && reads.size()>0){
				skipped+=reads.size();
				if(rosb!=null){
					rosb.add(new ArrayList<Read>(1), ln.id);
				}
				
				if(ros!=null){
					ros.add(new ArrayList<Read>(1), ln.id);
				}
				
				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln, ln.list.isEmpty());
			if(reads==null || reads.isEmpty()){
				ReadWrite.closeStreams(cris, ros, rosb);
				System.err.println("Skipped all of the reads.");
				System.exit(0);
			}
		}
		
		/* Create ProcessThreads */
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new ProcessThread(cris, ros, rosb));}
		for(ProcessThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(ProcessThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			readsKFiltered+=pt.readsKFilteredT;
			basesKFiltered+=pt.basesKFilteredT;
			readsQTrimmed+=pt.readsQTrimmedT;
			basesQTrimmed+=pt.basesQTrimmedT;
			readsKTrimmed+=pt.readsKTrimmedT;
			basesKTrimmed+=pt.basesKTrimmedT;
			readsQFiltered+=pt.readsQFilteredT;
			basesQFiltered+=pt.basesQFilteredT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris, ros, rosb);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	
	/**
	 * Loads kmers into a table.  Each thread handles all kmers X such that X%WAYS==tnum.
	 */
	private class LoadThread extends Thread{
		
		public LoadThread(final int tnum_){
			tnum=tnum_;
			map=keySets[tnum];
		}
		
		/**
		 * Get the next list of reads (or scaffolds) from the queue.
		 * @return List of reads
		 */
		private ArrayList<Read> fetch(){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			return list;
		}
		
		@Override
		public void run(){
			ArrayList<Read> reads=fetch();
			while(reads!=POISON){
				for(Read r : reads){
					assert(r.pairnum()==0);
					final Read r2=r.mate;

					final int rblen=(r==null || r.bases==null ? 0 : r.bases.length);
					final int rblen2=(r2==null || r2.bases==null ? 0 : r2.bases.length);

					added+=addToMap(r, rblen>10000000 ? k : rblen>1000000 ? 11 : rblen>100000 ? 2 : 0);
					if(r.mate!=null){
						added+=addToMap(r.mate, rblen2>10000000 ? k : rblen2>1000000 ? 11 : rblen2>100000 ? 2 : 0);
					}
				}
				reads=fetch();
			}
			
			if(map.canRebalance() && map.size()>2L*map.arrayLength()){
				map.rebalance();
			}
		}

		/**
		 * @param r The current read to process
		 * @param skip Number of bases to skip between kmers
		 * @return Number of kmers stored
		 */
		private long addToMap(Read r, int skip){
			skip=Tools.max(minSkip, Tools.min(maxSkip, skip));
			final byte[] bases=r.bases;
			final long shift=2*k;
			final long shift2=shift-2;
			final long mask=~((-1L)<<shift);
			final long kmask=kMasks[k];
			long kmer=0;
			long rkmer=0;
			long added=0;
			int len=0;

			if(bases!=null){
				refReadsT++;
				refBasesT+=bases.length;
			}
			if(bases==null || bases.length<k){return -1;}
			
			final int id=(Integer)r.obj;

			if(skip>1){ //Process while skipping some kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=Dedupe.baseToNumber[b];
					long x2=Dedupe.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=(rkmer>>>2)|(x2<<shift2);
					if(b=='N'){len=0;}else{len++;}
					if(verbose){System.err.println("Scanning1 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						refKmersT++;
						if(len%skip==0){
							final long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
							added+=addToMap(kmer, rkmer, k, extraBase, id, kmask);
							if(useShortKmers){
								if(i==k2){added+=addToMapRightShift(kmer, rkmer, id);}
								if(i==bases.length-1){added+=addToMapLeftShift(kmer, rkmer, extraBase, id);}
							}
						}
					}
				}
			}else{ //Process all kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=Dedupe.baseToNumber[b];
					long x2=Dedupe.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=(rkmer>>>2)|(x2<<shift2);
					if(b=='N'){len=0;}else{len++;}
					if(verbose){System.err.println("Scanning2 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						refKmersT++;
						final long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
						final long atm=addToMap(kmer, rkmer, k, extraBase, id, kmask);
						added+=atm;
//						assert(false) : atm+", "+map.contains(toValue(kmer, rkmer, kmask));
						if(useShortKmers){
							if(i==k2){added+=addToMapRightShift(kmer, rkmer, id);}
							if(i==bases.length-1){added+=addToMapLeftShift(kmer, rkmer, extraBase, id);}
						}
					}
				}
			}
			return added;
		}
		

		/**
		 * Adds short kmers on the left end of the read.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @return Number of kmers stored
		 */
		private long addToMapLeftShift(long kmer, long rkmer, final long extraBase, final int id){
			if(verbose){System.err.println("addToMapLeftShift");}
			long added=0;
			for(int i=k-1; i>=mink; i--){
				kmer=kmer&rightMasks[i];
				rkmer=rkmer>>>2;
				long x=addToMap(kmer, rkmer, i, extraBase, id, kMasks[i]);
				added+=x;
				if(verbose){
					if((toValue(kmer, rkmer, kMasks[i]))%WAYS==tnum){
						System.err.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added left-shift kmer "+AminoAcid.kmerToString(kmer&~kMasks[i], i)+"; value="+(toValue(kmer, rkmer, kMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+kMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
						System.err.println("i="+i+"; tnum="+tnum+"; Looking for left-shift kmer "+AminoAcid.kmerToString(kmer&~kMasks[i], i));
						final long value=toValue(kmer, rkmer, kMasks[i]);
						if(map.contains(value)){System.err.println("Found "+value);}
					}
				}
			}
			return added;
		}
		

		/**
		 * Adds short kmers on the right end of the read.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @return Number of kmers stored
		 */
		private long addToMapRightShift(long kmer, long rkmer, final int id){
			if(verbose){System.err.println("addToMapRightShift");}
			long added=0;
			for(int i=k-1; i>=mink; i--){
				long extraBase=kmer&3L;
				kmer=kmer>>>2;
				rkmer=rkmer&rightMasks[i];
//				assert(Long.numberOfLeadingZeros(kmer)>=2*(32-i)) : Long.numberOfLeadingZeros(kmer)+", "+i+", "+kmer+", "+kMasks[i];
//				assert(Long.numberOfLeadingZeros(rkmer)>=2*(32-i)) : Long.numberOfLeadingZeros(rkmer)+", "+i+", "+rkmer+", "+kMasks[i];
				long x=addToMap(kmer, rkmer, i, extraBase, id, kMasks[i]);
				added+=x;
				if(verbose){
					if((toValue(kmer, rkmer, kMasks[i]))%WAYS==tnum){
						System.err.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added right-shift kmer "+AminoAcid.kmerToString(kmer&~kMasks[i], i)+"; value="+(toValue(kmer, rkmer, kMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+kMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
						System.err.println("i="+i+"; tnum="+tnum+"; Looking for right-shift kmer "+AminoAcid.kmerToString(kmer&~kMasks[i], i));
						final long value=toValue(kmer, rkmer, kMasks[i]);
						if(map.contains(value)){System.err.println("Found "+value);}
					}
				}
			}
			return added;
		}
		
		
		/**
		 * Adds this kmer to the table, including any mutations implied by editDistance or hammingDistance.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @param kmask0
		 * @return Number of kmers stored
		 */
		private long addToMap(final long kmer, final long rkmer, final int len, final long extraBase, final int id, final long kmask0){
			
			assert(kmask0==kMasks[len]) : kmask0+", "+len+", "+kMasks[len]+", "+Long.numberOfTrailingZeros(kmask0)+", "+Long.numberOfTrailingZeros(kMasks[len]);
			
			if(verbose){System.err.println("addToMap_A; len="+len+"; kMasks[len]="+kMasks[len]);}
			assert((kmer&kmask0)==0);
			final long added;
			if(hammingDistance==0){
				final long key=toValue(kmer, rkmer, kmask0);
				if(key%WAYS!=tnum){return 0;}
				if(verbose){System.err.println("addToMap_B: "+AminoAcid.kmerToString(kmer&~kMasks[len], len)+" = "+key);}
				added=map.setIfNotPresent(key, id);
			}else if(editDistance>0){
//				long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
				added=mutate(kmer, rkmer, len, id, editDistance, extraBase);
			}else{
				added=mutate(kmer, rkmer, len, id, hammingDistance, -1);
			}
			if(verbose){System.err.println("addToMap added "+added+" keys.");}
			return added;
		}
		
		/**
		 * Mutate and store this kmer through 'dist' recursions.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @param dist Number of mutations
		 * @param extraBase Base added to end in case of deletions
		 * @return Number of kmers stored
		 */
		private long mutate(final long kmer, final long rkmer, final int len, final int id, final int dist, final long extraBase){
			long added=0;
			
			final long key=toValue(kmer, rkmer, kMasks[len]);
			
			if(verbose){System.err.println("mutate_A; len="+len+"; kmer="+kmer+"; rkmer="+rkmer+"; kMasks[len]="+kMasks[len]);}
			if(key%WAYS==tnum){
				if(verbose){System.err.println("mutate_B: "+AminoAcid.kmerToString(kmer&~kMasks[len], len)+" = "+key);}
				int x=map.setIfNotPresent(key, id);
				if(verbose){System.err.println("mutate_B added "+x+" keys.");}
				added+=x;
				assert(map.contains(key));
			}
			
			if(dist>0){
				final int dist2=dist-1;
				
				//Sub
				for(int j=0; j<4; j++){
					for(int i=0; i<len; i++){
						final long temp=(kmer&clearMasks[i])|setMasks[j][i];
						if(temp!=kmer){
							long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
							added+=mutate(temp, rtemp, len, id, dist2, extraBase);
						}
					}
				}
				
				if(editDistance>0){
					//Del
					if(extraBase>=0 && extraBase<=3){
						for(int i=1; i<len; i++){
							final long temp=(kmer&leftMasks[i])|((kmer<<2)&rightMasks[i])|extraBase;
							if(temp!=kmer){
								long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, -1);
							}
						}
					}

					//Ins
					final long eb2=kmer&3;
					for(int i=1; i<len; i++){
						final long temp0=(kmer&leftMasks[i])|((kmer&rightMasks[i])>>2);
						for(int j=0; j<4; j++){
							final long temp=temp0|setMasks[j][i-1];
							if(temp!=kmer){
								long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, eb2);
							}
						}
					}
				}

			}
			
			return added;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Number of kmers stored by this thread */
		public long added=0;
		/** Number of items encountered by this thread */
		public long refKmersT=0, refReadsT=0, refBasesT=0;
		/** Thread number; used to determine which kmers to store */
		public final int tnum;
		/** Buffer of input read lists */
		public final ArrayBlockingQueue<ArrayList<Read>> queue=new ArrayBlockingQueue<ArrayList<Read>>(32);
		
		/** Destination for storing kmers */
		private final AbstractKmerTable map;
		
	}
	
	/*--------------------------------------------------------------*/

	/**
	 * Matches read kmers against reference kmers, performs binning and/or trimming, and writes output. 
	 */
	private class ProcessThread extends Thread{
		
		/**
		 * Constructor
		 * @param cris_ Read input stream
		 * @param ros_ Unmatched read output stream (optional)
		 * @param rosb_ Matched read output stream (optional)
		 */
		public ProcessThread(ConcurrentReadStreamInterface cris_, RTextOutputStream3 ros_, RTextOutputStream3 rosb_){
			cris=cris_;
			ros=ros_;
			rosb=rosb_;
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			ArrayList<Read> bad=(rosb==null ? null : new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
			
			//While there are more reads lists...
			while(reads!=null && reads.size()>0){
				
				int removed=0;
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;

					final int initialLength1=(r1.bases==null ? 0 : r1.bases.length);
					final int initialLength2=(r2==null ? 0 : r2.bases==null ? 0 : r2.bases.length);

					final int minlen1=(int)Tools.max(initialLength1*minLenFraction, minReadLength);
					final int minlen2=(int)Tools.max(initialLength2*minLenFraction, minReadLength);
					
					if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}
					
					readsInT++;
					basesInT+=initialLength1;
					if(r2!=null){
						readsInT++;
						basesInT+=initialLength2;
					}
					
					boolean remove=false;
					
					//Determine whether to discard the reads based on average quality
					if(minAvgQuality>0){
						if(r1!=null && r1.quality!=null && r1.avgQuality()<minAvgQuality){r1.setDiscarded(true);}
						if(r2!=null && r2.quality!=null && r2.avgQuality()<minAvgQuality){r2.setDiscarded(true);}
					}
					//Determine whether to discard the reads based on the presence of Ns
					if(maxNs>=0){
						if(r1!=null && r1.countUndefined()>maxNs){r1.setDiscarded(true);}
						if(r2!=null && r2.countUndefined()>maxNs){r2.setDiscarded(true);}
					}

					if(removePairsIfEitherBad){remove=r1.discarded() || (r2!=null && r2.discarded());}
					else{remove=r1.discarded() && (r2==null || r2.discarded());}

					if(remove){
						if(r1!=null){
							basesQFilteredT+=r1.bases.length;
							readsQFilteredT++;
						}
						if(r2!=null){
							basesQFilteredT+=r2.bases.length;
							readsQFilteredT++;
						}
						if(bad!=null){bad.add(r1);}
					}else{
						//Process kmers
						
						if(ktrimLeft || ktrimRight || ktrimN){
							
							int rlen1=0, rlen2=0;
							int xsum=0;
							int rktsum=0;
							if(r1!=null){
								int x=ktrim(r1, keySets);
								xsum+=x;
								rktsum+=(x>0 ? 1 : 0);
								rlen1=r1.bases==null ? 0 : r1.bases.length;
								if(rlen1<minlen1){r1.setDiscarded(true);}
							}
							if(r2!=null){
								int x=ktrim(r2, keySets);
								xsum+=x;
								rktsum+=(x>0 ? 1 : 0);
								rlen2=r2.bases==null ? 0 : r2.bases.length;
								if(rlen2<minlen2){r2.setDiscarded(true);}
							}
							
							if((removePairsIfEitherBad && (r1.discarded() || (r2!=null && r2.discarded()))) || 
									(r1.discarded() && (r2==null || r2.discarded()))){
								if(!ktrimN){
									xsum+=(rlen1+rlen2);
									rktsum=(r1==null ? 0 : 1)+(r2==null ? 0 : 1);
								}
								remove=true;
								if(addTrimmedToBad && bad!=null){bad.add(r1);}
							}
							basesKTrimmedT+=xsum;
							readsKTrimmedT+=rktsum;
							
						}else{
							//Do kmer matching
							
							final int a=(kbig<=k ? countSetKmers(r1, keySets) : countSetKmersBig(r1, keySets));
							final int b=(kbig<=k ? countSetKmers(r2, keySets) : countSetKmersBig(r2, keySets));
							
							if(r1!=null && a>maxBadKmers){r1.setDiscarded(true);}
							if(r2!=null && b>maxBadKmers){r2.setDiscarded(true);}
							
							
							if((removePairsIfEitherBad && (r1.discarded() || (r2!=null && r2.discarded()))) || 
									(r1.discarded() && (r2==null || r2.discarded()))){
								remove=true;
								if(a>maxBadKmers){
									readsKFilteredT++;
									basesKFilteredT+=r1.bases.length;
								}
								if(b>maxBadKmers){
									readsKFilteredT++;
									basesKFilteredT+=r2.bases.length;
								}
								if(bad!=null){bad.add(r1);}
							}
//							assert(false) : "a="+a+", b="+b+", r1.discarded()="+r1.discarded()+", removePairsIfEitherBad="+removePairsIfEitherBad+
//								"\nremove="+remove+", readsKFilteredT="+readsKFilteredT;
						}
					}
					
					if(!remove){
						//Do quality trimming
						
						int rlen1=0, rlen2=0;
						if(r1!=null){
							if(qtrimLeft || qtrimRight){
								int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, 1);
								basesQTrimmedT+=x;
								readsQTrimmedT+=(x>0 ? 1 : 0);
							}
							rlen1=r1.bases==null ? 0 : r1.bases.length;
							if(rlen1<minlen1){r1.setDiscarded(true);}
						}
						if(r2!=null){
							if(qtrimLeft || qtrimRight){
								int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, 1);
								basesQTrimmedT+=x;
								readsQTrimmedT+=(x>0 ? 1 : 0);
							}
							rlen2=r2.bases==null ? 0 : r2.bases.length;
							if(rlen2<minlen2){r2.setDiscarded(true);}
						}
						
						//Discard reads if too short
						if((removePairsIfEitherBad && (r1.discarded() || (r2!=null && r2.discarded()))) || 
								(r1.discarded() && (r2==null || r2.discarded()))){
							basesQTrimmedT+=((r1.bases==null ? 0 : r1.bases.length)+(r2==null || r2.bases==null ? 0 : r2.bases.length));
							remove=true;
							if(addTrimmedToBad && bad!=null){bad.add(r1);}
						}
					}
					
					if(remove){
						//Evict read
						removed++;
						if(r2!=null){removed++;}
						reads.set(i, null);
//						System.err.println("X1\t"+removed);
					}else{
						//Track statistics
						
						if(r1!=null){
							readsOutT++;
							basesOutT+=r1.bases==null ? 0 : r1.bases.length;
						}
						if(r2!=null){
							readsOutT++;
							basesOutT+=r2.bases==null ? 0 : r2.bases.length;
						}
//						System.err.println("X2\t"+readsOutT);
					}
				}
				
				//Send matched list to matched output stream
				if(rosb!=null){
					rosb.add(bad, ln.id);
					bad.clear();
				}
				
				//Send unmatched list to unmatched output stream
				if(ros!=null){
					ros.add((removed>0 ? Tools.condenseNew(reads) : reads), ln.id); //Creates a new list if old one became empty, to prevent shutting down the cris.
				}
				
				//Fetch a new read list
				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln, ln.list.isEmpty());
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadStreamInterface cris;
		/** Output read streams */
		private final RTextOutputStream3 ros, rosb;
		
		private long readsInT=0;
		private long basesInT=0;
		private long readsOutT=0;
		private long basesOutT=0;
		
		private long readsQTrimmedT=0;
		private long basesQTrimmedT=0;
		private long readsQFilteredT=0;
		private long basesQFilteredT=0;

		private long readsKTrimmedT=0;
		private long basesKTrimmedT=0;
		private long readsKFilteredT=0;
		private long basesKFilteredT=0;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Trim a read to remove matching kmers and everything to their left or right.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of bases trimmed
	 */
	private final int ktrim(final Read r, final AbstractKmerTable[] sets){
		assert(ktrimLeft || ktrimRight || ktrimN);
		if(r==null || r.bases==null){return 0;}
		if(verbose){System.err.println("KTrimming read "+r.id);}
		final byte[] bases=r.bases, quals=r.quality;
		final int minlen=k-1;
		final int minlen2=(maskMiddle ? k/2 : k);
		final long shift=2*k;
		final long shift2=shift-2;
		final long mask=~((-1L)<<shift);
		final long kmask=kMasks[k];
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;
		
		int minLoc=999999999, minLocExclusive=999999999;
		int maxLoc=-1, maxLocExclusive=-1;
		final int initialLength=r.bases.length;
		
		if(bases==null || bases.length<k){return 0;}
		
		//Scan for normal kmers
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			long x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(b=='N' && forbidNs){len=0;}else{len++;}
			if(verbose){System.err.println("Scanning3 i="+i+", kmer="+kmer+", rkmer="+rkmer+", len="+len+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				final long value=toValue(kmer, rkmer, kmask);
				AbstractKmerTable set=sets[(int)(value%WAYS)];
				if(verbose){System.err.println("set.contains("+value+")="+set.contains(value));}
				if(set.contains(value)){
					minLoc=Tools.min(minLoc, i-k+1);
					assert(minLoc>=0);
					maxLoc=i;
					found++;
				}
			}
		}

		if(minLoc!=minLocExclusive){minLocExclusive=minLoc+k;}
		if(maxLoc!=maxLocExclusive){maxLocExclusive=maxLoc-k;}
		
		//If nothing was found, scan for short kmers.  Only used for trimming.
		if(useShortKmers && found==0){
			assert(!maskMiddle && middleMask==-1) : maskMiddle+", "+middleMask+", k="+", mink="+mink;
			
			//Look for short kmers on left side
			if(ktrimLeft || ktrimN){
				kmer=0;
				rkmer=0;
				len=0;
				for(int i=0; i<k && i<bases.length; i++){
					byte b=bases[i];
					long x=Dedupe.baseToNumber[b];
					long x2=Dedupe.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=rkmer|(x2<<(2*len));
					len++;
					if(verbose){System.err.println("Scanning4 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=mink){
						
						if(verbose){
							System.err.println("Looking for left kmer  "+AminoAcid.kmerToString(kmer, len));
							System.err.println("Looking for left rkmer "+AminoAcid.kmerToString(rkmer, len));
						}
						final long value=toValue(kmer, rkmer, kMasks[len]);
						AbstractKmerTable set=sets[(int)(value%WAYS)];
						if(set.contains(value)){
							if(verbose){System.err.println("Found "+value);}
							minLoc=0;
							minLocExclusive=Tools.min(minLocExclusive, i+1);
							maxLoc=Tools.max(maxLoc, i);
							maxLocExclusive=Tools.max(maxLocExclusive, 0);
							found++;
						}
					}
				}
			}

			//Look for short kmers on right side
			if(ktrimRight || ktrimN){
				kmer=0;
				rkmer=0;
				len=0;
				for(int i=bases.length-1; i>=0 && i>bases.length-k; i--){
					byte b=bases[i];
					long x=Dedupe.baseToNumber[b];
					long x2=Dedupe.baseToComplementNumber[b];
					kmer=kmer|(x<<(2*len));
					rkmer=((rkmer<<2)|x2)&mask;
					len++;
					if(verbose){System.err.println("Scanning5 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=mink){
						
						final long value=toValue(kmer, rkmer, kMasks[len]);
						if(verbose){System.err.println("Looking for right kmer "+AminoAcid.kmerToString(kmer&~kMasks[len], len)+"; value="+value+"; kmask="+kMasks[len]);}
						AbstractKmerTable set=sets[(int)(value%WAYS)];
						if(set.contains(value)){
							if(verbose){System.err.println("Found "+value);}
							minLoc=i;
							minLocExclusive=Tools.min(minLocExclusive, bases.length);
							maxLoc=bases.length-1;
							maxLocExclusive=Tools.max(maxLocExclusive, i-1);
							found++;
						}
					}
				}
			}
		}
		

		if(verbose){System.err.println("found="+found+", minLoc="+minLoc+", maxLoc="+maxLoc+", minLocExclusive="+minLocExclusive+", maxLocExclusive="+maxLocExclusive);}
		
		
		if(found==0){return 0;}
		
		if(ktrimN){ //Replace kmer hit zone with the trim symbol
			Arrays.fill(bases, minLoc, maxLoc+1, trimSymbol);
			if(quals!=null){Arrays.fill(quals, minLoc, maxLoc+1, (byte)0);}
			return maxLoc-minLoc+1;
		}
		
		if(ktrimLeft){ //Trim from the read start to the rightmost kmer base
			if(verbose){System.err.println("Left trimming to "+(ktrimExclusive ? maxLocExclusive+1 : maxLoc+1)+", "+0);}
			int x=TrimRead.trimToPosition(r, ktrimExclusive ? maxLocExclusive+1 : maxLoc+1, 0, 1);
			if(verbose){System.err.println("Trimmed "+x+" bases: "+new String(r.bases));}
			return x;
		}else{ //Trim from the leftmost kmer base to the read stop 
			assert(ktrimRight);
			if(verbose){System.err.println("Right trimming to "+0+", "+(ktrimExclusive ? minLocExclusive-1 : minLoc-1));}
			int x=TrimRead.trimToPosition(r, 0, ktrimExclusive ? minLocExclusive-1 : minLoc-1, 1);
			if(verbose){System.err.println("Trimmed "+x+" bases: "+new String(r.bases));}
			return x;
		}
	}
	
	/**
	 * Counts the number of kmer hits for a read.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of hits
	 */
	private final int countSetKmers(final Read r, final AbstractKmerTable sets[]){
		if(r==null || r.bases==null){return 0;}
		final byte[] bases=r.bases;
		final int minlen=k-1;
		final int minlen2=(maskMiddle ? k/2 : k);
		final long shift=2*k;
		final long shift2=shift-2;
		final long mask=~((-1L)<<shift);
		final long kmask=kMasks[k];
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;

		if(bases==null || bases.length<k){return -1;}
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			long x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(b=='N' && forbidNs){len=0;}else{len++;}
			if(verbose){System.err.println("Scanning6 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				final long key=toValue(kmer, rkmer, kmask);
				if(verbose){System.err.println("Testing key "+key);}
				AbstractKmerTable set=sets[(int)(key%WAYS)];
				final int id=set.getCount(key);
				if(id>0){
					if(verbose){System.err.println("Found = "+(found+1)+"/"+maxBadKmers);}
					if(found==maxBadKmers){
						int count=scaffoldCounts.addAndGet(id, 1);
						if(count<0){scaffoldCounts.set(id, Integer.MAX_VALUE);}
						if(hitCounts==null){
							return (found=found+1);
						}//Early exit, but prevents generation of histogram that goes over maxBadKmers+1.
					}
					found++;
					
//					if(found>maxBadKmers){
//						Integer id=set.get(key);
//						int count=scaffoldCounts.addAndGet(id, 1);
//						if(count<0){scaffoldCounts.set(id, Integer.MAX_VALUE);}
//						break;
//					}
				}
			}
		}
		
		if(hitCounts!=null){hitCounts.incrementAndGet(Tools.min(found, hitCounts.length()-1));}
		return found;
	}
	
	/** Estimates kmer hit counts for kmers longer than k using consecutive matches
	 * @param r
	 * @param sets
	 * @return Number of sets of consecutive hits of exactly length kbig
	 */
	private final int countSetKmersBig(final Read r, final AbstractKmerTable sets[]){
		assert(kbig>k);
		final int sub=kbig-k-1;
		assert(sub>=0) : kbig+", "+sub;
		if(r==null || r.bases==null){return 0;}
		final byte[] bases=r.bases;
		final int minlen=k-1;
		final int minlen2=(maskMiddle ? k/2 : k);
		final long shift=2*k;
		final long shift2=shift-2;
		final long mask=~((-1L)<<shift);
		final long kmask=kMasks[k];
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;

		int start=-1;
		int stop=-1;
		int id=-1;
		
		if(bases==null || bases.length<k){return -1;}

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			long x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(b=='N' && forbidNs){len=0;}else{len++;}
			if(verbose){System.err.println("Scanning7 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				final long key=toValue(kmer, rkmer, kmask);
				AbstractKmerTable set=sets[(int)(key%WAYS)];
				id=set.getCount(key);
				if(id>0){
					if(start==-1){start=i;}
					stop=i;
				}else{
					if(start>-1){
						int dif=stop-start-sub;
						stop=start=-1;
						if(dif>0){
							int old=found;
							found+=dif;
							if(found>maxBadKmers && old<=maxBadKmers){
								int count=scaffoldCounts.addAndGet(id, 1);
								if(count<0){scaffoldCounts.set(id, Integer.MAX_VALUE);}
								if(hitCounts==null){
									return found;
								}//Early exit, but prevents generation of histogram that goes over maxBadKmers+1.
							}
						}
					}
				}
			}
		}
		
		// This catches the case where valid kmers extend to the end of the read
		if(start>-1){
			int dif=stop-start-sub;
			stop=start=-1;
			if(dif>0){
				int old=found;
				found+=dif;
				if(found>maxBadKmers && old<=maxBadKmers){
					int count=scaffoldCounts.addAndGet(id, 1);
					if(count<0){scaffoldCounts.set(id, Integer.MAX_VALUE);}
				}
			}
		}
		
		if(hitCounts!=null){hitCounts.incrementAndGet(Tools.min(found, hitCounts.length()-1));}
		return found;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Object holding a String and a number, for tracking the number of kmer hits per scaffold.
	 */
	private static class StringNum implements Comparable<StringNum>{
		
		public StringNum(String s_, int n_){
			s=s_;
			n=n_;
		}
		public final int compareTo(StringNum o){
			if(n!=o.n){return o.n-n;}
			return s.compareTo(o.s);
		}
		public final boolean equals(StringNum o){
			return compareTo(o)==0;
		}
		public final String toString(){
			return s+"\t"+n;
		}
		
		/*--------------------------------------------------------------*/
		
		public final String s;
		public final int n;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Deprecated
	/** Now handled by Read class */
	private static final void toUpperCase(Read r){
		final byte[] bases1=r.bases, bases2=(r.mate==null ? null : r.mate.bases);
		if(bases1!=null){
			for(int i=0; i<bases1.length; i++){
				bases1[i]=(byte)Character.toUpperCase(bases1[i]);
			}
		}
		if(bases2!=null){
			for(int i=0; i<bases2.length; i++){
				bases2[i]=(byte)Character.toUpperCase(bases2[i]);
			}
		}
	}
	
	/** Print statistics about current memory use and availability */
	private static final void printMemory(){
		if(GC_BEFORE_PRINT_MEMORY){
			System.gc();
			System.gc();
		}
		Runtime rt=Runtime.getRuntime();
		long mmemory=rt.maxMemory()/1000000;
		long tmemory=rt.totalMemory()/1000000;
		long fmemory=rt.freeMemory()/1000000;
		long umemory=tmemory-fmemory;
		outstream.println("Memory: "+/*"max="+mmemory+"m, total="+tmemory+"m, "+*/"free="+fmemory+"m, used="+umemory+"m");
	}
	
	/**
	 * Transforms a kmer into a canonical value stored in the table.  Expected to be inlined.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @return Canonical value
	 */
	private final long toValue(long kmer, long rkmer, long lengthMask){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return (value&middleMask)|lengthMask;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;
	
	/** Initial size of data structures */
	private final int initialSize=128000;
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private final AbstractKmerTable[] keySets=new AbstractKmerTable[WAYS];
	/** A scaffold's name is stored at scaffoldNames.get(id).  
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	private final ArrayList<String> scaffoldNames=new ArrayList<String>();
	/** scaffoldCounts[id] stores the number of reads with kmer matches to that scaffold */ 
	private AtomicIntegerArray scaffoldCounts;
	/** hitCounts[x] stores the number of reads with exactly x kmer matches */
	private final AtomicLongArray hitCounts;
	/** Array of reference files from which to load kmers */
	private String[] ref=null;
	/** Array of literal strings from which to load kmers */
	private String[] literal=null;
	
	/** Input reads */
	private String in1=null, in2=null;
	/** Output reads (unmatched and at least minlen) */
	private String out1=null, out2=null;
	/** Output reads (matched or shorter than minlen) */
	private String outb1=null, outb2=null;
	/** Statistics output files */
	private String outstats=null, outduk=null, outrqc=null;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	private long maxReads=-1;
	/** Output reads in input order.  May reduce speed. */
	private boolean ORDERED=false;
	/** Attempt to match kmers shorter than normal k on read ends when doing kTrimming. */
	private boolean useShortKmers=false;
	/** Make the middle base in a kmer a wildcard to improve sensitivity */
	private boolean maskMiddle=true;
	
	/** Store reference kmers with up to this many substitutions */
	private int hammingDistance=0;
	/** Store reference kmers with up to this many edits (including indels) */
	private int editDistance=0;
	/** Never skip more than this many consecutive kmers when hashing reference. */
	private int maxSkip=99;
	/** Always skip at least this many consecutive kmers when hashing reference.
	 * Note that a skip of 1 means every kmer is used, 2 means every other, etc. */
	private int minSkip=1;
	
	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/
	
	long readsIn=0;
	long basesIn=0;
	long readsOut=0;
	long basesOut=0;
	
	long readsQTrimmed=0;
	long basesQTrimmed=0;
	long readsQFiltered=0;
	long basesQFiltered=0;
	
	long readsKTrimmed=0;
	long basesKTrimmed=0;
	long readsKFiltered=0;
	long basesKFiltered=0;
	
	long refReads=0;
	long refBases=0;
	long refKmers=0;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** Don't allow a read 'N' to match a reference 'A'.  
	 * Reduces sensitivity when hdist>0 or edist>0.  Default: false. */
	private final boolean forbidNs;
	/** AND bitmask with 0's at the middle base */ 
	private final long middleMask;
	/** Use HashForest data structure */
	private final boolean useForest;
	/** Use KmerTable data structure */
	private final boolean useTable;
	/** Use HashArray data structure (default) */
	private final boolean useArray;	
	
	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	/** Emulated kmer greater than k */
	private final int kbig;
	/** Shortest kmer to use for trimming */
	private final int mink;
//	/** OR bitmask for full-length kmers, with a single 1 bit immediately to the left of the kmer */
//	private final long kmask; //TODO: Eliminate this and just use kmasks array
	/** A read may contain up to this many kmers before being considered a match.  Default: 0 */
	private final int maxBadKmers;
	
	/** Quality-trim the left side */
	private final boolean qtrimLeft;
	/** Quality-trim the right side */
	private final boolean qtrimRight;
	/** Trim bases at this quality or below.  Default: 4 */
	private final byte trimq;
	/** Throw away reads below this average quality before trimming.  Default: 0 */
	private final byte minAvgQuality;
	/** Throw away reads containing more than this many Ns.  Default: -1 (disabled) */
	private final int maxNs;
	/** Throw away reads shorter than this after trimming.  Default: 20 */
	private final int minReadLength;
	/** Toss reads shorter than this fraction of initital length, after trimming */
	private final float minLenFraction;
	/** Trim matching kmers and all bases to the left */
	private final boolean ktrimLeft;
	/** Trim matching kmers and all bases to the right */
	private final boolean ktrimRight;
	/** Don't trim, but replace matching kmers with a symbol (default N) */
	private final boolean ktrimN;
	/** Exclude kmer itself when ktrimming */
	private final boolean ktrimExclusive;
	/** Replace bases covered by matched kmers with this symbol */
	private final byte trimSymbol;
	/** Output over-trimmed reads to outbad (outmatch).  If false, they are discarded. */
	private final boolean addTrimmedToBad;
	
	/** True iff java was launched with the -ea' flag */
	private final boolean EA;
	/** Skip this many initial input reads */
	private final long skipreads;

	/** Pairs go to outbad if either of them is bad, as opposed to requiring both to be bad.
	 * Default: true. */
	private final boolean removePairsIfEitherBad;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int VERSION=8;
	
	/** Number of tables (and threads, during loading) */ 
	private static final int WAYS=5; //123
	/** Verbose messages */
	public static final boolean verbose=false; //123
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.THREADS;
	/** Indicates end of input stream */
	private static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	/** Do garbage collection prior to printing memory usage */
	private static final boolean GC_BEFORE_PRINT_MEMORY=false;
	
	/** x&clearMasks[i] will clear base i */
	private static final long[] clearMasks;
	/** x|setMasks[i][j] will set base i to j */
	private static final long[][] setMasks;
	/** x&leftMasks[i] will clear all bases to the right of i (exclusive) */
	private static final long[] leftMasks;
	/** x&rightMasks[i] will clear all bases to the left of i (inclusive) */
	private static final long[] rightMasks;
	/** x|kMasks[i] will set the bit to the left of the leftmost base */
	private static final long[] kMasks;
	
	public static HashMap<String,String> RQC_MAP=null;
	
	/*--------------------------------------------------------------*/
	/*----------------      Static Initializers     ----------------*/
	/*--------------------------------------------------------------*/
	
	static{
		clearMasks=new long[32];
		leftMasks=new long[32];
		rightMasks=new long[32];
		kMasks=new long[32];
		setMasks=new long[4][32];
		for(int i=0; i<32; i++){
			clearMasks[i]=~(3L<<(2*i));
		}
		for(int i=0; i<32; i++){
			leftMasks[i]=((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			rightMasks[i]=~((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			kMasks[i]=((1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			for(long j=0; j<4; j++){
				setMasks[(int)j][i]=(j<<(2*i));
			}
		}
	}
	
}
