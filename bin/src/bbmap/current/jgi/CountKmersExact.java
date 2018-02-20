package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import kmer.AbstractKmerTable;
import kmer.HashArray;
import kmer.HashBuffer;
import kmer.HashForest;
import kmer.KCountArray;
import kmer.KmerCount7MTA;
import kmer.KmerTable;
import kmer.Primes;

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
 * @author Brian Bushnell
 * @date Nov 22, 2013
 *
 */
public class CountKmersExact {
	
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
		CountKmersExact cke=new CountKmersExact(args);
		
		///And run it
		cke.process();
	}
	
	/**
	 * Display usage information.
	 */
	private static void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("\njava -ea -Xmx20g -cp <path> jgi.CountKmersExact in=<input file>");
		outstream.println("\nOptional flags:");
		outstream.println("in=<file>          \tThe 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.");
		outstream.println("in2=<file>         \tUse this if 2nd read of pairs are in a different file.");
		outstream.println("out=<file>         \t(outnonmatch) The 'out=' flag is needed if the output file is not the second parameter.  'out=stdout' will pipe to standard out.");
		outstream.println("out2=<file>        \t(outnonmatch2) Use this to write 2nd read of pairs to a different file.");
		outstream.println("stats=<file>       \tWrite statistics about which contaminants were detected.");
		outstream.println("");
		outstream.println("threads=auto       \t(t) Set number of threads to use; default is number of logical processors.");
		outstream.println("overwrite=t        \t(ow) Set to false to force the program to abort rather than overwrite an existing file.");
		outstream.println("showspeed=t        \t(ss) Set to 'f' to suppress display of processing speed.");
		outstream.println("interleaved=auto   \t(int) If true, forces fastq input to be paired and interleaved.");
		outstream.println("k=28               \tKmer length used for finding contaminants.  Contaminants shorter than k will not be found.");
		outstream.println("minavgquality=0    \t(maq) Reads with average quality (before trimming) below this will be discarded.");
		outstream.println("touppercase=f      \t(tuc) Change all letters in reads and reference to upper-case.");
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
	public CountKmersExact(String[] args){
		for(String s : args){if(s.contains("standardout") || s.contains("stdout")){outstream=System.err;}}
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		FastaReadInputStream.SPLIT_READS=false;
		ByteFile.FORCE_MODE_BF2=Shared.THREADS>2;
		
		/* Initialize local variables with defaults */
		boolean setOut=false, setOutb=false, qtrimRight_=false, qtrimLeft_=false;
		boolean rcomp_=true;
		boolean useForest_=false, useTable_=false, useArray_=true;
		long skipreads_=0;
		int k_=28;
		int ways_=-1;
		byte qin=-1, qout=-1;
		
		byte trimq_=4;
		byte minAvgQuality_=0;
		int minReadLength_=20;
		int maxNs_=-1;
		boolean removePairsIfEitherBad_=true;
		int filterMax_=2;
		
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
			}else if(a.equals("stats")){
				outstats=b;
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
				
			}else if(a.equals("initialsize")){
				initialSize=Integer.parseInt(b);
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
			}else if(a.equals("buflen") || a.equals("bufflen") || a.equals("bufferlength")){
				buflen=Integer.parseInt(b);
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.equals("k")){
				assert(b!=null) : "\nThe k key needs an integer value from 1 to 31, such as k=28\n";
				k_=Integer.parseInt(b);
			}else if(a.equals("skipreads")){
				skipreads_=Long.parseLong(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.THREADS : Integer.parseInt(b));
			}else if(a.equals("minavgquality") || a.equals("maq")){
				minAvgQuality_=Byte.parseByte(b);
			}else if(a.equals("maxns")){
				maxNs_=Byte.parseByte(b);
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
//				verbose=Tools.parseBoolean(b);
			}else if(a.equals("rcomp")){
				rcomp_=Tools.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
			}else if(a.equals("fastaminlen") || a.equals("fastaminlength")){
				FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
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
			}else if(a.equals("prefilter")){
				if(b==null || b.length()<1 || !Character.isDigit(b.charAt(0))){
					prefilter=Tools.parseBoolean(b);
				}else{
					filterMax_=Integer.parseInt(b);
					prefilter=filterMax_>0;
				}
			}else if(i==0 && in1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				in1=args[i];
			}else if(i==1 && out1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				out1=args[i];
				setOut=true;
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
		if(ways_<1){
			long maxKmers=Runtime.getRuntime().maxMemory()/12;
			long minWays=Tools.min(10000, maxKmers/Integer.MAX_VALUE);
			ways_=(int)Tools.max(31, THREADS*4, minWays);
			ways_=(int)Primes.primeAtLeast(ways_);
			assert(ways_>0);
			System.err.println("ways="+ways_);
		}
		
		/* Set final variables; post-process and validate argument combinations */
		
		useForest=useForest_;
		useTable=useTable_;
		useArray=useArray_;
		rcomp=rcomp_;
		skipreads=skipreads_;
		trimq=trimq_;
		minAvgQuality=minAvgQuality_;
		minReadLength=minReadLength_;
		removePairsIfEitherBad=removePairsIfEitherBad_;
		maxNs=maxNs_;
		WAYS=ways_;
		filterMax=Tools.min(filterMax_, 0x7FFFFFFF);
		
		k=k_;
		k2=k-1;
		
		qtrimRight=qtrimRight_;
		qtrimLeft=qtrimLeft_;
		
		keySets=new AbstractKmerTable[WAYS];
		
		
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
		assert(!in1.equalsIgnoreCase(in2));
		assert(out1==null || !out1.equalsIgnoreCase(out2));
		assert(THREADS>0);

		assert(in1==null || in1.toLowerCase().startsWith("stdin") || in1.toLowerCase().startsWith("standardin") || new File(in1).exists()) : "Can't find "+in1;
		assert(in2==null || in2.toLowerCase().startsWith("stdin") || in2.toLowerCase().startsWith("standardin") || new File(in2).exists()) : "Can't find "+in2;
		
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
		Tools.testOutputFiles(overwrite, false, out1, out2,outstats);
		
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
		
		prefilterArray=null;
//		outstream.println();
		if(prefilter){
			
			long precells=-1;
			int prehashes=-1;
			int cbits=2;
			while(filterMax+1>=(1<<cbits)){cbits*=2;}
			byte minq=0;
			if(precells<1){
				Runtime rt=Runtime.getRuntime();
				final long memory=rt.maxMemory();
				long prebits=(long)(((memory-10)*8)*0.2f);
				precells=prebits/cbits;
			}
			if(prehashes<1){
				prehashes=2;
			}
			
			if(precells>100000){
				Timer ht=new Timer();
				ht.start();

				KmerCount7MTA.CANONICAL=true;
				ArrayList<String> extra=null;
				prefilterArray=KmerCount7MTA.makeKca(in1, in2, extra, k, cbits, 0, precells, prehashes, minq, true, maxReads, 1, 1, 1, 1, null);
				assert(filterMax<prefilterArray.maxValue);
				outstream.println("Made prefilter:   \t"+prefilterArray.toShortString(prehashes));
				double uf=prefilterArray.usedFraction();
				if(uf>0.6){
					outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" : 
						uf>0.7 ? "fairly" : "somewhat")+" full, which may reduce accuracy for kmers of depth under 3.  Ideal load is under 60% used." +
						"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; " +
							"or increase the values of the minprob flag to reduce spurious kmers.");
				}
				ht.stop();

				if(DISPLAY_PROGRESS){
					outstream.println("Prefilter time:\t"+ht);
					outstream.println("After prefilter:");
					printMemory();
					outstream.println();
				}
			}else{prefilter=false;}
		}
		
		/* Fill tables with kmers */
		long added=loadKmers(t);
		
		if(DISPLAY_PROGRESS){
			outstream.println("Final:");
			printMemory();
			outstream.println();
		}
		
		/* Write statistics to files */
		writeStats(System.nanoTime()-startTime);
		
		outstream.println("Input:                  \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		
		if(qtrimLeft || qtrimRight){
			outstream.println("QTrimmed:               \t"+readsTrimmed+" reads ("+String.format("%.2f",readsTrimmed*100.0/readsIn)+"%) \t"+
					basesTrimmed+" bases ("+String.format("%.2f",basesTrimmed*100.0/basesIn)+"%)");
		}
		if(minAvgQuality>0 || maxNs>=0){
			outstream.println("Low quality discards:   \t"+lowqReads+" reads ("+String.format("%.2f",lowqReads*100.0/readsIn)+"%) \t"+
					lowqBases+" bases ("+String.format("%.2f",lowqBases*100.0/basesIn)+"%)");
		}
		outstream.println("Result:                 \t"+readsOut+" reads ("+String.format("%.2f",readsOut*100.0/readsIn)+"%) \t"+
				basesOut+" bases ("+String.format("%.2f",basesOut*100.0/basesIn)+"%)");
		outstream.println("Unique Kmers:               \t"+added);
	}
	
	/**
	 * Write processing statistics in DUK's format.
	 * @param time Elapsed time, nanoseconds
	 */
	private void writeStats(long time){
		if(outstats==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outstats, overwrite, false, false);
		tsw.start();
		tsw.println(dukString(time));
		tsw.poisonAndWait();
	}
	
	/**
	 * Helper method; formats statistics to be duk-compatible
	 * @param time Elapsed time, nanoseconds
	 * @return duk output string
	 */
	private String dukString(long time){
		StringBuilder sb=new StringBuilder();
		sb.append("##INPUT PARAMETERS##\n");
		sb.append("#Query file:	"+in1+(in2==null ? "" : ","+in2)+"\n");
		sb.append("#Not matched reads file:	"+out1+(out2==null ? "" : ","+out2)+"\n");
		sb.append("#Output file (stats):	"+outstats+"\n");
		sb.append("#Mer size:	"+k+"\n");
		long size=0;
		for(AbstractKmerTable x : keySets){size+=x.size();}
		sb.append("#Quality trim:	"+((qtrimLeft || qtrimRight) ? trimq : "false")+"\n");
		sb.append("\n");
		
		sb.append("##REFERENCE STAT##\n");
//		sb.append("#Total Reads:	"+refReads+"\n");
//		sb.append("#Total Bases:	"+refBases+"\n");
//		sb.append("#Total kmers:	"+refKmers+"\n");
		sb.append("#Total stored kmers:	"+size+"\n");
		sb.append("\n");

		sb.append("## ELAPSED TIME##\n");
		sb.append("# Time:	"+String.format("%.2f", time/1000000000.0)+" seconds\n");
		sb.append("\n");

		sb.append("##QUERY FILE STAT##\n");
		sb.append("# Total number of reads:	"+readsIn+"\n");
		sb.append("\n");

		sb.append("##P-VALUE##\n");
		sb.append("#Avg number of Kmer for each read:	"+((basesIn/(Tools.max(readsIn, 1)))-k)+"\n");
//		sb.append("# P value for the given threshold 1 is 4.05231e-14\n"); //duk prints a P value; not sure what it means
		sb.append("\n");
		
		return sb.toString();
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	/**
	 * load reads into tables, using multiple ProcessThread.
	 * @param t
	 */
	private long loadKmers(Timer t){
		
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
		RTextOutputStream3 ros=null;
		if(out1!=null){
			final int buff=(!ORDERED ? 12 : Tools.max(32, 2*Shared.THREADS));
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, ORDERED);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, ORDERED);
			ros=new RTextOutputStream3(ff1, ff2, null, null, buff, null, true);
			ros.start();
		}
		if(ros!=null){
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
				
				if(ros!=null){
					ros.add(new ArrayList<Read>(1), ln.id);
				}
				
				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln, ln.list.isEmpty());
			if(reads==null || reads.isEmpty()){
				ReadWrite.closeStreams(cris, ros);
				System.err.println("Skipped all of the reads.");
				System.exit(0);
			}
		}
		
		/* Create ProcessThreads */
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new ProcessThread(cris));}
		for(ProcessThread pt : alpt){pt.start();}
		
		long added=0;
		
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
			added+=pt.added;
			
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			readsTrimmed+=pt.readsTrimmedT;
			basesTrimmed+=pt.basesTrimmedT;
			lowqReads+=pt.lowqReadsT;
			lowqBases+=pt.lowqBasesT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris, ros);
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
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
		public ProcessThread(ConcurrentReadStreamInterface cris_){
			cris=cris_;
			table=new HashBuffer(keySets, buflen);
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(reads!=null && reads.size()>0){
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;
					
					if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}
					
					readsInT++;
					basesInT+=(r1.bases==null ? 0 : r1.bases.length);
					if(r2!=null){
						readsInT++;
						basesInT+=(r2.bases==null ? 0 : r2.bases.length);
					}
					
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
					
					int rlen1=0, rlen2=0;
					if(r1!=null){
						if(qtrimLeft || qtrimRight){
							int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						rlen1=r1.bases==null ? 0 : r1.bases.length;
						if(rlen1<minReadLength){r1.setDiscarded(true);}
					}
					if(r2!=null){
						if(qtrimLeft || qtrimRight){
							int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						rlen2=r2.bases==null ? 0 : r2.bases.length;
						if(rlen2<minReadLength){r2.setDiscarded(true);}
					}

					if((removePairsIfEitherBad && (r1.discarded() || (r2!=null && r2.discarded()))) || 
							(r1.discarded() && (r2==null || r2.discarded()))){
						if(r1!=null){
							lowqBasesT+=r1.bases.length;
							lowqReadsT++;
						}
						if(r2!=null){
							lowqBasesT+=r2.bases.length;
							lowqReadsT++;
						}
					}else{
						//Process kmers

						if(r1!=null && !r1.discarded()){
							long temp=addKmersToTable(r1);
							added+=temp;
							if(verbose){System.err.println("Added "+temp);}
						}
						if(r2!=null && !r2.discarded()){
							long temp=addKmersToTable(r2);
							added+=temp;
							if(verbose){System.err.println("Added "+temp);}
						}
					}
				}
				
				//Fetch a new read list
				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln, ln.list.isEmpty());
			added+=table.flush();
		}
		
		

		
		/**
		 * Counts the number of kmer hits for a read.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return Number of hits
		 */
		private final int addKmersToTable(final Read r){
			if(r==null || r.bases==null){return 0;}
			final byte[] bases=r.bases;
			final long shift=2*k;
			final long shift2=shift-2;
			final long mask=~((-1L)<<shift);
			long kmer=0;
			long rkmer=0;
			int created=0;
			int len=0;

			if(bases==null || bases.length<k){return -1;}
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					final long key=toValue(kmer, rkmer);
					if(!prefilter || prefilterArray.read(key)>filterMax){
						int temp=table.incrementAndReturnNumCreated(key);
						created+=temp;
						if(verbose){System.err.println("Added "+temp);}
					}
				}
			}
			return created;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadStreamInterface cris;
		
		private final HashBuffer table;
		
		public long added=0;
		
		private long readsInT=0;
		private long basesInT=0;
		private long readsOutT=0;
		private long basesOutT=0;
		private long readsTrimmedT=0;
		private long basesTrimmedT=0;
		private long lowqReadsT=0;
		private long lowqBasesT=0;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public boolean dumpKmersAsText(String fname, int k){
		if(fname==null){return true;}
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, true);
		tsw.start();
		for(AbstractKmerTable set : keySets){
			set.dumpKmersAsText(tsw, k);
		}
		tsw.poisonAndWait();
		return tsw.errorState;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
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
	private final long toValue(long kmer, long rkmer){
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Has this class encountered errors while processing? */
	public boolean errorState=false;

	/** Use a count-min prefilter for low-depth kmers */
	public boolean prefilter=false;
	
	/** Initial size of data structures */
	private int initialSize=128000;
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private final AbstractKmerTable[] keySets;
	/** A scaffold's name is stored at scaffoldNames.get(id).  
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	private final ArrayList<String> scaffoldNames=new ArrayList<String>();
	
	private KCountArray prefilterArray=null;
	
	/** Input reads */
	private String in1=null, in2=null;
	/** Output reads */
	private String out1=null, out2=null;
	/** Statistics output file */
	private String outstats=null;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	private long maxReads=-1;
	/** Output reads in input order.  May reduce speed. */
	private boolean ORDERED=false;
	
	private int buflen=1000;
	
	long readsIn=0;
	long basesIn=0;
	long readsOut=0;
	long basesOut=0;
	long readsTrimmed=0;
	long basesTrimmed=0;
	long lowqReads=0;
	long lowqBases=0;
	
//	long refReads=0;
//	long refBases=0;
//	long refKmers=0;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of tables (and threads, during loading) */ 
	private final int WAYS;
	
	/** Filter kmers up to this level; don't store them in primary data structure */ 
	private final int filterMax;
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
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
	
	public static int VERSION=1;
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Verbose messages */
	public static final boolean verbose=false;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.THREADS;
	/** Indicates end of input stream */
	private static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	/** Do garbage collection prior to printing memory usage */
	private static final boolean GC_BEFORE_PRINT_MEMORY=false;
	
	/*--------------------------------------------------------------*/
	/*----------------      Static Initializers     ----------------*/
	/*--------------------------------------------------------------*/
	
//	static{
//	}
	

	
}
