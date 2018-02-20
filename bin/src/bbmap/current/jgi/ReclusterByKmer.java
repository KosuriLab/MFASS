package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.RTextOutputStream3;
import stream.Read;
import align2.ListNum;
import align2.Shared;
import align2.Tools;
import dna.AminoAcid;
import dna.Data;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;

/**
 * @author Brian Bushnell
 * @date Feb 7, 2014
 *
 */
public class ReclusterByKmer {
	
	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		ReclusterByKmer rbk=new ReclusterByKmer(args);
		rbk.process(t);
	}
	
	public ReclusterByKmer(String[] args){
		if(args==null || args.length==0){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether it was explicitly set.

		FastaReadInputStream.SPLIT_READS=false;
		stream.FastaReadInputStream.MIN_READ_LEN=1;
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.READ_BUFFER_NUM_BUFFERS=Tools.min(8, Shared.READ_BUFFER_NUM_BUFFERS);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		int k1_=12;
		int k2_=3;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("null") || a.equals(in2)){
				// do nothing
			}else if(a.equals("passes")){
				assert(false) : "'passes' is disabled.";
//				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
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
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("t") || a.equals("threads")){
				Shared.THREADS=Tools.max(Integer.parseInt(b), 1);
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("bf1")){
				ByteFile.FORCE_MODE_BF1=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF2=!ByteFile.FORCE_MODE_BF1;
			}else if(a.equals("bf2")){
				ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF1=!ByteFile.FORCE_MODE_BF2;
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				out1=b;
			}else if(a.equals("out2") || a.equals("output2")){
				out2=b;
			}else if(a.equals("extin")){
				extin=b;
			}else if(a.equals("extout")){
				extout=b;
			}else if(a.equals("trd") || a.equals("trc") || a.equals("trimreaddescription")){
				Shared.TRIM_READ_COMMENTS=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("fastareadlen") || a.equals("fastareadlength")){
				FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
				FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
			}else if(a.equals("fastaminread") || a.equals("fastaminlen") || a.equals("fastaminlength")){
				FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
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
			}else if(a.equals("tuc") || a.equals("touppercase")){
				Read.TO_UPPER_CASE=Tools.parseBoolean(b);
			}else if(a.equals("testinterleaved")){
				FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
				outstream.println("Set TEST_INTERLEAVED to "+FASTQ.TEST_INTERLEAVED);
				setInterleaved=true;
			}else if(a.equals("forceinterleaved")){
				FASTQ.FORCE_INTERLEAVED=Tools.parseBoolean(b);
				outstream.println("Set FORCE_INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				setInterleaved=true;
			}else if(a.equals("interleaved") || a.equals("int")){
				if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
				else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
					setInterleaved=true;
				}
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(out1==null && i==1 && !arg.contains("=")){
				out1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.THREADS>2){
//			if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null){
			if(out2!=null){
				printOptions();
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
			//out1="stdout";
			System.err.println("Warning: output destination not set; producing no output.  To print to standard out, set 'out=stdout.fq'");
		}
		
		if(!setInterleaved){
			assert(in1!=null && out1!=null) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, false, out1, out2)){
			throw new RuntimeException("\n\nOVERWRITE="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		
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
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

		/* Check for output file collisions */
		Tools.testOutputFiles(overwrite, false, out1, out2);

		k1=k1_;
		assert(k1>=-1 && k1<=15) : "k1 must lie between 1 and 15, inclusive (0 to disable)";
		k2=k2_;
		assert(k2>=-1 && k2<=6) : "k2 must lie between 1 and 6, inclusive (0 to disable)";

		arraylen1=(k1>0 ? maxCanonicalKmer(k1)+1 : 0);
		arraylen2=(k2>0 ? maxCanonicalKmer(k2)+1 : 0);
	}
	
	/*--------------------------------------------------------------*/
	
	/** TODO */
	public static void printOptions(){
		System.err.println("Usage information unavailable");
	}
	
	void process(Timer t){
		
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, colorspace, false, ffin1, ffin2, null, null);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}

		RTextOutputStream3 ros=null;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=new RTextOutputStream3(ffout1, ffout2, null, null, buff, null, false);
			ros.start();
		}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					{
						readsProcessed++;
						basesProcessed+=r1.bases==null ? 0 : r1.bases.length;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=r2.bases==null ? 0 : r2.bases.length;
					}
					
					
					boolean remove=false;
					if(remove){reads.set(idx, null);}
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();

		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private Cluster fetchCluster(int x){
		if(x>clusterList.size()){
			synchronized(clusterList){
				clusterList.ensureCapacity(2*x);
				for(int i=clusterList.size(); i<x; i++){
					clusterList.add(new Cluster(i));
				}
				clusterList.notifyAll();
			}
		}
		Cluster c=clusterList.get(x);
		while(c==null){
			synchronized(clusterList){
				c=clusterList.get(x);
				assert(c!=null);
			}
		}
		return c;
	}
	
	/**
	 * Generate kmer spectra for clusters
	 * @param t
	 */
	private void findKmerSpectra(Timer t){
		
		/* Create read input stream */
		final ConcurrentReadStreamInterface cris;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, false, ffin1, ffin2);
			Thread cristhread=new Thread(cris);
			cristhread.start();
			if(verbose){System.err.println("Started cris");}
		}
		
		/* Create ClusterThreads */
		ArrayList<ClusterThread> alct=new ArrayList<ClusterThread>(THREADS);
		
		for(int i=0; i<THREADS; i++){alct.add(new ClusterThread(i, CLUSTER_MODE_CREATE, -1, cris));}
		for(ClusterThread ct : alct){ct.start();}

		long readsIn=0, basesIn=0;
		
		/* Wait for threads to die, and gather statistics */
		for(ClusterThread ct : alct){
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					ct.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			readsIn+=ct.readsInT;
			basesIn+=ct.basesInT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris);
	}
	
	/**
	 * Assign reads to clusters using additional kmer information.
	 * @param t
	 */
	private void recluster(Timer t){
		
		/* Create read input stream */
		final ConcurrentReadStreamInterface cris;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, false, ffin1, ffin2);
			Thread cristhread=new Thread(cris);
			cristhread.start();
			if(verbose){System.err.println("Started cris");}
		}
		
		/* Create ClusterThreads */
		ArrayList<ClusterThread> alct=new ArrayList<ClusterThread>(THREADS);
		for(int i=0; i<THREADS; i++){alct.add(new ClusterThread(i, CLUSTER_MODE_RECLUSTER, ambigMode, cris));}
		for(ClusterThread ct : alct){ct.start();}
		
		long readsIn=0, basesIn=0;
		
		/* Wait for threads to die, and gather statistics */
		for(ClusterThread ct : alct){
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					ct.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			readsIn+=ct.readsInT;
			basesIn+=ct.basesInT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris);
	}
	
	/**
	 * @param bases
	 * @param object
	 * @param k
	 * @return
	 */
	public static int[] toKmerCounts(byte[] bases, Object object, int k) {
		// TODO Auto-generated method stub
		return null;
	}

	public static int[] toKmers(final byte[] bases, int[] array_, final int k){
		if(bases==null || bases.length<k){return null;}
		final int alen=bases.length-k+1;
		final int[] array=(array_!=null && array_.length==alen ? array_ : new int[alen]);
		
		final int shift=2*k;
		final int shift2=shift-2;
		final int mask=~((-1)<<shift);
		
		int kmer=0;
		int rkmer=0;
		int len=0;
		
		for(int i=0, j=0; i<bases.length; i++){
			byte b=bases[i];
			int x=Dedupe.baseToNumber[b];
			int x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
//			if(b=='N'){len=0;}else{len++;} //This version will transform 'N' into 'A'
			if(verbose){System.err.println("Scanning2 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k), Tools.min(i+1, k)));}
			if(len>=k){
				array[j]=Tools.min(kmer, rkmer);
				j++;
			}
		}
		
		Arrays.sort(array);
		return array;
	}
	
	public static int[] toKmerCounts(final byte[] bases, int[] array_, final int k, final int alen){
		if(bases==null || bases.length<k){return null;}
		final int[] array=(array_!=null && array_.length==alen ? array_ : new int[alen]);
		
		final int shift=2*k;
		final int shift2=shift-2;
		final int mask=~((-1)<<shift);
		
		int kmer=0;
		int rkmer=0;
		int len=0;
		
		for(int i=0, j=0; i<bases.length; i++){
			byte b=bases[i];
			int x=Dedupe.baseToNumber[b];
			int x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
//			if(b=='N'){len=0;}else{len++;} //This version will transform 'N' into 'A'
			if(verbose){System.err.println("Scanning2 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k), Tools.min(i+1, k)));}
			if(len>=k){
				array[Tools.min(kmer, rkmer)]++;
			}
		}
		
		Arrays.sort(array);
		return array;
	}
	
	public static int maxCanonicalKmer(int k){
		final int bits=2*k;
		final int max=(int)((1L<<bits)-1);
		int high=0;
		for(int kmer=0; kmer<=max; kmer++){
			int canon=Tools.min(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k));
			high=Tools.max(canon, high);
		}
		return high;
	}
	
	/**
	 * @param kmers Read kmers
	 * @param counts Cluster kmer counts
	 * @return Score
	 */
	private static final float andCount(int[] kmers, AtomicLongArray counts){
		int sum=0;
		for(int i=0; i<kmers.length; i++){
			int kmer=kmers[i];
			long count=counts.get(kmer);
			if(count>0){sum++;}
		}
		return sum/(float)kmers.length;
	}
	
	/**
	 * @param kmers Read kmers
	 * @param probs Cluster kmer frequencies
	 * @return Score
	 */
	private static final float innerProduct(int[] kmers, float[] probs){
		float sum=0;
		for(int kmer : kmers){
			if(kmer>=0){
				sum+=probs[kmer];
			}
		}
		return sum;
	}
	
	/**
	 * @param a Read kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	private static final float absDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			sum+=Tools.absdif((double)a[i], (double)b[i]);
		}

		return (float)sum;
	}
	
	/**
	 * @param a Read kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	private static final float rmsDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			double d=Tools.absdif((double)a[i], (double)b[i]);
			sum+=d*d;
		}

		return (float)Math.sqrt(sum/a.length);
	}
	
	/**
	 * @param a Read kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	private static final float ksFunction(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			double d=(double)a[i]*Math.log(a[i]/b[i]);
			sum+=d;
		}
		
		return (float)sum;
	}
	
	/*--------------------------------------------------------------*/
	
	private class ClusterThread extends Thread{
		
		public ClusterThread(int id_, int clusterMode_, int ambigMode_, ConcurrentReadStreamInterface cris_){
			id=id_;
			ambigMode=ambigMode_; 
			clusterMode=clusterMode_;
			cris=cris_;
			
			randy=(ambigMode==AMBIG_MODE_RAND) ? ThreadLocalRandom.current() : null;
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(reads!=null && reads.size()>0){
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					processRead(reads.get(i));
				}
				
				//Fetch a new read list
				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln, ln.list.isEmpty());
		}
		
		
		private void processRead(final Read r1){
			final Read r2=r1.mate;

			if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}

			readsInT++;
			basesInT+=(r1.bases==null ? 0 : r1.bases.length);
			if(r2!=null){
				readsInT++;
				basesInT+=(r2.bases==null ? 0 : r2.bases.length);
			}

			final ReadTag rt1=new ReadTag(r1);
			final ReadTag rt2=(r2==null ? null : new ReadTag(r2));

			r1.obj=rt1;
			if(r2!=null){r2.obj=rt2;}

			if(clusterMode==CLUSTER_MODE_CREATE){
				addToCluster(r1, r2, rt1, rt2);
			}else if(clusterMode==CLUSTER_MODE_RECLUSTER){
				reCluster(r1, r2, rt1, rt2);
			}else{
				throw new RuntimeException("Unknown mode "+clusterMode);
			}
		}
		
		private void addToCluster(Read r1, Read r2, ReadTag rt1, ReadTag rt2){
			final int cn1=rt1.cluster0;
			final int cn2=rt2==null ? cn1 : rt2.cluster0;
			if(cn1==cn2){
				Cluster c1=fetchCluster(cn1);
				c1.add(r1);
				c1.add(r2);
			}else{
				Cluster c1=fetchCluster(cn1);
				Cluster c2=fetchCluster(cn2);
				c1.add(r1);
				c1.add(r2);
				c2.add(r1);
				c2.add(r2);
			}
		}
		
		private void reCluster(Read r1, Read r2, ReadTag rt1, ReadTag rt2){

			assert(false) : "TODO";

			Cluster bestCluster1=null;
			Cluster bestCluster2=null;

			float bestScore1=-999999999, bestScore1_2=-999999999;
			float bestScore2=-999999999, bestScore2_1=-999999999;

			for(Cluster c : clusterList){

				float score1=c.score(r1);
				float score2=c.score(r2);

				if(bestCluster1==null || score1>bestScore1){
					bestCluster1=c;
					bestScore1=score1;
					bestScore1_2=score2;
				}
				if(bestCluster2==null || score2>bestScore2){
					bestCluster2=c;
					bestScore2=score2;
					bestScore2_1=score1;
				}
			}

			if(r2==null){
				rt1.cluster1=bestCluster1.id;
			}else if(bestCluster1==bestCluster2){
				rt1.cluster1=rt2.cluster1=bestCluster1.id;
			}else{
				assert(r1!=null && r2!=null && bestCluster1!=bestCluster2);

				float a=bestScore1+bestScore1_2;
				float b=bestScore2+bestScore2_1;

				if(ambigMode==AMBIG_MODE_BEST){
					if(a>=b){
						rt1.cluster1=rt2.cluster1=bestCluster1.id;
					}else{
						rt1.cluster1=rt2.cluster1=bestCluster2.id;
					}
				}else if(ambigMode==AMBIG_MODE_BOTH){
					assert(false) : "TODO";
				}else if(ambigMode==AMBIG_MODE_TOSS){
					rt1.cluster1=rt2.cluster1=-1;
				}else if(ambigMode==AMBIG_MODE_RAND){
					if(a<0 || b<0){
						float c=0-(Tools.min(a, b))*1.5f;
						a=a+c;
						b=a+c;
					}
					float coin=randy.nextFloat()*(a+b);
					if(coin<=a){
						rt1.cluster1=rt2.cluster1=bestCluster1.id;
					}else{
						rt1.cluster1=rt2.cluster1=bestCluster2.id;
					}
				}
			}

		}

		final int id;
		final int clusterMode;
		final int ambigMode;
		final ConcurrentReadStreamInterface cris;
		
		final ThreadLocalRandom randy;
		
		long readsInT;
		long basesInT;
		
	}
	
	/*--------------------------------------------------------------*/
	
	private class ReadTag{
		
		public ReadTag(Read r_){
			r=r_;
			strand=r.strand();

			int gcCount_=0;
			for(byte b : r.bases){
				if(b=='G' || b=='C'){
					gcCount_++;
				}
			}
			gcCount=gcCount_;
			
			processHeader(r.id);
		}
		
		private void processHeader(String s){
			assert(false) : "TODO";
			gc=-1;
			depth=-1;
			cluster0=-1;
		}

		Read r1(){
			return strand==0 ? r : r.mate;
		}
		
		Read r2(){
			return strand==1 ? r : r.mate;
		}
		
		ReadTag tag1(){
			return (ReadTag)r1().obj;
		}
		
		ReadTag tag2(){
			Read r2=r2();
			return r2==null ? null : (ReadTag)r2.obj;
		}
		
//		private int[] toKmers(final int k){
//			return ReclusterByKmer.toKmers(r.bases, null, k);
//		}
		
		int[] kmerArray1(){
			if(kmerArray1==null){kmerArray1=ReclusterByKmer.toKmers(r.bases, null, k1);}
			return kmerArray1;
		}
		
		int[] kmerArray2(){
			if(kmerArray2==null){kmerArray2=ReclusterByKmer.toKmerCounts(r.bases, null, k2);}
			return kmerArray2;
		}
		
		float[] kmerFreq2(){
			if(kmerFreq2==null){
				int[] counts=kmerArray2();
				if(counts!=null){
					long sum=Tools.sum(counts);
					kmerFreq2=new float[counts.length];
					float extra=(0.05f/counts.length);
					float mult=0.95f/sum;
					for(int i=0; i<counts.length; i++){
						kmerFreq2[i]=counts[i]*mult+extra;
					}
				}
			}
			return kmerFreq2;
		}
		
		/** Sorted long kmers */
		private int[] kmerArray1;
		
		/** Canonically-ordered short kmer counts */ 
		private int[] kmerArray2;
		
		private float[] kmerFreq2;
		
		final Read r;
		final byte strand;
		final int gcCount;
		
		int depth;
		int cluster0=-1; //initial cluster
		int cluster1=-1; //final cluster

		float gc;
		
	}
	
	/*--------------------------------------------------------------*/
	
	private class Cluster{
		
		public Cluster(int id_){
			id=id_;
			
			kmerArray1=new AtomicLongArray(arraylen1);
			kmerProbArray1=new float[arraylen1];

			kmerArray2=new AtomicLongArray(arraylen2);
			kmerProbArray2=new float[arraylen2];
		}

		public void recalculate(){
			gc=(float)(gcCount.doubleValue()/baseCount.doubleValue());

			if(k1>0){
				long kmerCount=0;
				for(int i=0; i<arraylen1; i++){
					kmerCount+=kmerArray1.get(i);
				}
				double extra=(0.05/arraylen1);
				double mult=(0.95/kmerCount);
				for(int i=0; i<arraylen1; i++){
					kmerProbArray1[i]=(float)(kmerArray1.get(i)*mult+extra);
				}
			}
			if(k2>0){
				long kmerCount=0;
				for(int i=0; i<arraylen2; i++){
					kmerCount+=kmerArray2.get(i);
				}
				double extra=(0.05/arraylen2);
				double mult=(0.95/kmerCount);
				for(int i=0; i<arraylen2; i++){
					kmerProbArray2[i]=(float)(kmerArray2.get(i)*mult+extra);
				}
			}
		}
		
		public void resetAtomics(){
			for(int i=0; i<arraylen1; i++){
				kmerArray1.set(i, 0);
			}
			for(int i=0; i<arraylen2; i++){
				kmerArray2.set(i, 0);
			}
			depthsum1.set(0);
			depthsum2.set(0);
			readCount.set(0);
			baseCount.set(0);
			gcCount.set(0);
		}
		
		public void add(Read r){
			if(r==null){return;}
			ReadTag rt=(ReadTag)r.obj;
			assert(rt!=null);
			final byte[] bases=r.bases;
			
			readCount.incrementAndGet();
			baseCount.addAndGet(bases.length);
			gcCount.addAndGet(rt.gcCount);
			
			if(rt.strand==0){
				depthsum1.addAndGet(rt.depth);
			}else{
				depthsum2.addAndGet(rt.depth);
			}
			
			if(k1>0){
				int[] kmers=rt.kmerArray1();
				int kmer=-1, run=0;
				for(int i=0; i<kmers.length; i++){
					int x=kmers[i];
					if(x==kmer){
						run++;
					}else{
						if(run>0){kmerArray1.addAndGet(kmer, run);}
						kmer=x;
						run=1;
					}
				}
				if(run>0){kmerArray1.addAndGet(kmer, run);}
			}

			if(k2>0){
				int[] kmers=rt.kmerArray2();
				for(int kmer=0; kmer<kmers.length; kmer++){
					int x=kmers[kmer];
					if(x>0){kmerArray2.addAndGet(kmer, x);}
				}
			}
		}
		
		/**
		 * @param r1
		 * @return
		 */
		public float score(Read r) {
			if(r==null){return 0;}
			return r.mate==null ? scoreSingle(r) : scorePaired(r);
		}
		
		/**
		 * @param r1
		 * @return
		 */
		public float scoreSingle(Read r) {
			if(r==null){return 0;}
			ReadTag rt=(ReadTag)r.obj;
			
			assert(false) : "TODO";
			float depthScore=scoreDepthSingle(rt);
			float gcScore=scoreGcSingle(rt);
			float kmerScore=scoreKmer1(rt);
			assert(false);
			float depthWeight=.2f;
			float gcWeight=.2f;
			float kmerWeight=.6f;
			
			return depthWeight*depthScore+gcWeight*gcScore+kmerWeight*kmerScore;
		}
		
		/**
		 * @param rt
		 * @return
		 */
		private float scoreKmer1(ReadTag rt) {
			int[] kmers=rt.kmerArray1();
			
			float score=0;
			if(scoreMode1==SCORE_MODE_AND){
				float f=andCount(kmers, kmerArray1);
				assert(false);
			}else if(scoreMode1==SCORE_MODE_MULT){
				float f=innerProduct(kmers, kmerProbArray1);
				assert(false);
			}else{
				throw new RuntimeException(""+scoreMode1);
			}
			
			return score;
		}
		
		/**
		 * @param rt
		 * @return
		 */
		private float scoreKmer2(ReadTag rt) {
			int[] kmers=rt.kmerArray2();
			float[] probs=rt.kmerFreq2();
			
			float score=0;
			if(scoreMode2==SCORE_MODE_AND){
				float f=andCount(kmers, kmerArray2);
				assert(false);
			}else if(scoreMode2==SCORE_MODE_MULT){
				float f=innerProduct(kmers, kmerProbArray2);
				assert(false);
			}if(scoreMode2==SCORE_MODE_DIF){
				float f=absDif(probs, kmerProbArray2);
				assert(false);
			}else if(scoreMode2==SCORE_MODE_RMS){
				float f=rmsDif(probs, kmerProbArray2);
				assert(false);
			}else if(scoreMode2==SCORE_MODE_KS){
				float f=ksFunction(probs, kmerProbArray2);
				assert(false);
			}else{
				throw new RuntimeException(""+scoreMode2);
			}
			
			return score;
		}

		/**
		 * @param rt
		 * @return
		 */
		private float scoreGcSingle(ReadTag rt) {
			assert(false) : "TODO";
			// TODO Auto-generated method stub
			return 0;
		}

		/**
		 * @param rt
		 * @return
		 */
		private float scoreDepthSingle(ReadTag rt) {
			assert(false) : "TODO";
			// TODO Auto-generated method stub
			return 0;
		}

		/**
		 * @param r1
		 * @return
		 */
		public float scorePaired(Read r) {
			assert(false) : "TODO";
			if(r==null){return 0;}
			ReadTag rt=(ReadTag)r.obj;
			
//			ReadTag rt1=rt.r
			
			return 0;
		}
		
		public final int id;
		
		public float gc;
		public int depth1, depth2;
		
		final AtomicLongArray kmerArray1;
		final float[] kmerProbArray1;
		
		final AtomicLongArray kmerArray2;
		final float[] kmerProbArray2;
		
		final AtomicLong depthsum1=new AtomicLong(0);
		final AtomicLong depthsum2=new AtomicLong(0);
		
		final AtomicLong readCount=new AtomicLong(0);
		final AtomicLong baseCount=new AtomicLong(0);
//		final AtomicLong kmerCount=new AtomicLong(0);
		final AtomicLong gcCount=new AtomicLong(0);
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	public boolean errorState=false;
	
	private final ArrayList<Cluster> clusterList=new ArrayList<Cluster>(256);
	
	/** 'big' kmer */
	public final int k1;
	/** 'small' kmer */
	public final int k2;

	public final int arraylen1;
	public final int arraylen2;
	
	private String in1=null;
	private String in2=null;

	private String out1=null;
	private String out2=null;
	
	private String extin=null;
	private String extout=null;
	
	private boolean overwrite=false;
	private boolean colorspace=false;
	
	private long maxReads=-1;
	
	private byte qin=-1;
	private byte qout=-1;

	private int scoreMode1=SCORE_MODE_MULT;
	private int scoreMode2=SCORE_MODE_RMS;
	private int ambigMode=AMBIG_MODE_RAND;
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	
	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	private PrintStream outstream=System.err;
	
	private int THREADS=Shared.THREADS;
	
	/*--------------------------------------------------------------*/

	public static boolean verbose=false;

	public static final int CLUSTER_MODE_CREATE=0;
	public static final int CLUSTER_MODE_RECLUSTER=1;
	public static final int CLUSTER_MODE_REFINE=2;

	public static final int SCORE_MODE_DIF=0;
	public static final int SCORE_MODE_RMS=1;
	public static final int SCORE_MODE_AND=2;
	public static final int SCORE_MODE_MULT=3;
	public static final int SCORE_MODE_KS=4;
	
	public static final int AMBIG_MODE_BEST=0;
	public static final int AMBIG_MODE_BOTH=1;
	public static final int AMBIG_MODE_TOSS=2;
	public static final int AMBIG_MODE_RAND=3;
	
}
