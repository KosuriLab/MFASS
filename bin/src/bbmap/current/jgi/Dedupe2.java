package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import stream.ConcurrentCollectionReadInputStream;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;

import align2.ListNum;
import align2.LongM;
import align2.ReadLengthComparator;
import align2.Shared;
import align2.Tools;
import dna.AminoAcid;
import dna.Data;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;

/**
 * Uses a hashset for lower memory during exact match detection.
 * This makes it probablistic.
 * @author Brian Bushnell
 * @date Jul 18, 2013
 *
 */
public final class Dedupe2 {
	
	public static void main(String[] args){
		
		if(args==null || args.length==0 || (args.length==1 && 
				(args[0].equalsIgnoreCase("-h") || args[0].equals("-help") || args[0].equals("--help") || args[0].equals("-?") || args[0].equals("?")))){
			printOptions();
			System.exit(0);
		}
		Dedupe2 dd=new Dedupe2(args);
		dd.process();
	}
	
	private static void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("\njava -ea -Xmx106g -cp <path> jgi.Dedupe2 <input file> <output file>");
		outstream.println("\nOptional flags:");
		outstream.println("in=<file>          \tThe 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.");
		outstream.println("out=<file>         \tThe 'out=' flag is needed if the output file is not the second parameter.  'out=stdout' will pipe to standard out.");
		outstream.println("showspeed=t        \tSet to 'f' to suppress display of processing speed.");
		outstream.println("minscaf=0          \tIgnore scaffolds shorter than this.");
		outstream.println("testrc=t           \t(trc) Test reverse-complements as well as normal orientation.");
		outstream.println("testmatch=t        \t(tm) Test for exact matches of scaffolds.");
		outstream.println("testcontainment=t  \t(tc) Test for full containments of scaffolds.");
		outstream.println("storename=t        \t(sn) Store scaffold names (set false to save memory).");
		outstream.println("storequality=t     \t(sq) Store quality values for fastq assemblies (set false to save memory).");
		outstream.println("exact=t            \t(ex) Only allow exact symbol matches.  When false, an 'N' will match any symbol.");
		outstream.println("uniquenames=t      \t(un) Ensure all output scaffolds have unique names.  Uses more memory.");
		outstream.println("maxedits=0         \t(e) Absorb contained sequences with up to this many mismatches (subs only, no indels).");
		outstream.println("minidentity=100    \t(mid) Absorb contained sequences with percent identity of at least this (subs only, no indels).");
		outstream.println("k=31               \tKmer length used for finding containments.  Containments shorter than k will not be found.");
		outstream.println("ziplevel=2         \tSet to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.");
		outstream.println("sort=f             \tsort output by scaffold length (otherwise it will be random).\n" +
						  "                   \t'a' for ascending, 'd' for descending, 'f' for false (no sorting).");
		outstream.println("minlengthpercent=0 \t(mlp) Smaller contig must be at least this percent of larger contig's length to be absorbed." +
				  		  "                   \tThis option is only for compatibility with vmatch and changing it is not recommended, " +
				  		  "                   \tas it may cause nondeterministic output.");
	}
	
	public Dedupe2(String[] args){
		for(String s : args){if(s.contains("standardout") || s.contains("stdout")){outstream=System.err;}}
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		ReadWrite.ZIPLEVEL=2;
		//ReadWrite.USE_UNPIGZ=true;
		FastaReadInputStream.SPLIT_READS=false;
		boolean setOut=false;

		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("in")){
				if(b.indexOf(',')>=0 && !new File(b).exists()){
					in=b.split(",");
				}else{
					in=new String[] {b};
				}
			}else if(a.equals("out")){
				out=b;
				setOut=true;
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
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
			}else if(a.equals("sort")){
				if(b==null){sort=true;}
				else if(b.equalsIgnoreCase("a")){
					sort=true;
					ascending=true;
				}else if(b.equalsIgnoreCase("d")){
					sort=true;
					ascending=false;
				}else{
					sort=Tools.parseBoolean(b);
				}
			}else if(a.equals("trc") || a.equals("testrc")){
				ignoreReverseComplement=!Tools.parseBoolean(b);
			}else if(a.equals("tc") || a.equals("testcontainment") || a.equals("containment")){
				testContainment=Tools.parseBoolean(b);
			}else if(a.equals("tm") || a.equals("testmatch")){
				testMatch=Tools.parseBoolean(b);
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				k2=k-1;
				assert(k>0 && k<32) : "k must be between 1 and 31; default is 31, and lower values are slower.";
			}else if(a.equals("minscaf")){
				MINSCAF=FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(a.equals("mlp") || a.equals("minlengthpercent")){
				minLengthPercent=Float.parseFloat(b);
			}else if(a.equals("e") || a.equals("maxedits")){
				maxEdits=Integer.parseInt(b);
			}else if(a.equals("mid") || a.equals("minidentity")){
				minIdentity=Float.parseFloat(b);
				minIdentityMult=(minIdentity==100f ? 0 : (100f-minIdentity)/100f);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Integer.parseInt(b);
			}else if(a.equals("showspeed")){
				showSpeed=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("contigbreak") || (arg.contains("=") && (a.equals("n") || a.equals("-n")))){
				maxNs=Integer.parseInt(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("storename") || a.equals("sn")){
				storeName=Tools.parseBoolean(b);
			}else if(a.equals("storesuffix") || a.equals("ss")){
				storeSuffix=Tools.parseBoolean(b);
			}else if(a.equals("storequality") || a.equals("sq")){
				storeQuality=Tools.parseBoolean(b);
			}else if(a.equals("exact") || a.equals("ex")){
				exact=Tools.parseBoolean(b);
			}else if(a.equals("uniquenames") || a.equals("un")){
				uniqueNames=Tools.parseBoolean(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
			}else if(i==0 && in==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				String c=args[i];
				if(c.indexOf(',')>=0 && !new File(c).exists()){
					in=c.split(",");
				}else{
					in=new String[] {c};
				}
			}else if(i==1 && out==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				out=args[i];
				setOut=true;
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(maxEdits>0 || minIdentity<100){storeSuffix=true;}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		for(int i=0; i<in.length; i++){
			if(in[i].equalsIgnoreCase("stdin") && !new File(in[i]).exists()){in[i]="stdin.fa";}
		}
		
//		assert(false) : Arrays.toString(in);
		
		if(!setOut){out="stdout.fa";}
		else if("stdout".equalsIgnoreCase(out) || "standarddout".equalsIgnoreCase(out)){
			out="stdout.fa";
			outstream=System.err;
		}
		if(!Tools.canWrite(out, overwrite)){throw new RuntimeException("Output file "+out+" already exists, and overwrite="+overwrite);}
		
		for(int i=0; i<in.length; i++){
			assert(!in[i].equalsIgnoreCase(out));
		}
		if(testContainment){affixMap=new HashMap<LongM, ArrayList<Unit>>(4000000);}
//		assert(false) : testContainment+", "+(affixMap==null);
	}
	
	public void process(){
		
		Timer t=new Timer();
		t.start();
		
		boolean dq0=FASTQ.DETECT_QUALITY;
		boolean ti0=FASTQ.TEST_INTERLEAVED;
		int rbl0=Shared.READ_BUFFER_LENGTH;
		FASTQ.DETECT_QUALITY=false;
		FASTQ.TEST_INTERLEAVED=false;
		Shared.READ_BUFFER_LENGTH=16;
		
		process2();
		
		FASTQ.DETECT_QUALITY=dq0;
		FASTQ.TEST_INTERLEAVED=ti0;
		Shared.READ_BUFFER_LENGTH=rbl0;
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:   \t\t\t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException("Dedupe2 terminated in an error state; the output may be corrupt.");
		}
	}
	
	private static final void printMemory(){
		long mmemory=Runtime.getRuntime().maxMemory()/1000000;
		long tmemory=Runtime.getRuntime().totalMemory()/1000000;
		long fmemory=Runtime.getRuntime().freeMemory()/1000000;
		long umemory=tmemory-fmemory;
		outstream.println("Memory: max="+mmemory+", total="+tmemory+", free="+fmemory+", used="+umemory);
	}
	
	public void process2(){
		
		final TextStreamWriter tsw=(out==null ? null : new TextStreamWriter(out, overwrite, false, true));
		
		Timer t=new Timer();
		t.start();
		outstream.println("Initial:");
		printMemory();
		outstream.println();
		
		{
			crisa=new ConcurrentGenericReadInputStream[in.length];
			multipleInputFiles=crisa.length>1;
			for(int i=0; i<in.length; i++){
				final ConcurrentReadStreamInterface cris;
				{
					FileFormat ff1=FileFormat.testInput(in[i], FileFormat.FASTQ, null, !multipleInputFiles, true);
					cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, null);
					Thread th=new Thread(cris);
					th.start();
				}
				crisa[i]=cris;
			}

			ArrayList<HashThread2> alht=new ArrayList<HashThread2>(THREADS);
			for(int i=0; i<THREADS; i++){alht.add(new HashThread2(true, testContainment, false, testMatch));}
			for(HashThread2 ht : alht){ht.start();}
			for(HashThread2 ht : alht){
				while(ht.getState()!=Thread.State.TERMINATED){
					try {
						ht.join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				matches+=ht.matchesT;
				collisions+=ht.collisionsT;
				containments+=ht.containmentsT;
				containmentCollisions+=ht.containmentCollisionsT;
				basecontainments+=ht.basecontainmentsT;
				basematches+=ht.basematchesT;
				addedToMain+=ht.addedToMainT;
				readsProcessed+=ht.readsProcessedT;
				basesProcessed+=ht.basesProcessedT;
			}
			alht.clear();
			for(ConcurrentReadStreamInterface cris : crisa){
				errorState|=ReadWrite.closeStream(cris);
			}
			crisa=null;
		}
		synchronized(this){
			t.stop();
			outstream.println("Finished exact matches.  Time: "+t);
			printMemory();
			outstream.println();
			Tools.pause(800);
			t.start();
		}
		if(testContainment){
			ArrayList<Read> list=new ArrayList<Read>((int)addedToMain);
			for(Unit u : codeMap){
					if(u.valid()){list.add(u.r);}
			}
			
//			if(minLengthPercent>0){
//				if(verbose){System.err.println("Sorting.");}
//				Collections.sort(list, ReadLengthComparator.comparator);
//				Collections.reverse(list);
//				assert(list.isEmpty() || list.get(0).bases.length<=list.get(list.size()-1).bases.length) : 
//					list.get(0).bases.length+", "+list.get(list.size()-1).bases.length;
//			}
			
			synchronized(this){
				t.stop();
				outstream.println("Allocated list.  Time: "+t);
				printMemory();
				outstream.println();
				Tools.pause(800);
				t.start();
			}
			crisa=new ConcurrentCollectionReadInputStream[] {new ConcurrentCollectionReadInputStream(list, null, maxReads)};
			Thread cristhread=new Thread(crisa[0]);
			cristhread.start();
			
			ArrayList<HashThread2> alht=new ArrayList<HashThread2>(THREADS);
			for(int i=0; i<THREADS; i++){alht.add(new HashThread2(false, false, true, false));}
			
			for(HashThread2 ht : alht){ht.start();}
			for(HashThread2 ht : alht){
				while(ht.getState()!=Thread.State.TERMINATED){
					try {
						ht.join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
//				matches+=ht.matchesT;
//				collisions+=ht.collisionsT;
				containments+=ht.containmentsT;
				containmentCollisions+=ht.containmentCollisionsT;
				basecontainments+=ht.basecontainmentsT;
//				basematches+=ht.basematchesT;
//				addedToMain+=ht.addedToMainT;
//				readsProcessed+=ht.readsProcessedT;
//				basesProcessed+=ht.basesProcessedT;
			}
			alht.clear();
			for(ConcurrentReadStreamInterface cris : crisa){
				errorState|=ReadWrite.closeStream(cris);
			}
			crisa=null;
			affixMap.clear();
			affixMap=null;
		}
		synchronized(this){
			t.stop();
			outstream.println("Finished containment.  Time: "+t);
			printMemory();
			outstream.println();
			Tools.pause(800);
			t.start();
		}
		
		outstream.println("Input:                  \t"+readsProcessed+" reads \t\t"+basesProcessed+" bases.");
		outstream.println("Duplicates:             \t"+matches+" reads ("+String.format("%.2f",matches*100.0/readsProcessed)+"%) \t"+
				basematches+" bases ("+String.format("%.2f",basematches*100.0/basesProcessed)+"%)     \t"+collisions+" collisions.");
		outstream.println("Containments:           \t"+containments+" reads ("+String.format("%.2f",containments*100.0/readsProcessed)+"%) \t"+
				basecontainments+" bases ("+String.format("%.2f",basecontainments*100.0/basesProcessed)+"%)    \t"+containmentCollisions+" collisions.");
//		outstream.println("Result:                 \t"+(addedToMain-containments)+" reads \t\t"+(basesProcessed-basematches-basecontainments)+" bases.");
		
		long outReads=(addedToMain-containments);
		long outBases=(basesProcessed-basematches-basecontainments);
		outstream.println("Result:                 \t"+outReads+" reads ("+String.format("%.2f",outReads*100.0/readsProcessed)+"%) \t"+
				outBases+" bases ("+String.format("%.2f",outBases*100.0/basesProcessed)+"%)");
		
		outstream.println("");
		
		if(tsw!=null){
			writeOutput(tsw, t);
		}
		
	}
	
	private static ArrayList<Read> addToArray(HashSet<Unit> codeMap, boolean sort, boolean ascending, boolean clear, long outNum){
		assert(outNum<=Integer.MAX_VALUE);
		if(verbose){System.err.println("Making list.");}
		ArrayList<Read> list=new ArrayList<Read>((int)outNum);
		if(verbose){System.err.println("Adding.");}
		for(Unit u : codeMap){
			if(u.valid()){list.add(u.r);}
		}
		if(clear){codeMap.clear();}
		
		if(sort){
			if(verbose){System.err.println("Sorting.");}
			Collections.sort(list, ReadLengthComparator.comparator);
			if(ascending){
				Collections.reverse(list);
				assert(list.isEmpty() || list.get(0).bases.length<=list.get(list.size()-1).bases.length) : 
					list.get(0).bases.length+", "+list.get(list.size()-1).bases.length;
			}else{
				assert(list.isEmpty() || list.get(0).bases.length>=list.get(list.size()-1).bases.length) : 
					list.get(0).bases.length+", "+list.get(list.size()-1).bases.length;
			}
		}
		assert(list.size()==outNum) : list.size()+", "+outNum;
		return list;
	}
	
	private void writeOutput(TextStreamWriter tsw, Timer t){
		
		ArrayList<Read> list=addToArray(codeMap, sort, ascending, true, addedToMain-containments);
		codeMap=null;
		
		if(sort){
			synchronized(this){
				t.stop();
				outstream.println("Sorted output.  Time: "+t);
				printMemory();
				outstream.println();
				Tools.pause(800);
				t.start();
			}
		}
		
		writeOutput(tsw, list);
		
		synchronized(this){
			t.stop();
			outstream.println("Printed output.  Time: "+t);
			printMemory();
			outstream.println();
			Tools.pause(800);
			t.start();
		}
	}
	

	
	private void writeOutput(TextStreamWriter tsw, ArrayList<Read> list){
		
		if(verbose){System.err.println("Writing.");}
		tsw.start();
		
		HashSet<String> names=((uniqueNames && storeName) ? 
				new HashSet<String>(Tools.min(Integer.MAX_VALUE, Tools.max((int)addedToMain, (int)(addedToMain*1.35)))) : null);
		long rid=0;
		for(int x=0; x<list.size(); x++){
			Read r=list.get(x);
			list.set(x, null);
			if(multipleInputFiles){
				r.numericID=(rid++);
			}
			if(names!=null){
				String name=(r.id==null ? ""+r.numericID : r.id);
				if(names.contains(name)){
					for(long i=0; i<Integer.MAX_VALUE; i++){
						String name2=name+"_dd"+i;
						if(!names.contains(name2)){
							r.id=name2;
							names.add(name2);
							break;
						}
					}
				}else{
					names.add(name);
				}
			}
			tsw.println(r);
		}
		tsw.poisonAndWait();
		
	}
	
	
	public static long hash(byte[] bases){
		long code=bases.length;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int mode=(int)(code&31);
			code=code^hashcodes[b][mode];
			code=Long.rotateLeft(code, 1);
		}
		return code;
	}
	
	
	public static long hashReversed(byte[] bases){
		long code=bases.length;
		for(int i=bases.length-1; i>=0; i--){
			byte b=bases[i];
			b=AminoAcid.baseToComplementExtended[b];
			int mode=(int)(code&31);
			code=code^hashcodes[b][mode];
			code=Long.rotateLeft(code, 1);
		}
		return code;
	}
	
	
	public static boolean isCanonical(byte[] bases){
		if(ignoreReverseComplement || bases==null || bases.length==0){return true;}
		final int lim=(bases.length+1)/2;
		for(int i=0, j=bases.length-1; i<lim; i++, j--){
			byte a=bases[i], b=AminoAcid.baseToComplementExtended[bases[j]];
			if(a<b){return true;}
			if(b<a){return false;}
		}
		assert((bases.length&1)==0 || bases[lim]==AminoAcid.baseToComplementExtended[bases[lim]]) : 
			bases.length+", "+bases[lim]+", "+(char)bases[lim]; //palindrome test
		return true; //palindrome
	}
	
	
	private static synchronized long[][] makeCodes(int symbols, int modes){
		Random randy=new Random(1);
		long[][] r=new long[symbols][modes];
		for(int i=0; i<symbols; i++){
			for(int j=0; j<modes; j++){
				r[i][j]=randy.nextLong();
			}
		}
		return r;
	}
	
	/** Handles IUPAC codes */
	private static synchronized long[][] makeCodes2(int modes){
		long[][] r0=makeCodes(26, modes);
		long[][] r=new long[Tools.max('Z','z')+1][];
		for(int i=0; i<26; i++){
			char c=(char)('A'+i);
			r[c]=r[Character.toLowerCase(c)]=r0[i];
		}
		return r;
	}
	
	
	/** 
	 * Creates Unit objects or uses ones already attached to reads.
	 * Places them in local storage and percolates them to shared storage (codeMap), removing exact duplicates.
	 * Also hashes tips and places these in shared affixMap.
	 * Looks for containments in the affix map.
	 * @author Brian Bushnell
	 * @date Jul 24, 2013
	 *
	 */
	private final class HashThread2 extends Thread{
		
		public HashThread2(boolean addToCodeMap_, boolean addToAffixMap_, boolean findContainments_, boolean findMatches_){
			addToCodeMap=addToCodeMap_; 
			addToAffixMap=addToAffixMap_; 
			findContainments=findContainments_;
			findMatches=findMatches_;
			tid=getTid();
//			if(findContainments){containmentMapT=new HashMap<LongM, ArrayList<Unit>>(threadMaxReadsToBuffer*8);}
			crisq=new ArrayDeque<ConcurrentReadStreamInterface>(crisa.length);
			for(int i=0; i<crisa.length; i++){
				crisq.add(crisa[(i+tid)%crisa.length]);
			}
		}
		
		public void run(){
			
			ConcurrentReadStreamInterface cris=crisq.poll();

			while(cris!=null){
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				//			long xx=0;
				while(reads!=null && reads.size()>0){

					for(Read r : reads){
						if(r.bases!=null && r.bases.length>=MINSCAF){
							assert(r.mate==null);
							if(!storeName){r.id=null;}
							if(!storeQuality){r.quality=null;}
							readsProcessedT++;
							//					xx++;
							//					outstream.println("Processing read "+r.id+", "+xx);
							basesProcessedT+=r.bases==null ? 0 : r.bases.length;

//							final long code;
//							final Unit u;
//							if(r.obj==null){
//								final boolean canonical=isCanonical(r.bases);
//								code=(canonical ? hash(r.bases) : hashReversed(r.bases));
//								u=(r.obj!=null ? (Unit)r.obj : new Unit(r, canonical, code));
//								u=(r.obj!=null ? (Unit)r.obj : new Unit(r));
//								r.obj=u;
//							}else{
//								u=(Unit)r.obj;
//								code=u.code;
//							}
//							assert(u.r==r && r.obj==u);

							final Unit u=(r.obj!=null ? (Unit)r.obj : new Unit(r));
							assert(u.r==r && (r.obj==u || r.obj==null));
							final long code=u.code1;
							r.obj=u;
							assert(u.r==r && r.obj==u);

							//					if(verbose){System.err.println("Generated "+code+" for sequence "+new String(r.bases, 0, Tools.min(40, r.bases.length)));}

							if(addToCodeMap){
								boolean b=codeMapT.add(u);
								if(b){
									basesStoredT+=r.bases.length;
								}else{
									matchesT++;
									basematchesT+=r.bases.length;
								}
							}

							if(findContainments){
								//						System.err.println("\naffixMap:\n"+affixMap+"\n");
								int x=findContainments(u);
							}

							//					if(verbose){System.err.println("mapT.size="+mapT.size()+", basesStoredT="+basesStoredT);}
						}
					}

					if(codeMapT!=null && (codeMapT.size()>threadMaxReadsToBuffer || basesStoredT>threadMaxBasesToBuffer)){
						assert(addToCodeMap);
						long added=mergeMaps();
						addedToMainT+=added;
					}

					cris.returnList(ln, ln.list.isEmpty());
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				cris.returnList(ln, ln.list.isEmpty());
				if(codeMapT!=null && !codeMapT.isEmpty()){
					long added=mergeMaps();
					addedToMainT+=added;
				}
				cris=crisq.poll();
			}
			
			codeMapT=null;
		}
		
		private int findContainments(final Unit u){
			if(minLengthPercent<=0 && maxEdits<=0 && minIdentity>=100 && !u.valid()){return 0;}
			final byte[] bases=u.bases();
			final int minlen=k-1;
			final long shift=2*k;
			final long shift2=shift-2;
			final long mask=~((-1L)<<shift);
			long kmer=0;
			long rkmer=0;
			int hits=0;
			int currentContainments=0;
			
			if(bases==null || bases.length<k){return -1;}
			
//			final int stop=bases.length-k;
			final LongM key=new LongM();
			
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=baseToNumber[b];
				long x2=baseToRcompNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(verbose){System.err.println("Scanning i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(i>=minlen){
					key.set(Tools.max(kmer, rkmer)); //Canonical
					ArrayList<Unit> list=affixMap.get(key);
					if(list!=null){
						for(Unit u2 : list){
							if(u!=u2 && !u.equals(u2)){
								if(u2.valid()){
									hits++;
									if(verbose){
//										System.err.println("Found potential containment at i="+i+", key="+key.value()+", pre="+u2.prefix+
//												", suf="+u2.suffix+", tcan="+u2.tipCanonical()+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i, k)));
										System.err.println("Found potential containment at i="+i+", key="+key.value()+", pre="+u2.prefix+
												", suf="+u2.suffix+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i, k)));
									}
									if(u.contains(u2, i, key)){
										synchronized(u2){
											if(u2.valid()){
												currentContainments++;
												basecontainmentsT+=u2.length();
												u2.setInvalid();
											}
										}
										
										if(verbose){System.err.println("Added containment "+u2);}
									}
								}
							}
						}
					}
				}
			}
//			assert(false) : hits+", "+currentContainments+", "+basecontainments+"\n"+containmentMapT+"\n";
			
			containmentCollisionsT+=(hits-currentContainments);
//			outstream.println("hits="+hits+", currentContainments="+currentContainments);
			containmentsT+=currentContainments;
			return hits;
		}
		
		private long mergeMaps(){
			if(verbose){System.err.println("Merging.");}
			long novelReads=0;
			long collisionReads=0;
			
			synchronized(codeMap){
				for(Unit u : codeMapT){
					if(codeMap.contains(u)){
						//do nothing
						matchesT++;
						basematchesT+=u.length();
					}else{
						codeMap.add(u);
						addedList.add(u);
						novelReads++;
					}
				}
			}
			
			if(verbose){System.err.println("Novel reads = "+novelReads+", conflicts = "+0);}
			
			if(verbose){System.err.println("Done Merging.");}
			if(verbose){System.err.println("mapT.size="+codeMapT.size()+", basesStoredT="+basesStoredT);}
			
			codeMapT.clear();
			
			if(!addedList.isEmpty()){
				if(addToAffixMap){
					assert(affixMap!=null);
					synchronized(affixMap){
						for(Unit u : addedList){
							if(u.prefix!=-1 || u.prefix!=u.suffix){
								LongM p=new LongM(u.prefix, false);
								if(affixMap.containsKey(p)){
									affixMap.get(p).add(u);
								}else{
									ArrayList<Unit> alu=new ArrayList<Unit>(2);
									alu.add(u);
									affixMap.put(p, alu);
								}
							}
							if(storeSuffix && u.prefix!=u.suffix){
								LongM p=new LongM(u.suffix, false);
								if(affixMap.containsKey(p)){
									affixMap.get(p).add(u);
								}else{
									ArrayList<Unit> alu=new ArrayList<Unit>(2);
									alu.add(u);
									affixMap.put(p, alu);
								}
							}
						}
					}
				}
			}
			
			addedList.clear();
			basesStoredT=0;
			return collisionReads+novelReads;
		}
		
		private int getTid(){
			synchronized(HashThread2.class){
				int x=tcount;
				tcount++;
				return x;
			}
		}
		
		private HashSet<Unit> codeMapT=new HashSet<Unit>(threadMaxReadsToBuffer*8);
		private ArrayList<Unit> addedList=new ArrayList<Unit>(threadMaxReadsToBuffer);
		
		long matchesT=0;
		long basematchesT=0;
		long basecontainmentsT=0;
		long collisionsT=0;
		long containmentsT=0;
		long containmentCollisionsT=0;
		long basesStoredT=0;
		long addedToMainT=0;
		long readsProcessedT=0;
		long basesProcessedT=0;
		
		private final boolean addToCodeMap;
		private final boolean addToAffixMap;
		private final boolean findContainments;
		private final boolean findMatches;
		private final int tid;
		private final ArrayDeque<ConcurrentReadStreamInterface> crisq;
	}
	
	public static boolean equalsRC(byte[] a, byte[] b){
		if(a==b){return true;}
		if(a==null || b==null){return false;}
		if(a.length!=b.length){return false;}

		boolean ca=isCanonical(a);
		boolean cb=isCanonical(b);
		
		if(ca==cb){
			for(int i=0; i<a.length; i++){
				final byte aa=a[i], bb=b[i];
				if(aa!=bb && aa!='N' && bb!='N'){return false;}
			}
		}else{
			for(int i=0, j=b.length-1; i<a.length; i++, j--){
				final byte aa=a[i], bb=AminoAcid.baseToComplementExtended[b[j]];
				if(aa!=bb && aa!='N' && bb!='N'){return false;}
			}
		}
		return true;
	}
	
	public static boolean equalsRC(Unit ua, Unit ub){
		return ua.code1==ub.code1 && ua.code2==ub.code2 && (ua.canonical()==ub.canonical() ? (ua.prefix==ub.prefix && ua.suffix==ub.suffix) : 
			 (ua.prefix==ub.suffix && ua.suffix==ub.prefix)) && compareRC(ua, ub)==0;
	}
	
	//TODO
	//This is really for sorting by length.
	public static int compareRC(Unit ua, Unit ub){
		if(ua.length()!=ub.length()){return ua.length()>ub.length() ? 1 : -1;}
		if(ua.r==null || ub.r==null){
			if(ua.canonical()){
				if(ub.canonical()){
					if(ua.prefix!=ub.prefix){return ua.prefix>ub.prefix ? 1 : -1;}
					if(ua.suffix!=ub.suffix){return ua.suffix>ub.suffix ? 1 : -1;}
				}else{
					if(ua.prefix!=ub.suffix){return ua.prefix>ub.suffix ? 1 : -1;}
					if(ua.suffix!=ub.prefix){return ua.suffix>ub.prefix ? 1 : -1;}
				}
			}else{
				if(ub.canonical()){
					if(ua.suffix!=ub.prefix){return ua.suffix>ub.prefix ? 1 : -1;}
					if(ua.prefix!=ub.suffix){return ua.prefix>ub.suffix ? 1 : -1;}
				}else{
					if(ua.suffix!=ub.suffix){return ua.suffix>ub.suffix ? 1 : -1;}
					if(ua.prefix!=ub.prefix){return ua.prefix>ub.prefix ? 1 : -1;}
				}
			}
			if(ua.code1!=ub.code1){return ua.code1>ub.code1 ? 1 : -1;}
			if(ua.code2!=ub.code2){return ua.code2>ub.code2 ? 1 : -1;}
		}
		final byte[] a=ua.r.bases, b=ub.r.bases;
		if(a==b){return 0;}
		if(a==null || b==null){return a==null ? -1 : 1;}
		
		if(ua.canonical()==ub.canonical()){
			for(int i=0; i<a.length; i++){
				final byte aa=a[i], bb=b[i];
				if(aa!=bb && aa!='N' && bb!='N'){return aa-bb;}
			}
		}else{
			for(int i=0, j=b.length-1; i<a.length; i++, j--){
				final byte aa=a[i], bb=AminoAcid.baseToComplementExtended[b[j]];
				if(aa!=bb && aa!='N' && bb!='N'){return aa-bb;}
			}
		}
		return 0;
	}
	
	private static long hashTip(byte[] bases, boolean prefix, int k){
		if(bases==null || bases.length<k){return -1;}

		final long shift=2*k;
		final long shift2=shift-2;
		final long mask=~((-1L)<<shift);
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		final int start=(prefix ? 0 : bases.length-k);
		final int stop=start+k;

//		if(verbose){
//			System.err.println("\n"+new String(bases));
//			System.err.println("prefix="+prefix+", start="+start+", stop="+stop);
////			System.err.print(new String(bases));
//		}
		for(int i=start; i<stop; i++){
			byte b=bases[i];
//			if(verbose){System.err.print((char)b);}
			long x=baseToNumber[b];
			long x2=baseToRcompNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			len++;
		}
		if(verbose){System.err.println(new String(bases, start, k)+" = "+Tools.max(kmer, rkmer));}
		assert(len==k) : len+","+k;
		return Tools.max(kmer, rkmer);
	}
	
	private static final int calcMaxEdits(int maxEdits, float minIdentityMult, int len){
		return minIdentityMult==0 ? maxEdits : Tools.max(maxEdits, (int)Math.round(len*minIdentityMult));
	}
	
	
	private class Unit implements Comparable<Unit>{
		
		public Unit(Read r_){
			this(r_, isCanonical(r_.bases));
		}

		public Unit(Read r_, boolean canonical_){
//			this(r_, canonical_, canonical_ ? hash(r_.bases) : hashReversed(r_.bases));
			this(r_, canonical_, hash(r_.bases), hashReversed(r_.bases));
		}
		
		public Unit(Read r_, boolean canonical_, long codeF_, long codeR_){
			r=r_;
			code1=Tools.min(codeF_, codeR_);
			code2=Tools.max(codeF_, codeR_);
			long f=r.bases.length;
			prefix=hashTip(r.bases, true, k);
			suffix=hashTip(r.bases, false, k);
			if(canonical_){f|=CANON_MASK;}
			flags=f;
			assert(canonical()==canonical_);
			assert(length()==r.bases.length);
		}
		
		/**
		 * @param u2
		 * @param loc
		 * @param key
		 * @return
		 */
		public boolean contains(Unit u2, int loc, LongM key) {
			if(verbose){System.err.println("Considering key "+key+", unit "+u2);}
			if(minLengthPercent>0 && (u2.length()*100f/length())<minLengthPercent){return false;}
			assert(u2.code1!=code1 || u2.code2!=code2 || u2.length()!=length() || 
					(canonical()==u2.canonical() ? (u2.prefix!=prefix && u2.suffix!=suffix) : (u2.prefix!=suffix && u2.suffix!=prefix))) : 
						"Collision? \n"+this+"\n"+u2+"\n"+r+"\n"+u2.r;
			if(key.value()==u2.prefix){
				if(verbose){System.err.println("Containment A");}
				return containsForward(u2, loc-k2) || containsReverseRC(u2, loc);
			}
			if(verbose){System.err.println("Containment B");}
			assert(key.value()==u2.suffix);
			return containsReverse(u2, loc) || containsForwardRC(u2, loc-k2);
		}
		
		private boolean containsForward(Unit u2, int start) {
			if(start+u2.length()>length()){return false;}
//			if(true){return false;}
			if(u2.r!=null){
				final byte[] a=bases(), b=u2.bases();
				int mismatches=0, maxMismatches=calcMaxEdits(maxEdits, minIdentityMult, b.length);
				
				for(int i=start, j=0; j<b.length; i++, j++){
					byte aa=a[i];
					byte bb=b[j];
					if(aa!=bb){
						if(exact || (AminoAcid.isFullyDefined(aa) && AminoAcid.isFullyDefined(bb))){
							if((mismatches=mismatches+1)>maxMismatches){return false;}
						}
					}
				}
				return true;
			}else{
				assert(false) : "TODO: Verify by hashing and checking both tips";
				return false;
			}
		}
		
		private boolean containsForwardRC(Unit u2, int start) {
			if(ignoreReverseComplement){return false;}
			if(start+u2.length()>length()){return false;}
//			if(true){return false;}
			if(u2.r!=null){
				final byte[] a=bases(), b=u2.bases();
				int mismatches=0, maxMismatches=calcMaxEdits(maxEdits, minIdentityMult, b.length);

				for(int i=start, j=b.length-1; j>=0; i++, j--){
					byte aa=a[i];
					byte bb=AminoAcid.baseToComplementExtended[b[j]];
					if(aa!=bb){
						if(exact || (AminoAcid.isFullyDefined(aa) && AminoAcid.isFullyDefined(bb))){
							if((mismatches=mismatches+1)>maxMismatches){return false;}
						}
					}
				}
				return true;
			}else{
				assert(false) : "TODO: Verify by hashing and checking both tips";
				return false;
			}
		}
		
		private boolean containsReverse(Unit u2, int start) {
			if(start+1<u2.length()){return false;}
//			if(true){return false;}
			if(u2.r!=null){
				final byte[] a=bases(), b=u2.bases();
				int mismatches=0, maxMismatches=calcMaxEdits(maxEdits, minIdentityMult, b.length);

				for(int i=start, j=b.length-1; j>=0; i--, j--){
					byte aa=a[i];
					byte bb=b[j];
					if(aa!=bb){
						if(exact || (AminoAcid.isFullyDefined(aa) && AminoAcid.isFullyDefined(bb))){
							if((mismatches=mismatches+1)>maxMismatches){return false;}
						}
					}
				}
				return true;
			}else{
				assert(false) : "TODO: Verify by hashing and checking both tips";
				return false;
			}
		}
		
		private boolean containsReverseRC(Unit u2, int start) {
			if(ignoreReverseComplement){return false;}
			if(start+1<u2.length()){return false;}
//			if(true){return false;}
			if(u2.r!=null){
				final byte[] a=bases(), b=u2.bases();
				int mismatches=0, maxMismatches=calcMaxEdits(maxEdits, minIdentityMult, b.length);

				for(int i=start, j=0; j<b.length; i--, j++){
					byte aa=a[i];
					byte bb=AminoAcid.baseToComplementExtended[b[j]];
					if(aa!=bb){
						if(exact || (AminoAcid.isFullyDefined(aa) && AminoAcid.isFullyDefined(bb))){
							if((mismatches=mismatches+1)>maxMismatches){return false;}
						}
					}
				}
				return true;
			}else{
				assert(false) : "TODO: Verify by hashing and checking both tips";
				return false;
			}
		}
		
		/**
		 * @param u2
		 * @param loc
		 * @param key
		 * @return
		 */
		public boolean overlaps(Unit u2, int loc, LongM key) {
			assert(false) : "TODO";
			if(verbose){System.err.println("Considering key "+key+", unit "+u2);}
			if(minLengthPercent>0){
				final int len1=length(), len2=u2.length();
				if(Tools.min(len1, len2)*100f/Tools.max(len1, len2)<minLengthPercent){return false;}
			}
			assert(u2.code1!=code1 || u2.code2!=code2 || u2.length()!=length() || 
					(canonical()==u2.canonical() ? (u2.prefix!=prefix && u2.suffix!=suffix) : (u2.prefix!=suffix && u2.suffix!=prefix))) : 
						"Collision? \n"+this+"\n"+u2+"\n"+r+"\n"+u2.r;
			if(key.value()==u2.prefix){
				if(verbose){System.err.println("Overlap A");}
				if(overlapsForward(u2, loc-k2)){
					if(verbose){System.err.println("Overlap AF");}
					return true;
				}
				if(overlapsReverseRC(u2, loc)){
					if(verbose){System.err.println("Overlap AR");}
					return true;
				}
			}else{
				if(verbose){System.err.println("Overlap B");}
				assert(key.value()==u2.suffix);
				if(overlapsForwardRC(u2, loc-k2)){
					if(verbose){System.err.println("Containment BF");}
					return true;
				}
				if(overlapsReverse(u2, loc)){
					if(verbose){System.err.println("Containment BR");}
					return true;
				}
			}
			return false;
		}
		
		private boolean overlapsForward(Unit u2, int start) {
			
			final int len1=length(), len2=u2.length();
			final int overlapLength=Tools.min(len1-start, len2);
			
			if(overlapLength<minOverlap){return false;}
			if(minOverlapPercent>0f && (overlapLength*100f/Tools.min(len1, len2))<minLengthPercent){return false;}
			
			final byte[] a=bases(), b=u2.bases();
			assert(a!=null && b!=null) : "Null bases for "+code1+" or "+u2.code1;
			int mismatches=0, maxMismatches=calcMaxEdits(maxEdits, minIdentityMult, overlapLength);
			
			for(int i=start, j=0; j<overlapLength; i++, j++){
				byte aa=a[i];
				byte bb=b[j];
				if(aa!=bb){
					if(exact || (AminoAcid.isFullyDefined(aa) && AminoAcid.isFullyDefined(bb))){
						if((mismatches=mismatches+1)>maxMismatches){return false;}
					}
				}
			}
			return true;
		}
		
		private boolean overlapsForwardRC(Unit u2, int start) {
			if(ignoreReverseComplement){return false;}
			final int len1=length(), len2=u2.length();
			final int overlapLength=Tools.min(len1-start, len2);
			
			if(overlapLength<minOverlap){return false;}
			if(minOverlapPercent>0f && (overlapLength*100f/Tools.min(len1, len2))<minLengthPercent){return false;}
			
			final byte[] a=bases(), b=u2.bases();
			assert(a!=null && b!=null) : "Null bases for "+code1+" or "+u2.code1;
			int mismatches=0, maxMismatches=calcMaxEdits(maxEdits, minIdentityMult, b.length);
			
			for(int i=start, j=b.length-1, lim=start+overlapLength; i<lim; i++, j--){
				byte aa=a[i];
				byte bb=AminoAcid.baseToComplementExtended[b[j]];
				if(aa!=bb){
					if(exact || (AminoAcid.isFullyDefined(aa) && AminoAcid.isFullyDefined(bb))){
						if((mismatches=mismatches+1)>maxMismatches){return false;}
					}
				}
			}
			return true;
		}
		
		private boolean overlapsReverse(Unit u2, int start) {
			
			final int len1=length(), len2=u2.length();
			final int overlapLength=Tools.min(len1-start, len2);
			
			if(overlapLength<minOverlap){return false;}
			if(minOverlapPercent>0f && (overlapLength*100f/Tools.min(len1, len2))<minLengthPercent){return false;}
			
			final byte[] a=bases(), b=u2.bases();
			assert(a!=null && b!=null) : "Null bases for "+code1+" or "+u2.code1;
			int mismatches=0, maxMismatches=calcMaxEdits(maxEdits, minIdentityMult, b.length);

			for(int i=start, j=b.length-1; i>=0 && j>=0; i--, j--){
				byte aa=a[i];
				byte bb=b[j];
				if(aa!=bb){
					if(exact || (AminoAcid.isFullyDefined(aa) && AminoAcid.isFullyDefined(bb))){
						if((mismatches=mismatches+1)>maxMismatches){return false;}
					}
				}
			}
			return true;
		}
		
		private boolean overlapsReverseRC(Unit u2, int start) {
			if(ignoreReverseComplement){return false;}
			final int len1=length(), len2=u2.length();
			final int overlapLength=Tools.min(len1-start, len2);
			
			if(overlapLength<minOverlap){return false;}
			if(minOverlapPercent>0f && (overlapLength*100f/Tools.min(len1, len2))<minLengthPercent){return false;}

			final byte[] a=bases(), b=u2.bases();
			assert(a!=null && b!=null) : "Null bases for "+code1+" or "+u2.code1;
			int mismatches=0, maxMismatches=calcMaxEdits(maxEdits, minIdentityMult, b.length);

			for(int i=start, j=0; j<overlapLength; i--, j++){
				byte aa=a[i];
				byte bb=AminoAcid.baseToComplementExtended[b[j]];
				if(aa!=bb){
					if(exact || (AminoAcid.isFullyDefined(aa) && AminoAcid.isFullyDefined(bb))){
						if((mismatches=mismatches+1)>maxMismatches){return false;}
					}
				}
			}
			return true;
		}
		
		@Override
		public int compareTo(Unit b) {
			return compareRC(this, b);
		}
		
		public boolean equals(Object b){return equals((Unit)b);}
		public boolean equals(Unit b){return equalsRC(this, b);}
		
		@Override
		public int hashCode(){
			return (int)((code1^(code1>>>32))&0xFFFFFFFFL);
		}
		
		private synchronized void setInvalid(){
			assert(valid());
			flags|=VALID_MASK;
		}
		
		public byte[] bases(){return r==null ? null : r.bases;}

		public String name(){return r!=null ? r.id : null /*code+""*/;}
		public String toString(){return "("+code1+","+code2+","+length()+","+prefix+","+suffix+","+(canonical()?"c":"nc")/*+","+(tipCanonical()?"tc":"ntc")*/+")";}
		
		
		public final Read r;
		public final long code1;
		public final long code2;
		public final long prefix;
		public final long suffix;
//		private boolean valid=true;
		
		private long flags;
		/** True if the original read orientation was canonical */
		public final boolean canonical(){return (CANON_MASK&flags)!=0;}
		/** True if the original read orientation was canonical */
		public final boolean valid(){return (VALID_MASK&flags)==0;}
//		/** True if the original read tip orientation was canonical, false if the prefix and suffix were swapped */
//		public final boolean tipCanonical(){return (CANON_TIP_MASK&flags)!=0;}
		/** True if the original read tip orientation was canonical, false if the prefix and suffix were swapped */
		public final int length(){return (int)(LEN_MASK&flags);}

		private static final long LEN_MASK=0x7FFFFFFFL;	
		private static final long CANON_MASK=(1L<<33);		
		private static final long VALID_MASK=(1L<<34);
	}
	
	private ConcurrentReadStreamInterface crisa[];
	
	private String[] in=null;
	private String out=null;
	private int maxNs=-1;
	private long maxReads=-1;
	public boolean errorState=false;
	boolean sort=false;
	boolean ascending=true;
	boolean testContainment=true;
	boolean testMatch=true;
	boolean storeSuffix=false;
	boolean storeName=true;
	boolean storeQuality=true;
	boolean exact=true;
	boolean uniqueNames=true;
	private boolean multipleInputFiles=false;
	
	int maxEdits=0;
	float minIdentity=100;
	float minIdentityMult=0;
	float minLengthPercent=0;
	int minOverlap=200;
	float minOverlapPercent=0;

	long readsProcessed=0;
	long basesProcessed=0;
	long collisions=0;
	long containments=0;
	long containmentCollisions=0;
	long matches=0;
	long basematches=0;
	long basecontainments=0;
	long addedToMain=0;
	
	int k=31;
	int k2=k-1;
	
	private static int tcount=0;
	
	private HashSet<Unit> codeMap=new HashSet<Unit>(4000000);
	private HashMap<LongM, ArrayList<Unit>> affixMap=null;

	private static final long[][] hashcodes=makeCodes2(32);
	public static final byte[] baseToNumber=new byte[128];
	public static final byte[] baseToRcompNumber=new byte[128];
	public static final byte[] baseToRcomp=new byte[128];
	private static PrintStream outstream=System.err;
	public static boolean overwrite=false;
	public static boolean showSpeed=true;
	public static boolean verbose=false;
	public static boolean ignoreReverseComplement=false;
	public static int MINSCAF=0;
	public static int THREADS=Shared.THREADS;
	public static int threadMaxReadsToBuffer=4000;
	public static int threadMaxBasesToBuffer=32000000;
	public static boolean PAUSE=false;
	
	static{//All others are 0
		baseToNumber['A']=baseToNumber['a']=0;
		baseToNumber['C']=baseToNumber['c']=1;
		baseToNumber['G']=baseToNumber['g']=2;
		baseToNumber['T']=baseToNumber['t']=3;
		
		baseToRcompNumber['A']=baseToRcompNumber['a']=3;
		baseToRcompNumber['C']=baseToRcompNumber['c']=2;
		baseToRcompNumber['G']=baseToRcompNumber['g']=1;
		baseToRcompNumber['T']=baseToRcompNumber['t']=0;
	}
	
}
