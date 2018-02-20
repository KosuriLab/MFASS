package jgi;

import java.util.ArrayList;
import java.util.Arrays;

import kmer.KCountArray;
import kmer.KmerCount4;
import kmer.KmerCount6MT;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.Read;

import dna.AminoAcid;
import dna.Timer;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;

import align2.ListNum;
import align2.Tools;

/**
 * @author Brian Bushnell
 * @date Jul 31, 2012
 *
 */
public class FilterReads {
	
	public static void main(String[] args){
		System.err.println("Executing "+(new Object() { }.getClass().getEnclosingClass().getName())+" "+Arrays.toString(args)+"\n");
		Timer t=new Timer();
		t.start();
		String fname1=args[0];
		String fname2=("null".equals(args[1]) ? null : args[1]);
		
		String outName=args[2];
		assert(!fname1.equalsIgnoreCase(outName));
		assert(fname2==null || !fname2.equalsIgnoreCase(outName));
//		FASTQ.TEST_INTERLEAVED=false;
//		FASTQ.PARSE_CUSTOM=false;
		
		long readsOut=-1;
		int sections=1;
		int hashes=1;
		int passes=1;
		int matrixbits=33;
		
		ArrayList<String> extraFiles=new ArrayList<String>();
		
		for(int i=3; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=(split.length>1 ? split[1] : "true");
			
			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("reads")){
				readsOut=Long.parseLong(b);
			}else if(a.equals("readlen")){
				MIN_LEN=Integer.parseInt(b);
				MIN_LEN_2=(MIN_LEN*2+2)/3;
			}else if(a.contains("testinterleaved")){
				FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
			}else if(a.contains("forceinterleaved")){
				FASTQ.FORCE_INTERLEAVED=Tools.parseBoolean(b);
			}else if(a.contains("testkmer") || a.contains("kmertest")){
				TEST_KMERS=Tools.parseBoolean(b);
			}else if(a.startsWith("rcomp") || a.startsWith("reverse")){
				RCOMP=Tools.parseBoolean(b);
			}else if(a.equals("k") || a.equals("kmer") || a.equals("kmerlen")){
				K=Integer.parseInt(b);
			}else if(a.equals("cbits")){
				CBITS=Integer.parseInt(b);
			}else if(a.equals("hashes")){
				hashes=Integer.parseInt(b);
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}else if(a.equals("matrixbits")){
				matrixbits=Integer.parseInt(b);
			}else if(a.equals("minq") || a.equals("minquality")){
				MIN_QUALITY_READ=Byte.parseByte(b);
				MIN_QUALITY_KMER=Byte.parseByte(b);
			}else if(a.startsWith("minqr") || a.startsWith("minqualityr")){
				MIN_QUALITY_READ=Byte.parseByte(b);
			}else if(a.startsWith("minqk") || a.startsWith("minqualityk")){
				MIN_QUALITY_KMER=Byte.parseByte(b);
			}else if(a.startsWith("minavgq")){
				MIN_AVG_QUALITY=Byte.parseByte(b);
			}else if(a.startsWith("maxerr")){
				MAX_EXPECTED_ERRORS=Float.parseFloat(b);
			}else if(a.equals("lowthresh") || a.equals("lowkmerthresh")){
				LOW_KMER_THRESH=Integer.parseInt(b);
			}else if(a.startsWith("extrafile")){
				String[] sp2=b.split(",");
				for(String s : sp2){
					extraFiles.add(s);
				}
			}else if(a.startsWith("parsecustom")){
				FASTQ.PARSE_CUSTOM=Tools.parseBoolean(b);
			}else if(a.equals("sections")){
				sections=Integer.parseInt(b);
				assert(sections>0);
			}else{
				throw new RuntimeException(args[i]+"\n"+a+"\n"+b);
			}
		}
		if(readsOut<1){readsOut=Long.MAX_VALUE;}
		if(sections>1){assert(outName.contains("SECTION#"));}
		
		KCountArray kca=null;
		if(TEST_KMERS){
			Timer ht=new Timer();
			ht.start();
			
			long maxreads=(readsOut==Long.MAX_VALUE ? -1 : Tools.max(readsOut, readsOut*8*sections));
			KmerCount6MT.maxReads=maxreads;
			int kbits=2*K;
			matrixbits=Tools.min(kbits, matrixbits);
			
			kca=KmerCount6MT.makeKca(fname1, fname2, extraFiles, K, CBITS, 0, matrixbits, hashes, MIN_QUALITY_KMER, RCOMP, maxreads, passes, 1, 1, 2);

			KmerCount6MT.printStatistics(kca);
			MASK=~((-1L)<<(kbits));
			ht.stop();
			System.out.println("Finished hashing "+KmerCount4.keysCounted+" kmers.");
			System.out.println("Hashing time: "+ht);
			
		}
		
		final ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(fname1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(fname2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(-1, false, true, ff1, ff2);
			Thread th=new Thread(cris);
			th.start();
		}
		
//		assert(false) : cris.getClass();
		long kept=0;
		
		boolean sf=outName.contains("SECTION#");
		TextStreamWriter tswAll1=sf ? new TextStreamWriter(outName.replace("SECTION#", "ALL").replaceFirst("#", "1"), true, false, true) : null;
		TextStreamWriter tswAll2=sf && cris.paired() ? new TextStreamWriter(outName.replace("SECTION#", "ALL").replaceFirst("#", "2"), true, false, true) : null;
		TextStreamWriter tswBad1=sf ? new TextStreamWriter(outName.replace("SECTION#", "BAD").replaceFirst("#", "1"), true, false, true) : null;
		TextStreamWriter tswBad2=sf && cris.paired() ? new TextStreamWriter(outName.replace("SECTION#", "BAD").replaceFirst("#", "2"), true, false, true) : null;
		
		if(tswAll1!=null){
			(tswAll1).start();
			tswAll1.print("#"+Read.header()+"\n");
		}
		if(tswAll2!=null){
			(tswAll2).start();
			tswAll2.print("#"+Read.header()+"\n");
		}
		if(tswBad1!=null){
			(tswBad1).start();
			tswBad1.print("#"+Read.header()+"\n");
		}
		if(tswBad2!=null){
			(tswBad2).start();
			tswBad2.print("#"+Read.header()+"\n");
		}
		
		for(int i=1; i<=sections; i++){
			long x=process(cris, outName.replace("SECTION#", i+""), readsOut, kca, tswAll1, tswAll2, tswBad1, tswBad2);
			kept+=x;
			if(x<1){break;}
		}
		ReadWrite.closeStream(cris);

		if(tswAll1!=null){tswAll1.poison();}
		if(tswAll2!=null){tswAll2.poison();}
		if(tswBad1!=null){tswBad1.poison();}
		if(tswBad2!=null){tswBad2.poison();}
		
		System.out.println("Removed reads:                  \t"+removed+String.format("    \t%.2f%%", removed*100d/(removed+kept)));
		System.out.println("Removed short reads:            \t"+removedShort);
		System.out.println("Removed low avg quality reads:  \t"+removedAvgQuality);
		System.out.println("Removed low min quality reads:  \t"+removedMinQuality);
		System.out.println("Removed N-containing reads:     \t"+removedNocall);
		System.out.println("Removed error-prone reads:      \t"+removedErrors);
		System.out.println("Removed low-kmer reads:         \t"+removedKmer);
		
		t.stop();
		System.out.println("Time:\t"+t);
	}
	
	
	public static long process(ConcurrentReadStreamInterface cris, String outfile, final long readsOut, final KCountArray kca2, 
			final TextStreamWriter tswAll1, final TextStreamWriter tswAll2, final TextStreamWriter tswBad1, final TextStreamWriter tswBad2){

		
		ListNum<Read> ln=cris.nextList();
		if(ln==null || ln.list==null || ln.list.isEmpty()){return 0;}
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		TextStreamWriter tsw1=new TextStreamWriter(outfile.replaceFirst("#", "1"), true, false, true);
		tsw1.start();
		tsw1.print("#"+Read.header()+"\n");
		TextStreamWriter tsw2=cris.paired() ? new TextStreamWriter(outfile.replaceFirst("#", "2"), true, false, true) : null;
		if(tsw2!=null){
			assert(!tsw1.fname.equalsIgnoreCase(tsw2.fname));
			tsw2=(tsw2);
			tsw2.start();
			tsw2.print("#"+Read.header()+"\n");
		}
		
		
		long readsKept=0;
		
		while(reads!=null && reads.size()>0 && readsKept<readsOut){
			for(Read r : reads){
				
				Read r1x=null;  //Copy of untrimmed read
//				r1x=new Read(r.bases, r.chrom, r.start, r.stop, r.id, r.quality, r.numericID, r.flags);
//				if(r.mate!=null){r1x.mate=new Read(r.mate.bases, r.mate.chrom, r.mate.start, r.mate.stop, r.mate.id, r.mate.quality, r.mate.numericID, r.mate.flags);}
				
				trim(r);
				if(r.mate!=null){trim(r.mate);}
				if(passesFilter(r, kca2)){
					StringBuilder sb=r.toText(false).append('\n');
					tsw1.print(sb);
					if(tswAll1!=null){tswAll1.print(sb);}
					if(r.mate!=null){
						sb=r.mate.toText(false).append('\n');
						if(tsw2==null){tsw1.print(sb);}
						else{tsw2.print(sb);}
						if(tswAll2==null){if(tswAll1!=null){tswAll1.print(sb);}}
						else{tswAll2.print(sb);}
					}
					readsKept++;
//					System.out.print(", "+readsOut);
				}else if(r1x!=null){
					StringBuilder sb=r1x.toText(false).append('\n');
					tswBad1.print(sb);
					if(r1x.mate!=null){
						sb=r1x.mate.toText(false).append('\n');
						if(tswBad2==null){tswBad1.print(sb);}
						else{tswBad2.print(sb);}
					}
				}else{
					
					if(r.bases!=null && r.bases.length>=MIN_LEN_2 && (r.mate==null || (r.mate.bases!=null && r.bases.length>=MIN_LEN_2))){
						r.setDiscarded(false);
						assert(r!=null);
						assert(r.toText(false)!=null);
						StringBuilder sb=r.toText(false).append('\n');
						if(tswBad1!=null){tswBad1.print(sb);}
						if(r.mate!=null){
							r.mate.setDiscarded(false);
							sb=r.mate.toText(false).append('\n');
							if(tswBad2==null){if(tswBad1!=null){tswBad1.print(sb);}}
							else{tswBad2.print(sb);}
						}
					}
				}
				if(readsKept>=readsOut){break;}
			}
			cris.returnList(ln, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
//		System.out.println("Returning last list.");
		cris.returnList(ln, ln.list.isEmpty());
//		System.out.println("Closing.");
//		cris.close();
//		System.out.println("Closed.");
		
//		ln=cris.nextList();
//		for(int i=0; i<100; i++){
//			ln=cris.nextList();
//			ln.list.clear();
//			System.out.println("Returned extra list.");
//			cris.returnList(ln, true);
//		}
		
		tsw1.poison();
		if(tsw2!=null){tsw2.poison();}
		
		if(tsw1.isAlive()){
			try {
				synchronized(tsw1){
					tsw1.join();
				}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(tsw2!=null && tsw2.isAlive()){
			try {
				synchronized(tsw2){
					tsw2.join();
				}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		return readsKept;
	}

	
	public static void trim(Read r) {
		assert(r.bases.length>=MIN_LEN) : r.bases.length;
		byte[] bases=r.bases;
		byte[] quals=r.quality;
		if(bases.length<MIN_LEN_2){
//			r.bases=null;
//			r.quality=null;
//			r.match=null;
			assert(false);
			r.setDiscarded(true);
			return;
		}
		
		int left=0;
		int right=bases.length-1;
		int safe=0;
		final int minsafe=8;
		
		assert(r.bases.length>=MIN_LEN) : r.bases.length;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			byte q=quals[i];
			if(b=='N' || q<MIN_QUALITY_READ){
				safe=0;
				left=i+1;
			}else{
				safe++;
				if(safe>=minsafe){break;}
			}
		}
		
		safe=0;
		
		for(int i=bases.length-1; i>0; i--){
			byte b=bases[i];
			byte q=quals[i];
			if(b=='N' || q<MIN_QUALITY_READ){
				safe=0;
				right=i-1;
			}else{
				safe++;
				if(safe>=minsafe){break;}
			}
		}
		
		if(right-left+1<MIN_LEN_2)
		
		if(left==0 && right==bases.length-1){
			//do nothing
			return;
		}else if(right-left+1<MIN_LEN_2){
			r.match=null;
			r.bases=null;
			r.quality=null;
			r.setMapped(false);
			r.setDiscarded(true);
//			r.clearSite();
			return;
		}else if(right-left+1<MIN_LEN){
			r.setDiscarded(true);
		}
		r.bases=Arrays.copyOfRange(bases, left, right+1);
		r.quality=Arrays.copyOfRange(quals, left, right+1);
		r.mapLength=r.bases.length;
		assert(r.bases.length>=MIN_LEN || r.discarded());
	}


	public static boolean passesFilter(Read r, KCountArray kca2){
		boolean b=passesFilter2(r, kca2) && (r.mate==null || passesFilter2(r.mate, kca2));
		if(!b){removed++;}
		return b;
	}
	
	public static boolean passesFilter2(Read r, KCountArray kca2){
		
		if(r.bases==null || r.bases.length<MIN_LEN){removedShort++; return false;}
		if(r.avgQuality()<MIN_AVG_QUALITY){removedAvgQuality++; return false;}
		final int minq=r.minQualityLastNBases(r.bases.length);
		if(minq<MIN_QUALITY_READ){removedMinQuality++; return false;}
		final int n=r.countNocalls();
		if(n>MAX_N){removedNocall++; return false;}
		final float errors=r.expectedErrors();
		if(errors>MAX_EXPECTED_ERRORS){removedErrors++; return false;}
		
		if(r.discarded()){removedShort++; return false;} //Removed for some unknown reason; probably too short.
		
		if(TEST_KMERS && kca2!=null){
			int len=0;
			long kmer=0;
			int low=0, low2=0, high=0;
			final int lim2=Tools.min(LOW_KMER_THRESH+3, kca2.maxValue-1);
			byte[] bases=r.bases;
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				int x=AminoAcid.baseToNumber[b];
				if(x<0){//'N' encountered.
					len=0;
					kmer=0;
				}else{
					kmer=((kmer<<2)|x)&MASK;
					len++;
					if(len>=K){
						int count=kca2.read(kmer);
						if(count<=LOW_KMER_THRESH){low++;}
						if(count<=lim2){low2++;}
						else{high++;}
					}
				}
			}
			if(low>=MIN_LOW_KMERS_TO_TOSS){removedKmer++; return false;}
			if(high<1){removedKmer++; return false;} //Too many 'N' to make a useful read
			if(low2>0 && (r.expectedErrors()>0.20f || n>0 || minq<=14)){
				return false;
			}
		}
		
		return true;
	}
	

	public static int MIN_LEN=60;
	public static int MIN_LEN_2=45;
	public static int MAX_N=1;
	public static int MIN_AVG_QUALITY=19;
	
	/** Quality for trimming/keeping reads */
	public static byte MIN_QUALITY_READ=6;
	
	/** Quality for generating trusted kmers */
	public static byte MIN_QUALITY_KMER=11;
	
	public static float MAX_EXPECTED_ERRORS=1.2f;
	public static boolean TEST_KMERS=true;
	public static int K=17;
	public static int CBITS=2;
	private static long MASK=0;
	private static boolean RCOMP=true;
	
	/** Kmers with this frequency or less are considered "low" */
	public static int LOW_KMER_THRESH=1;
	/** A read must have at least this many low kmers to toss it.  Generally, a read with an error should have multiple low kmers, unless they are at the tip. */
	public static int MIN_LOW_KMERS_TO_TOSS=2;
	
	//TODO: RARE stuff is not implemented 
	
	/** Kmers with at least RARE_KMER_MIN_THRESH and at most RARE_KMER_MAX_THRESH are considered "rare" */
	public static int RARE_KMER_MIN_THRESH=8;
	/** Kmers with at least RARE_KMER_MIN_THRESH and at most RARE_KMER_MAX_THRESH are considered "rare" */
	public static int RARE_KMER_MAX_THRESH=32;
	/** A read must have at least this many more rare kmers than low kmers to retain it despite being suspect. */
	public static int MIN_RARE_KMERS_TO_KEEP=1;
	public static boolean RETAIN_RARE_KMERS=false;

	private static long removedShort=0;
	private static long removedAvgQuality=0;
	private static long removedMinQuality=0;
	private static long removedNocall=0;
	private static long removedErrors=0;
	private static long removedKmer=0;
	private static long removed=0;
	
	
}
