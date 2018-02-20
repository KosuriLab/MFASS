package kmer;

import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import jgi.ErrorCorrect;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;

import align2.ListNum;
import align2.Shared;
import align2.Tools;
import dna.AminoAcid;
import dna.Timer;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 *
 */
public class KmerCount7MTA {
	
public static void main(String[] args){
		
		Timer t=new Timer();
		t.start();
		
		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		int k=14;
		int cbits=16;
		int gap=0;
		int matrixbits=-1;
		int hashes=1;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=(split.length>1 ? split[1] : "true");

			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.startsWith("gap")){
				gap=Integer.parseInt(b);
			}else if(a.startsWith("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.startsWith("matrixbits")){
				matrixbits=Integer.parseInt(b);
			}else if(a.startsWith("hashes")){
				hashes=Integer.parseInt(b);
			}else if(a.equals("canonical")){
				CANONICAL=Tools.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		int kbits=Tools.min(2*k, 62);
		if(matrixbits<0){
			matrixbits=kbits;
		}
		matrixbits=Tools.min(kbits, matrixbits);
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			FastaReadInputStream.TARGET_READ_LEN=300000000;
			FastaReadInputStream.MIN_READ_LEN=k;
		}else{
			FASTQ.PARSE_CUSTOM=false;
//			assert(false) : FASTQ.PARSE_CUSTOM;
		}
		
		KCountArray count=KCountArray.makeNew(1L<<kbits, 1L<<matrixbits, cbits, gap, hashes);
		try {
			count=countFastq(fname1, fname2, k, cbits, gap, true, count);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		count.shutdown();
		
//		verbose=true;
		
		t.stop();
		System.out.println("Finished counting; time = "+t);
		
		printStatistics(count);
		
	}

	public static void printStatistics(KCountArray count){
		long[] freq=count.transformToFrequency();

//		System.out.println(count+"\n");
//		System.out.println(Arrays.toString(freq)+"\n");
		
		long sum=sum(freq);
		System.out.println("Kmer fraction:");
		int lim1=8, lim2=16;
		for(int i=0; i<lim1; i++){
			String prefix=i+"";
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+String.format("%.3f%%   ",(100l*freq[i]/(double)sum))+"\t"+freq[i]);
		}
		while(lim1<=freq.length){
			int x=0;
			for(int i=lim1; i<lim2; i++){
				x+=freq[i];
			}
			String prefix=lim1+"-"+(lim2-1);
			if(lim2>=freq.length){prefix=lim1+"+";}
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+String.format("%.3f%%   ",(100l*x/(double)sum))+"\t"+x);
			lim1*=2;
			lim2=min(lim2*2, freq.length);
		}
		
		long sum2=sum-freq[0];
		long x=freq[1];
		System.out.println();
		System.out.println("Keys Counted:  \t         \t"+keysCounted);
		System.out.println("Unique:        \t         \t"+sum2);
		System.out.println("Avg Sites/Key: \t         \t"+String.format("%.3f    ",(keysCounted*1d/sum2)));
		System.out.println();
		System.out.println("Singleton:     \t"+String.format("%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
		x=sum2-x;
		System.out.println("Useful:        \t"+String.format("%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles, int k, int cbits){
		return makeKca(fname1, fname2, extraFiles, k, cbits, 0, Tools.min(2*k, 35), 1, minQuality, true, maxReads, 1, 1, 1, 2);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles, 
			int k, int cbits, int gap, int matrixbits, int hashes, int minqual, boolean rcomp, long maxreads){
		assert(matrixbits<63);
		return makeKca(fname1, fname2, extraFiles, k, cbits, gap, matrixbits, hashes, minqual, rcomp, maxreads, 1, 1, 1, 2);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles, 
			int k, int cbits, int gap, int matrixbits, int hashes, int minqual, boolean rcomp, long maxreads, int passes, int stepsize, int thresh1, int thresh2){
		assert(matrixbits<63);
		return makeKca(fname1, fname2, extraFiles, 
				k, cbits, gap, 1L<<matrixbits, hashes, minqual, rcomp, maxreads, passes, stepsize, thresh1, thresh2, null);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles, 
			int k, int cbits, int gap, long cells, int hashes, int minqual, boolean rcomp, long maxreads, int passes, int stepsize, int thresh1, int thresh2){
		return makeKca(fname1, fname2, extraFiles, 
				k, cbits, gap, cells, hashes, minqual, rcomp, maxreads, passes, stepsize, thresh1, thresh2, null);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles, 
			int k, int cbits, int gap, long cells, int hashes, int minqual, boolean rcomp, long maxreads, int passes, int stepsize, int thresh1, int thresh2,
			KCountArray prefilter){
		final int kbits=Tools.min(2*k, 62);
//		verbose=true;
		if(verbose){System.err.println("Making kca from ("+fname1+", "+fname2+")\nk="+k+", gap="+gap+", cells="+Tools.toKMG(cells)+", cbits="+cbits);}
		
		boolean oldsplit=FastaReadInputStream.SPLIT_READS;
		long oldmax=maxReads;
		byte oldq=minQuality;
		maxReads=maxreads;
		minQuality=(byte)minqual;
		//	System.out.println("kbits="+(kbits)+" -> "+(1L<<kbits)+", matrixbits="+(matrixbits)+" -> "+(1L<<matrixbits)+", cbits="+cbits+", gap="+gap+", hashes="+hashes);
		KCountArray kca=KCountArray.makeNew(1L<<kbits, cells, cbits, gap, hashes, prefilter);
		
//		System.out.println("a");
		{//For processing input lists
			ArrayList<String> extra2=null;
			if(fname1!=null && fname1.contains(",")){
				String[] s=fname1.split(",");
				if(extra2==null){extra2=new ArrayList<String>();}
				for(int i=1; i<s.length; i++){extra2.add(s[i]);}
				fname1=s[0];
			}
			if(fname2!=null && fname2.contains(",")){
				String[] s=fname2.split(",");
				if(extra2==null){extra2=new ArrayList<String>();}
				for(int i=1; i<s.length; i++){extra2.add(s[i]);}
				fname2=s[0];
			}
			if(extra2!=null){
				if(extraFiles!=null){
					for(String s : extraFiles){
						extra2.add(s);
					}
				}
				extraFiles=extra2;
			}
		}
//		System.out.println("b");
		
		if(extraFiles!=null){
			for(String s : extraFiles){
				if(fileIO.FileFormat.hasFastaExtension(s)){
					FastaReadInputStream.SPLIT_READS=false;
				}
			}
		}
		
//		System.out.println("c");
		if(passes==1){
//			System.out.println("c1");
			try {
				countFastq(fname1, fname2, k, cbits, gap, rcomp, kca);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					try {
						countFastq(s, null, k, cbits, gap, rcomp, kca);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			kca.shutdown();

		}else{
//			System.out.println("c2");
			assert(passes>1);
			KCountArray trusted=null;
			for(int i=1; i<passes; i++){
				boolean conservative=i>2;// /*or, alternately, (trusted==null || trusted.capacity()>0.3)
				int step=(stepsize==1 ? 1 : stepsize+i%2);
				//			if(!conservative){step=(step+3)/4;}
				if(!conservative){step=Tools.min(3, (step+3)/4);}

				try {
					countFastq(fname1, fname2, k, cbits, true, kca, trusted, maxreads, thresh1, step, conservative);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(extraFiles!=null){
					maxReads=-1;
					for(String s : extraFiles){
						try {
							countFastq(s, null, k, cbits, true, kca, trusted, maxreads, thresh1, step, conservative);
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
				kca.shutdown();
				
				System.out.println("Trusted:   \t"+kca.toShortString());
				trusted=kca;
				kca=KCountArray.makeNew(1L<<kbits, cells, cbits, gap, hashes, prefilter);

			}

			try {
				countFastq(fname1, fname2, k, cbits, true, kca, trusted, maxreads, thresh2, stepsize, true);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					try {
						countFastq(s, null, k, cbits, true, kca, trusted, maxreads, thresh2, stepsize, true);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			kca.shutdown();
		}
//		System.out.println("d");
		minQuality=oldq;
		maxReads=oldmax;
		FastaReadInputStream.SPLIT_READS=oldsplit;
		
		
		return kca;
	}
	
	public static KCountArray countFasta(String fname, int k, int cbits, int gap){
		assert(k>=1 && k<20);
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		final long cells=mask+1;
		if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
		final KCountArray count=KCountArray.makeNew(cells, cbits, gap);
		
		TextFile tf=new TextFile(fname, false, false);
		
		long kmer=0; //current kmer
		int len=0;  //distance since last contig start or ambiguous base
		
		String s=null;
		for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.charAt(0)!='>'){
				for(int i=0; i<s.length(); i++){
					char c=s.charAt(i);
					int x=AminoAcid.baseToNumber[c];
					if(x<0){
						len=0;
						kmer=0;
					}else{
						kmer=((kmer<<2)|x)&mask;
						len++;
						if(len>=k && (!CANONICAL || KCountArray.isCanonical(kmer, k))){
							keysCounted++;
							count.increment(kmer);
						}
					}
				}
			}
		}
		return count;
	}
	
	public static KCountArray countFastq(String reads1, String reads2, int k, int cbits, int gap, boolean rcomp, KCountArray count) throws Exception{
		assert(k>=1 && (count!=null || k<20));
		final int kbits=Tools.min(2*k, 62);
		final long mask=~((-1L)<<(kbits));
//		System.err.println("countFastq...  making a new cris");
		if(count==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			count=KCountArray.makeNew(cells, cbits, gap);
		}
		assert(gap==count.gap);
		
		final ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, ff2);
			Thread th=new Thread(cris);
			th.start();
		}
		
		assert(cris!=null) : reads1;
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
//		countFastq(cris, k, rcomp, count);
//		assert(false) : THREADS;
		CountFastqThread[] cta=new CountFastqThread[THREADS];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountFastqThread(cris, k, rcomp, count);
			cta[i].start();
		}
//		System.out.println("~1");
		for(int i=0; i<cta.length; i++){
//			System.out.println("~2");
			CountFastqThread ct=cta[i];
			synchronized(ct){
//				System.out.println("~3");
				while(ct.getState()!=State.TERMINATED){
//					System.out.println("~4");
					try {
						ct.join(2000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
//					System.out.println("~5");
				}
			}
		}
//		System.out.println("~6");
		
		ReadWrite.closeStream(cris);
		if(verbose){System.err.println("Closed stream");}
		if(verbose){System.err.println("Processed "+readsProcessed+" reads.");}

		
		return count;
	}
	

	
	

	
	public static KCountArray countFastq(final String reads1, final String reads2, final int k, final int cbits, final boolean rcomp, 
			KCountArray count, final KCountArray trusted, final long maxReads, final int thresh, final int detectStepsize, final boolean conservative)
			throws Exception{
		
		assert(k>=1 && (count!=null || k<20));
		final int kbits=Tools.min(2*k, 62);
		final long mask=~((-1L)<<(kbits));
		
//		System.out.println("k="+k+", kbits="+kbits+", mask="+Long.toHexString(mask)+", thresh="+thresh);
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);
		
//		verbose=true;
		
		if(count==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			count=KCountArray.makeNew(cells, cbits, 0);
		}
		
		final ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, ff2);
			Thread th=new Thread(cris);
			th.start();
		}
		
		assert(cris!=null) : reads1;
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		

//		countFastq(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);

//		assert(false) : THREADS;
		CountFastqThread[] cta=new CountFastqThread[THREADS];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountFastqThread(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);
			cta[i].start();
		}
		
		for(int i=0; i<cta.length; i++){
			CountFastqThread ct=cta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		cris.close();
		if(verbose){System.err.println("Closed stream");}
		
//		System.out.println("*** after ***");
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);
		
		return count;
	}
	
	public static long[] transformToFrequency(int[] count){
		long[] freq=new long[2000];
		int max=freq.length-1;
		for(int i=0; i<count.length; i++){
			int x=count[i];
			x=min(x, max);
			freq[x]++;
		}
		return freq;
	}
	
	public static long sum(int[] array){
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}
	
	public static long sum(long[] array){
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}
	
	private static class CountFastqThread extends Thread{
		
		CountFastqThread(final ConcurrentReadStreamInterface cris_, final int k_, final boolean rcomp_, final KCountArray count_){
			this(cris_, k_, rcomp_, count_, null, 2, 1, true);
		}
		
		CountFastqThread(final ConcurrentReadStreamInterface cris_, final int k_, final boolean rcomp_, 
		final KCountArray count_, final KCountArray trusted_, final int thresh_, final int detectStepsize_, final boolean conservative_){
			cris=cris_;
			k=k_;
			rcomp=rcomp_;
			count=count_;
			trusted=trusted_;
			thresh=thresh_;
			detectStepsize=detectStepsize_;
			conservative=conservative_;
		}
		
		public void run(){
//			System.out.println("Running");
			if(trusted==null){
				countFastq(cris, k, rcomp, count);
			}else{
				countFastq(cris, k, rcomp, count, trusted, thresh, detectStepsize,  conservative);
			}
			
			synchronized(getClass()){
				keysCounted+=keysCountedLocal;
				readsProcessed+=readsProcessedLocal;

				if(verbose){System.err.println(keysCounted+", "+keysCountedLocal);}
				if(verbose){System.err.println(readsProcessed+", "+readsProcessedLocal);}
			}
		}
		
		
		
		private final void countFastq(ConcurrentReadStreamInterface cris, int k, boolean rcomp, KCountArray count){
			assert(k>=1 && count!=null);
			
//			System.out.println("Waiting for list");
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
//			System.out.println("Got list: "+(ln==null ? "null" : ln.id)+", "+(ln==null || ln.list==null ? "null" : ln.list.size()));

			long[] array=null;
			if(count.gap==0){
				final int kbits=Tools.min(2*k, 62);
				final long mask=~((-1L)<<(kbits));
				
				
				while(reads!=null && reads.size()>0){
					//System.err.println("reads.size()="+reads.size());
					for(Read r : reads){
						
						if(readsamplerate<2 || r.numericID%readsamplerate==0){
							readsProcessedLocal++;
//							addRead(r, count, k, mask, rcomp);
//							if(r.mate!=null){
//								addRead(r.mate, count, k, mask, rcomp);
//							}
							array=addRead_Advanced(r, count, k, mask, array);
						}
//						System.out.println(r);
//						System.out.println("kmers hashed: "+keysCountedLocal);
					}
					//System.err.println("returning list");
					cris.returnList(ln, ln.list.isEmpty());
					//System.err.println("fetching list");
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
			}else{
				final int k1=(k+1)/2;
				final int k2=k/2;
				final int kbits1=2*k1;
				final int kbits2=2*k2;
				final int gap=count.gap;
				final long mask1=~((-1L)<<(kbits1));
				final long mask2=~((-1L)<<(kbits2));
				while(reads!=null && reads.size()>0){
					//System.err.println("reads.size()="+reads.size());
					for(Read r : reads){
						
						if(readsamplerate<2 || r.numericID%readsamplerate==0){
							readsProcessedLocal++;
							addReadSplit(r, count, k1, k2, mask1, mask2, gap, rcomp);
							if(r.mate!=null){
								addReadSplit(r.mate, count, k1, k2, mask1, mask2, gap, rcomp);
							}
						}

					}
					//System.err.println("returning list");
					cris.returnList(ln, ln.list.isEmpty());
					//System.err.println("fetching list");
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
			}
			
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln, ln==null ? true : ln.list.isEmpty());
			if(verbose){System.err.println("Returned list");}
		}
		
		


		private final void countFastq(final ConcurrentReadStreamInterface cris, final int k, final boolean rcomp, 
				final KCountArray count, final KCountArray trusted, final int thresh, final int detectStepsize, final boolean conservative){
			if(count.gap>0){countFastqSplit(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);}
			assert(k>=1 && (count!=null || k<20));
			final int kbits=Tools.min(2*k, 62);
			final long mask=~((-1L)<<(kbits));
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			long[] array=null;
			while(reads!=null && reads.size()>0){
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					
					Read r2=r.mate;
					{
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r, trusted, k, thresh, detectStepsize) : 
								ErrorCorrect.detectTrusted(r, trusted, k, thresh, detectStepsize));
//							System.out.println("\n"+toString(bs, r.bases.length));
//							System.out.println(new String(r.bases));
							if(bs!=null){
								for(int i=bs.nextClearBit(0); i<r.bases.length; i=bs.nextClearBit(i+1)){
									r.bases[i]='N';
									if(r.quality!=null){r.quality[i]=0;}
								}
							}
//							System.out.println(new String(r.bases));
//							System.out.println("used = "+String.format("%.3f%%",count.usedFraction()*100));
//							System.out.println("used = "+((KCountArray4)count).cellsUsed());
//							if(bs.length()<r.bases.length){r=null;}
						}
//						if(r!=null){addRead(r, count, k, mask, rcomp);}
					}
					if(r2!=null){
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r2, trusted, k, thresh, detectStepsize) : 
								ErrorCorrect.detectTrusted(r2, trusted, k, thresh, detectStepsize));
							if(bs!=null){
								for(int i=bs.nextClearBit(0); i<r2.bases.length; i=bs.nextClearBit(i+1)){
									r2.bases[i]='N';
									if(r2.quality!=null){r2.quality[i]=0;}
								}
							}
						}
					}
					array=addRead_Advanced(r, count, k, mask, array);

				}
				//System.err.println("returning list");
				cris.returnList(ln, ln.list.isEmpty());
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln, ln.list.isEmpty());
			if(verbose){System.err.println("Returned list");}
		}
		
		
		private void countFastqSplit(final ConcurrentReadStreamInterface cris, final int k, final boolean rcomp, 
				final KCountArray count, final KCountArray trusted, final int thresh, final int detectStepsize, final boolean conservative){
			assert(false) : cris.paired();
			assert(count.gap>0);
			assert(k<32 && k>=1 && (count!=null || k<20));
			final int kbits=2*k;
			final long mask=~((-1L)<<(kbits));
			

			final int k1=(k+1)/2;
			final int k2=k/2;
			final int kbits1=2*k1;
			final int kbits2=2*k2;
			final int gap=count.gap;
			final long mask1=~((-1L)<<(kbits1));
			final long mask2=~((-1L)<<(kbits2));
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(reads!=null && reads.size()>0){
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					
					Read r2=r.mate;
					{
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r, trusted, k, thresh, detectStepsize) : 
								ErrorCorrect.detectTrusted(r, trusted, k, thresh, detectStepsize));
//							System.out.println("\n"+toString(bs, r.bases.length));
//							System.out.println(new String(r.bases));
							for(int i=bs.nextClearBit(0); i<r.bases.length; i=bs.nextClearBit(i+1)){
								r.bases[i]='N';
								r.quality[i]=0;
							}
//							System.out.println(new String(r.bases));
//							System.out.println("used = "+String.format("%.3f%%",count.usedFraction()*100));
//							System.out.println("used = "+((KCountArray4)count).cellsUsed());
//							if(bs.length()<r.bases.length){r=null;}
						}
//						if(r!=null){addRead(r, count, k, mask, rcomp);}
						
						addReadSplit(r, count, k1, k2, mask1, mask2, gap, rcomp);
					}
					if(r2!=null){
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r2, trusted, k, thresh, detectStepsize) : 
								ErrorCorrect.detectTrusted(r2, trusted, k, thresh, detectStepsize));
							for(int i=bs.nextClearBit(0); i<r2.bases.length; i=bs.nextClearBit(i+1)){
								r2.bases[i]='N';
								r2.quality[i]=0;
							}
						}
						addReadSplit(r2, count, k1, k2, mask1, mask2, gap, rcomp);
					}

				}
				//System.err.println("returning list");
				cris.returnList(ln, ln.list.isEmpty());
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln, ln.list.isEmpty());
			if(verbose){System.err.println("Returned list");}
		}
		
		
		/**
		 * Hash a read's kmers into the KCountArray.
		 * Advanced mode processes paired reads together and sorts kmers to eliminate spurious duplicates.
		 * @param r1
		 * @param count
		 * @param k
		 * @param mask
		 * @param rcomp
		 */
		private final long[] addRead_Advanced(Read r1, final KCountArray count, final int k, final long mask, long[] array){
			assert(count.gap==0) : "Gapped: TODO";
			if(PREJOIN && r1.mate!=null && r1.insert()>0){
				r1.mate.reverseComplement();
				r1=r1.joinRead();
			}
			Read r2=r1.mate;
			int len1=r1.bases==null ? 0 : Tools.max(0, r1.bases.length-k+1);
			int len2=(r2==null || r2.bases==null) ? 0 : Tools.max(0, r2.bases.length-k+1);
			int len=len1+len2;
			if(len<1){return array;}
			if(array==null || array.length!=len){array=new long[len];}
			Arrays.fill(array, -1);
			fillKmerArray(r1, k, mask, array, 0, len1);
			if(r2!=null){fillKmerArray(r2, k, mask, array, len1, len);}
			if(KEEP_DUPLICATE_KMERS){
				for(long kmer : array){
					if(kmer!=-1){
						keysCountedLocal++;
						count.increment(kmer);
					}
				}
			}else{
				Arrays.sort(array);
				long prev=-1;
				for(int i=0; i<array.length; i++){
					long kmer=array[i];
					if(kmer!=prev){
						keysCountedLocal++;
						count.increment(kmer);
						prev=kmer;
					}
				}
			}
			return array;
		}
		
		private final void fillKmerArray(Read r, final int k, final long mask, final long[] array, final int start, final int stop){
			if(k>31){
				fillKmerArrayLong(r, k, array, start, stop);
				return;
			}
			assert(count.gap==0);
			assert(k<32);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(array!=null);
			
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			
			if(bases==null || bases.length<k+count.gap){return;}
			
			assert(kmersamplerate<2) : "samplerate: TODO";
			
			for(int pass=0; pass<2; pass++){
				int len=0;
				int idx=(pass==0 ? start-k+1 : stop+k-2);
				long kmer=0;
				float prob=1;
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					int x=AminoAcid.baseToNumber[b];

					byte q;
					if(quals==null){
						q=50;
					}else{
						q=quals[i];
						prob=prob*align2.QualityTools.PROB_CORRECT[q];
						if(len>k){
							byte oldq=quals[i-k];
							prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
						}
					}

					if(x<0 || q<minQuality){
						len=0;
						kmer=0;
						prob=1;
					}else if(prob<minProb){
						//do nothing
					}else{
						kmer=((kmer<<2)|x)&mask;
						len++;
						if(len>=k && (kmersamplerate<2 || i%kmersamplerate==0)){
							//						System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
							
//							assert(array[idx]==-1 || array[idx]==AminoAcid.reverseComplementBinaryFast(kmer, k)) : 
//								"\npass="+pass+", start="+start+", stop="+stop+", i="+i+", idx="+idx+"\n"+
//								"array[idx]="+array[idx]+", kmer="+kmer+", rcomp="+AminoAcid.reverseComplementBinaryFast(kmer, k)+
//								"\n"+Arrays.toString(array); //TODO: Remove slow assertion
//							System.out.println("\npass="+pass+", i="+i+", idx="+idx+"\nComparing "+kmer+"\tto  \t"+array[idx]);
							array[idx]=Tools.max(array[idx], kmer);
//							System.out.println("Kept      "+array[idx]+"\trcomp\t"+AminoAcid.reverseComplementBinaryFast(kmer, k));
						}
					}
					if(pass==0){idx++;}else{idx--;}
				}
//				System.out.println(Arrays.toString(array));
				r.reverseComplement();
			}
		}
		
		private final void fillKmerArrayLong(Read r, final int k, final long[] array, final int start, final int stop){
			assert(k>31);
			assert(count.gap==0);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(array!=null);
			
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			
			final int tailshift=k%32;
			final int tailshiftbits=tailshift*2;
			final long mask=Long.MAX_VALUE;
			
			if(bases==null || bases.length<k+count.gap){return;}
			
			assert(kmersamplerate<2) : "samplerate: TODO";
			
			for(int pass=0; pass<2; pass++){
				int len=0;
				int idx=(pass==0 ? start-k+1 : stop+k-2);
				long kmer=0;
				float prob=1;
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					int x=AminoAcid.baseToNumber[b];
					
					byte q;
					if(quals==null){
						q=50;
					}else{
						q=quals[i];
						prob=prob*align2.QualityTools.PROB_CORRECT[q];
						if(len>k){
							byte oldq=quals[i-k];
							prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
						}
					}

					if(x<0 || q<minQuality){
						len=0;
						kmer=0;
						prob=1;
					}else if(prob<minProb){
						//do nothing
					}else{
						kmer=Long.rotateLeft(kmer, 2);
						kmer=kmer^x;
						len++;

						if(len>=k){
							if(len>k){
								long x2=AminoAcid.baseToNumber[bases[i-k]];
								kmer=kmer^(x2<<tailshiftbits);
							}
							if((kmersamplerate<2 || i%kmersamplerate==0)){
								array[idx]=mask&Tools.max(array[idx], kmer);
							}
						}
					}
					if(pass==0){idx++;}else{idx--;}
				}
//				System.out.println(Arrays.toString(array));
				r.reverseComplement();
			}
		}
		
		private final void addReadSplit(Read r, final KCountArray count, final int k1, final int k2, final long mask1, final long mask2, final int gap, boolean rcomp){
			assert(false) : "TODO";
			if(PREJOIN && r.mate!=null && r.insert()>0){
				if(verbose){System.err.println("Prejoining "+r.numericID+" at "+r.insert());}
				r.mate.reverseComplement();
				r=r.joinRead();
			}
			
			int len=0;
			int shift=k2*2;
			long kmer1=0;
			long kmer2=0;
			float prob=1;
			byte[] bases=r.bases;
			byte[] quals=r.quality;
			
			assert(kmer1>=kmer2);
			
//			assert(false) : k1+", "+k2+", "+mask1+", "+mask2+", "+gap;
			
			if(verbose){System.err.println("Hashing read "+r.numericID+"; loop limits "+(k1+gap)+"-"+(bases.length));}
			for(int i=0, j=i+k1+gap; j<bases.length; i++, j++){
				int x1=AminoAcid.baseToNumber[bases[i]];
				int x2=AminoAcid.baseToNumber[bases[j]];
				
				byte q1, q2;
				if(quals==null){
					q1=50;
					q2=50;
				}else{
					q1=quals[i];
					q2=quals[j];
					prob=prob*align2.QualityTools.PROB_CORRECT[q1]*align2.QualityTools.PROB_CORRECT[q2];
					if(len>k){
						byte oldq1=quals[i-k1];
						byte oldq2=quals[j-k2];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq1]*align2.QualityTools.PROB_CORRECT_INVERSE[oldq2];
					}
				}
				
				if(x1<0 || x2<0 || q1<minQuality || q2<minQuality){
					len=0;
					kmer1=0;
					kmer2=0;
					prob=1;
				}else if(prob<minProb){
					//do nothing
				}else{
					kmer1=((kmer1<<2)|x1)&mask1;
					kmer2=((kmer2<<2)|x2)&mask2;
					len++;
					if(len>=k1){
						
						keysCountedLocal++;
//						System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
						
						long key=(kmer1<<shift)|kmer2;
//						System.err.println(Long.toHexString(key));
						
						if(verbose){System.err.println("Hashing key "+Long.toHexString(key)+" at length "+len);}
						count.increment(key);
//						count.increment(kmer);
						
						
//						System.out.println(" -> "+count.read(kmer));
//						System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//						array[(int)kmer]++;
//						System.out.println(" -> "+array[(int)kmer]+"\n");
//						assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
					}
				}
			}
			if(rcomp){
				r.reverseComplement();
				addReadSplit(r, count, k1, k2, mask1, mask2, gap, false);
			}
		}

		private final ConcurrentReadStreamInterface cris;
		private final int k;
		private final boolean rcomp; 
		private final KCountArray count;
		private final KCountArray trusted;
		private final int thresh;
		private final int detectStepsize;
		private final boolean conservative;
		private long keysCountedLocal=0;
		private long readsProcessedLocal=0;
	}
	
	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}
	
	public static byte minQuality=6;
	public static long readsProcessed=0;
	public static long maxReads=-1;
	public static int kmersamplerate=1;
	public static int readsamplerate=1;
	public static int BUFFERLEN=500;
	
	public static float minProb=0.5f;

	public static long keysCounted=0;
	
	public static int THREADS=Shared.THREADS;
	public static boolean verbose=false;
	public static boolean PREJOIN=false;
	public static boolean CANONICAL=false;
	public static boolean KEEP_DUPLICATE_KMERS=false;
	
}
