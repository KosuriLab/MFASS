package kmer;

import java.util.ArrayList;
import java.util.BitSet;

import jgi.ErrorCorrect;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;

import align2.ListNum;
import dna.AminoAcid;
import dna.Timer;
import fileIO.FileFormat;
import fileIO.TextFile;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 *
 */
public class KmerCount6 {
	
public static void main(String[] args){
		
		Timer t=new Timer();
		t.start();
		
		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		int k=14;
		int cbits=16;
		int gap=0;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=(split.length>1 ? split[1] : "true");

			if(a.equals("null")){
				// do nothing
			}else if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.startsWith("gap")){
				gap=Integer.parseInt(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		KCountArray count=null;
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			FastaReadInputStream.TARGET_READ_LEN=300000000;
			FastaReadInputStream.MIN_READ_LEN=k;
		}else{
			FASTQ.PARSE_CUSTOM=false;
//			assert(false) : FASTQ.PARSE_CUSTOM;
		}
		
		count=countFastq(fname1, fname2, k, cbits, gap, true, null);
		
		
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
						if(len>=k){
							keysCounted++;
							count.increment(kmer);
						}
					}
				}
			}
		}
		return count;
	}
	
	public static KCountArray countFastq(String reads1, String reads2, int k, int cbits, int gap, boolean rcomp, KCountArray count){
		assert(k<32 && k>=1 && (count!=null || k<20));
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		
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
		System.err.println("Started cris");
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		countFastq(cris, k, rcomp, count);
		
		cris.close();
		if(verbose){System.err.println("Closed stream");}
		if(verbose){System.err.println("Processed "+readsProcessed+" reads.");}

		
		return count;
	}
	
	
	public static void countFastq(ConcurrentReadStreamInterface cris, int k, boolean rcomp, KCountArray count){
		assert(k<32 && k>=1 && (count!=null || k<20));
		
		assert(count!=null);
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		
		if(count.gap==0){
			final int kbits=2*k;
			final long mask=~((-1L)<<(kbits));
			
			
			while(reads!=null && reads.size()>0){
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					readsProcessed++;

					addRead(r, count, k, mask, rcomp);
					if(r.mate!=null){
						addRead(r.mate, count, k, mask, rcomp);
					}

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
					readsProcessed++;

					addReadSplit(r, count, k1, k2, mask1, mask2, gap, rcomp);
					if(r.mate!=null){
						addReadSplit(r.mate, count, k1, k2, mask1, mask2, gap, rcomp);
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
		cris.returnList(ln, ln.list.isEmpty());
		if(verbose){System.err.println("Returned list");}
	}
	
	

	
	public static KCountArray countFastq(final String reads1, final String reads2, final int k, final int cbits, final boolean rcomp, 
			KCountArray count, final KCountArray trusted, final long maxReads, final int thresh, final int detectStepsize, final boolean conservative){
		
		assert(k<32 && k>=1 && (count!=null || k<20));
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		
//		System.out.println("k="+k+", kbits="+kbits+", mask="+Long.toHexString(mask)+", thresh="+thresh);
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);
		
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
		System.err.println("Started cris");
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		countFastq(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);
		
		cris.close();
		if(verbose){System.err.println("Closed stream");}
		
//		System.out.println("*** after ***");
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);
		
		return count;
	}
	
	

	
	public static void countFastq(final ConcurrentReadStreamInterface cris, final int k, final boolean rcomp, 
			final KCountArray count, final KCountArray trusted, final int thresh, final int detectStepsize, final boolean conservative){
		
		assert(k<32 && k>=1 && (count!=null || k<20));
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		
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
//						System.out.println("\n"+toString(bs, r.bases.length));
//						System.out.println(new String(r.bases));
						for(int i=bs.nextClearBit(0); i<r.bases.length; i=bs.nextClearBit(i+1)){
							r.bases[i]='N';
							r.quality[i]=0;
						}
//						System.out.println(new String(r.bases));
//						System.out.println("used = "+String.format("%.3f%%",count.usedFraction()*100));
//						System.out.println("used = "+((KCountArray4)count).cellsUsed());
//						if(bs.length()<r.bases.length){r=null;}
					}
//					if(r!=null){addRead(r, count, k, mask, rcomp);}
					addRead(r, count, k, mask, rcomp);
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
					addRead(r2, count, k, mask, rcomp);
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
	
	
	
	public static void addRead(final Read r, final KCountArray count, final int k, final long mask, boolean rcomp){
		int len=0;
		long kmer=0;
		byte[] bases=r.bases;
		byte[] quals=r.quality;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0 || (quals!=null && quals[i]<minQuality)){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				if(len>=k){
					keysCounted++;
//					System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
					count.increment(kmer);
//					System.out.println(" -> "+count.read(kmer));
//					System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//					array[(int)kmer]++;
//					System.out.println(" -> "+array[(int)kmer]+"\n");
//					assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
				}
			}
		}
		if(rcomp){
			r.reverseComplement();
			addRead(r, count, k, mask, false);
		}
	}
	
	public static void addReadSplit(final Read r, final KCountArray count, final int k1, final int k2, final long mask1, final long mask2, final int gap, boolean rcomp){
		int len=0;
		int shift=k2*2;
		long kmer1=0;
		long kmer2=0;
		byte[] bases=r.bases;
		byte[] quals=r.quality;
		
		assert(kmer1>=kmer2);
		
//		assert(false) : k1+", "+k2+", "+mask1+", "+mask2+", "+gap;
		
		for(int i=0, j=i+k1+gap; j<bases.length; i++, j++){
			int x1=AminoAcid.baseToNumber[bases[i]];
			int x2=AminoAcid.baseToNumber[bases[j]];
			if(x1<0 || x2<0 || (quals!=null && (quals[i]<minQuality || quals[j]<minQuality))){
				len=0;
				kmer1=0;
				kmer2=0;
			}else{
				kmer1=((kmer1<<2)|x1)&mask1;
				kmer2=((kmer2<<2)|x2)&mask2;
				len++;
				if(len>=k1){
					keysCounted++;
//					System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
					
					long key=(kmer1<<shift)|kmer2;
//					System.err.println(Long.toHexString(key));
					count.increment(key);
//					System.out.println(" -> "+count.read(kmer));
//					System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//					array[(int)kmer]++;
//					System.out.println(" -> "+array[(int)kmer]+"\n");
//					assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
				}
			}
		}
		if(rcomp){
			r.reverseComplement();
			addReadSplit(r, count, k1, k2, mask1, mask2, gap, false);
		}
	}
	
	public static void addReadSplit(final byte[] bases, final KCountArray count, final int k1, final int k2, final long mask1, final long mask2, final int gap, boolean rcomp){
		int len=0;
		int shift=k2*2;
		long kmer1=0;
		long kmer2=0;
		byte[] quals=null;
		
		assert(kmer1>=kmer2);
		
//		assert(false) : k1+", "+k2+", "+mask1+", "+mask2+", "+gap;
		
		for(int i=0, j=i+k1+gap; j<bases.length; i++, j++){
			int x1=AminoAcid.baseToNumber[bases[i]];
			int x2=AminoAcid.baseToNumber[bases[j]];
			if(x1<0 || x2<0 || (quals!=null && (quals[i]<minQuality || quals[j]<minQuality))){
				len=0;
				kmer1=0;
				kmer2=0;
			}else{
				kmer1=((kmer1<<2)|x1)&mask1;
				kmer2=((kmer2<<2)|x2)&mask2;
				len++;
				if(len>=k1){
					keysCounted++;
//					System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
					
					long key=(kmer1<<shift)|kmer2;
					System.out.println(Long.toHexString(kmer1));
					System.out.println(Long.toHexString(kmer2));
					System.out.println(Long.toHexString(key));
					count.increment(key);
//					System.out.println(" -> "+count.read(kmer));
//					System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//					array[(int)kmer]++;
//					System.out.println(" -> "+array[(int)kmer]+"\n");
//					assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
				}
			}
		}
		if(rcomp){
			AminoAcid.reverseComplementBasesInPlace(bases);
			addReadSplit(bases, count, k1, k2, mask1, mask2, gap, false);
		}
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
	
	private static ThreadLocal<long[]> buffer=new ThreadLocal<long[]>();
	
	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}
	
	public static boolean verbose=false;
	public static byte minQuality=9;
	public static long readsProcessed=0;
	public static long maxReads=-1;
	public static int bufferlen=1000;

	public static long keysCounted=0;
	
}
