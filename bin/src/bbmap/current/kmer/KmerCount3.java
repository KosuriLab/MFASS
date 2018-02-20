package kmer;

import java.util.ArrayList;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FastaReadInputStream;
import stream.FastqReadInputStream;
import stream.Read;

import align2.ListNum;
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
public class KmerCount3 {
	
public static void main(String[] args){
		
		Timer t=new Timer();
		t.start();
		
		String fname1=args[0];
		String fname2=(args.length>3 || args[1].contains(".") ? args[1] : null);
		int k=Integer.parseInt(args[args.length-2]);
		int cbits=Integer.parseInt(args[args.length-1]);
		
		KCountArray2 count=null;
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			FastaReadInputStream.TARGET_READ_LEN=300000000;
			FastaReadInputStream.MIN_READ_LEN=k;
		}
		count=countFastq(fname1, fname2, k, cbits);
		
		
		t.stop();
		System.out.println("Finished counting; time = "+t);
		
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
	}
	
	public static KCountArray2 countFasta(String fname, int k, int cbits){
		assert(k>=1 && k<20);
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		final long cells=mask+1;
		if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
		final KCountArray2 count=new KCountArray2(cells, cbits);
		
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
							count.increment(kmer, 1);
						}
					}
				}
			}
		}
		return count;
	}
	
	public static KCountArray2 countFastq(String reads1, String reads2, int k, int cbits){
		assert(k>=1 && k<20);
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		final long cells=mask+1;
		if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
		final KCountArray2 count=new KCountArray2(cells, cbits);
		
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
		System.err.println("Paired: "+paired);
		
		long kmer=0; //current kmer
		int len=0;  //distance since last contig start or ambiguous base
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(paired==(r.mate!=null));
			}
			
			while(reads!=null && reads.size()>0){
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					readsProcessed++;
					
					len=0;
					kmer=0;
					byte[] bases=r.bases;
					byte[] quals=r.quality;
					for(int i=0; i<bases.length; i++){
						byte b=bases[i];
						int x=AminoAcid.baseToNumber[b];
						if(x<0 || quals[i]<minQuality){
							len=0;
							kmer=0;
						}else{
							kmer=((kmer<<2)|x)&mask;
							len++;
							if(len>=k){
//								System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
								count.increment(kmer, 1);
//								System.out.println(" -> "+count.read(kmer));
//								System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//								array[(int)kmer]++;
//								System.out.println(" -> "+array[(int)kmer]+"\n");
//								assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
							}
						}
					}
					
					if(r.mate!=null){
						len=0;
						kmer=0;
						bases=r.mate.bases;
						quals=r.mate.quality;
						for(int i=0; i<bases.length; i++){
							byte b=bases[i];
							int x=AminoAcid.baseToNumber[b];
							if(x<0 || quals[i]<minQuality){
								len=0;
								kmer=0;
							}else{
								kmer=((kmer<<2)|x)&mask;
								len++;
								if(len>=k){
									count.increment(kmer, 1);
								}
							}
						}
					}
					
				}
				//System.err.println("returning list");
				cris.returnList(ln, ln.list.isEmpty());
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			System.err.println("Finished reading");
			cris.returnList(ln, ln.list.isEmpty());
			System.err.println("Returned list");
			ReadWrite.closeStream(cris);
			System.err.println("Closed stream");
			System.err.println("Processed "+readsProcessed+" reads.");
		}
		
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
	
	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}
	
	public static boolean verbose=false;
	public static byte minQuality=5;
	public static long readsProcessed=0;
	public static long maxReads=10000000000L;
	
}
