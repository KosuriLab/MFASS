package stream;

import java.util.ArrayList;
import java.util.Arrays;

import align2.Shared;
import align2.Tools;

import dna.Data;
import fileIO.TextFile;

public class FastaStream {
	
	public static void main(String[] args){
		FastaStream qs=new FastaStream(args[0]);
		for(int i=0; i<9000; i++){
			byte[][] next=qs.next();
			System.out.println(new String(next[0]));
			System.out.println(new String(next[1]));
			System.out.println();
		}
	}
	
	public FastaStream(String fname){
		
		if(!fileIO.FileFormat.hasFastaExtension(fname) && !fname.startsWith("stdin")){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+fname);
		}
		tf=new TextFile(fname, false, false);
		
		if(SKIP_INITIAL_READS){
			while(consumed<INITIAL_READS_TO_SKIP){nextBlock();}
		}
	}
	
	public boolean hasMore() {
		//System.err.println("Calling hasMore()");
		if(buffer==null || next>=buffer.length){
//			System.err.println("Attempting to fill buffer.");
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.length);
	}
	
	public byte[][] next() {
		//System.err.println("Calling next()");
		if(!hasMore()){return null;}
		byte[][] r=buffer[next];
		buffer[next]=null;
		next++;
		consumed++;
		return r;
	}
	
	public synchronized byte[][][] nextBlock() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.length){fillBuffer();}
		byte[][][] r=buffer;
		buffer=null;
		if(r!=null && r.length==0){r=null;}
		consumed+=(r==null ? 0 : r.length);
		return r;
	}
	
	private synchronized void fillBuffer(){
		//System.err.println("Calling fillBuffer()");
		
		assert(buffer==null || next>=buffer.length);
		
		buffer=null;
		next=0;
		
		buffer=toPairs(tf, BUF_LEN);
		if(buffer.length<BUF_LEN){tf.close();}
		
		generated+=buffer.length;
	}
	
	private static byte[][][] toPairs(TextFile tf, int maxToReturn){
		//System.err.println("Calling toPairs()");
		String s=null;
		ArrayList<byte[][]> list=new ArrayList<byte[][]>(Tools.min(16384, maxToReturn));
		
		String[] pair=new String[2];
		
		int cntr=0;
		int added=0;
		final byte dot='.';
		final byte N='N';
		
		for(s=tf.nextLine(); s!=null && added<maxToReturn; s=tf.nextLine()){
			while(s!=null && s.startsWith("#")){s=tf.nextLine();}
			if(s==null){break;}
//			System.err.println("Processing "+s);
			
			pair[cntr]=s;
			cntr++;
			if(cntr==2){
				assert(pair[0].startsWith(">"));
				
				//TODO Note: These assertions are only for SOLiD colorspace.
				assert(pair[1].charAt(0)=='T' || pair[1].charAt(0)=='G' || 
						pair[1].charAt(0)=='.' || Character.isDigit(pair[1].charAt(0))) : pair[1];
				assert(pair[1].charAt(1)=='.' || Character.isDigit(pair[1].charAt(1))) : pair[1];
				
				byte[][] fixed=new byte[2][];
				fixed[0]=Arrays.copyOfRange(pair[0].getBytes(), 1, pair[0].length());
				fixed[1]=pair[1].getBytes();
				for(int i=0; i<fixed[1].length; i++){
					if(fixed[1][i]==dot){
						fixed[1][i]=N;
					}else if(Character.isDigit(fixed[1][i])){
						fixed[1][i]-='0';
					}
				}
				
				
				list.add(fixed);
				cntr=0;
				added++;
				
				if(added>=maxToReturn){break;}
				
//				System.out.println(r.chrom+", "+r.strand+", "+r.loc);
//				assert(false);
			}
		}
		
//		for(int i=0; i<12000; i++){tf.nextLine();}
		
		assert(list.size()<=maxToReturn);
		return list.toArray(new byte[0][][]);
	}
	
	public void close(){
		tf.close();
	}

	private byte[][][] buffer;
	
	private int next=0;

	private final TextFile tf;

	public static final int BUF_LEN=Shared.READ_BUFFER_LENGTH;

	public long generated=0;
	public long consumed=0;
	
	public static boolean SKIP_INITIAL_READS=QualStream.SKIP_INITIAL_READS;
	public static final int INITIAL_READS_TO_SKIP=QualStream.INITIAL_READS_TO_SKIP;

}
