package stream;

import java.util.ArrayList;
import java.util.Arrays;

import align2.Shared;
import align2.Tools;

import dna.Data;
import fileIO.TextFile;

public class QualStream {
	
	public static void main(String[] args){
		QualStream qs=new QualStream(args[0]);
		for(int i=0; i<9000; i++){
			byte[][] next=qs.next();
			System.out.println(new String(next[0]));
			for(int j=0; j<next[1].length; j++){
				System.out.print(next[1][j]);
				System.out.print(" ");
			}
			System.out.println("\n");
		}
	}
	
	public QualStream(String fname){
		
		boolean flag=(fname.contains(".qual"));
		if(!flag){
			System.err.println("Warning: Did not find expected fastq file extension for filename "+fname);
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
		
		for(s=tf.nextLine(); s!=null && added<maxToReturn; s=tf.nextLine()){
			while(s!=null && s.startsWith("#")){s=tf.nextLine();}
			if(s==null){break;}
//			System.err.println("Processing "+s);
			
			pair[cntr]=s;
			cntr++;
			if(cntr==2){
				assert(pair[0].startsWith(">"));
				assert(pair[1].startsWith("-") || Character.isDigit(pair[1].charAt(0)));
				
//				int spaces=0;
//				for(int i=0; i<pair[1].length(); i++){
//					if(pair[1].charAt(i)==' '){
//						spaces++;
//					}
//				}
//				byte[] numbers=new byte[spaces+1];
				
				byte[] numbers=toNumbers(pair[1]);
				
				byte[][] fixed=new byte[2][];
				fixed[0]=Arrays.copyOfRange(pair[0].getBytes(), 1, pair[0].length());
				fixed[1]=numbers;
				
				
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
	
	/** TODO Speed up by removing the "split" call. */
	public static byte[] toNumbers(String s){
		//System.err.println("Calling toNumbers()");
		String[] split=s.split(" ");
		byte[] out=new byte[split.length];
		for(int i=0; i<split.length; i++){
//			out[i]=(byte)(Byte.parseByte(split[i])+34);
			out[i]=(byte)(Byte.parseByte(split[i]));
		}
		return out;
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
	
	public static boolean SKIP_INITIAL_READS=false;
	public static final int INITIAL_READS_TO_SKIP=5000000;

}
