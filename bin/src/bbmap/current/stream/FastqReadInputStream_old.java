package stream;

import java.util.ArrayList;

import align2.Shared;

import dna.Data;
import fileIO.FileFormat;
import fileIO.TextFile;

public class FastqReadInputStream_old extends ReadInputStream {
	
	public static void main(String[] args){
		
		FASTQ.PARSE_CUSTOM=false;
		
		FastqReadInputStream_old fris=new FastqReadInputStream_old(args[0], false, true);
		
		Read r=fris.next();
		System.out.println(r.toText(false));
		
	}
	
	public FastqReadInputStream_old(String fname, boolean colorspace_, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess_, false), colorspace_);
	}

	
	public FastqReadInputStream_old(FileFormat ff, boolean colorspace_){

		colorspace=colorspace_;
		
		stdin=ff.stdio();
		if(!ff.fastq()){
			System.err.println("Warning: Did not find expected fastq file extension for filename "+ff.name());
		}
		
		if(FASTQ.PARSE_CUSTOM){
			try {
				String s[]=ff.name().split("_");
//				maxSnps=toNumber(s[3]);
//				maxInss=toNumber(s[4]);
//				maxDels=toNumber(s[5]);
//				maxSubs=toNumber(s[6]);
				
				s=s[7].split("\\.");
				
				s=s[0].split("-");
				
//				minChrom=Gene.toChromosome(s[0]);
//				maxChrom=Gene.toChromosome(s[1]);

			} catch (Exception e) {
				// TODO Auto-generated catch block
				//			e.printStackTrace();
				if(Data.WINDOWS){System.out.println("Note: Filename indicates non-synthetic data, but FASTQ.PARSE_CUSTOM="+FASTQ.PARSE_CUSTOM);}
			}
		}
		
		tf=new TextFile(ff.name(), false, false);
		interleaved=((tf.is==System.in || stdin) ? FASTQ.FORCE_INTERLEAVED : FASTQ.isInterleaved(tf.name));
		
	}

	@Override
	public void start() {
//		if(cris!=null){new Thread(cris).start();}
	}
	
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.length){
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.length);
	}

	@Override
	public Read next() {
		if(!hasMore()){return null;}
		Read r=buffer[next];
		buffer[next]=null;
		next++;
		consumed++;
		return r;
	}
	
	@Override
	public synchronized Read[] nextBlock() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.length){fillBuffer();}
		Read[] r=buffer;
		buffer=null;
		if(r!=null && r.length==0){r=null;}
		consumed+=(r==null ? 0 : r.length);
//		System.err.println(hashCode()+" produced "+r[0].numericID);
		return r;
	}
	
	@Override
	public synchronized ArrayList<Read> nextList() {
		return toList(nextBlock());
	}
	public final boolean preferArrays(){return true;}
	
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.length);
		
		buffer=null;
		next=0;
		
		buffer=FASTQ.toReads(tf, BUF_LEN, colorspace, nextReadID, interleaved);
		nextReadID+=buffer.length;
		if(buffer.length<BUF_LEN){tf.close();}
		
		generated+=buffer.length;
	}
	
	public boolean close(){
		return tf.close();
	}

	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		tf.reset();
	}

	@Override
	public boolean paired() {
		return interleaved;
	}

	private Read[] buffer=null;
	private int next=0;
	
	private final TextFile tf;
	private final boolean interleaved;

	private final int BUF_LEN=Shared.READ_BUFFER_LENGTH;

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public final boolean colorspace;
	public final boolean stdin;

}
