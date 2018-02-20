package stream;

import java.util.ArrayList;

import align2.Shared;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;

public class FastqReadInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		
		FASTQ.PARSE_CUSTOM=false;
		
		FastqReadInputStream fris=new FastqReadInputStream(args[0], false, true);
		
		Read r=fris.next();
		System.out.println(r.toText(false));
		
	}
	
	public FastqReadInputStream(String fname, boolean colorspace_, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess_, false), colorspace_);
	}

	
	public FastqReadInputStream(FileFormat ff, boolean colorspace_){
		if(verbose){System.err.println("FastqReadInputStream("+ff+", "+colorspace_+")");}

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
				
//				s=s[8].split("\\.");
//				
//				s=s[0].split("-");
				
				if(s.length!=8 && s.length!=9){
					if(Data.WINDOWS){System.err.println("Note: Filename indicates non-synthetic data, but FASTQ.PARSE_CUSTOM="+FASTQ.PARSE_CUSTOM);}
				}
				
//				minChrom=Gene.toChromosome(s[0]);
//				maxChrom=Gene.toChromosome(s[1]);

			} catch (Exception e) {
				// TODO Auto-generated catch block
				//			e.printStackTrace();
				if(Data.WINDOWS){System.err.println("Note: Filename indicates non-synthetic data, but FASTQ.PARSE_CUSTOM="+FASTQ.PARSE_CUSTOM);}
			}
		}
		
		tf=ByteFile.makeByteFile(ff, false);
		interleaved=((tf.is()==System.in || stdin) ? FASTQ.FORCE_INTERLEAVED : FASTQ.isInterleaved(tf.name()));
//		assert(false) : interleaved;
	}

	@Override
	public void start() {
//		if(cris!=null){new Thread(cris).start();}
	}
	
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.size()){
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}

	@Override
	public Read next() {
		if(!hasMore()){return null;}
		Read r=buffer.set(next, null);
		next++;
		consumed++;
		return r;
	}
	
	@Override
	public synchronized Read[] nextBlock() {
		ArrayList<Read> x=nextList();
		if(x==null){return null;}
		return x.toArray(new Read[x.size()]);
	}
	
	@Override
	public synchronized ArrayList<Read> nextList() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> list=buffer;
		buffer=null;
		if(list!=null && list.size()==0){list=null;}
		consumed+=(list==null ? 0 : list.size());
//		System.err.println(hashCode()+" produced "+r[0].numericID);
		return list;
	}
	public final boolean preferArrays(){return false;}
	
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=FASTQ.toReadList(tf, BUF_LEN, colorspace, nextReadID, interleaved);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){tf.close();}
		
		generated+=bsize;
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in FastqReadInputStream.");
			}
		}
	}
	
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
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
	public boolean paired() {return interleaved;}
	
	/** Return true if this stream has detected an error */
	public boolean errorState(){return errorState || FASTQ.errorState();}

	private ArrayList<Read> buffer=null;
	private int next=0;
	
	private final ByteFile tf;
	private final boolean interleaved;

	private final int BUF_LEN=Shared.READ_BUFFER_LENGTH;
	private final long MAX_DATA=Shared.READ_BUFFER_MAX_DATA; //TODO - lot of work for unlikely case of super-long fastq reads.  Must be disabled for paired-ends.

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public final boolean colorspace;
	public final boolean stdin;
	public static boolean verbose=false;

}
