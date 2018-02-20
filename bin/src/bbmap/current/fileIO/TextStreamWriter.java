package fileIO;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import stream.Read;

import dna.Data;



public class TextStreamWriter extends Thread {

//	public TextStreamWriter(String fname_, boolean overwrite_){
//		this(fname_, overwrite_, false, false, -1);
//	}
//	
//	public TextStreamWriter(String fname_, boolean overwrite_, boolean append_){
//		this(fname_, overwrite_, append_, false, -1);
//	}
	
	public TextStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_){
		this(fname_, overwrite_, append_, allowSubprocess_, 0);
	}
	
	public TextStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_, int format){
		this(FileFormat.testOutput(fname_, FileFormat.TEXT, format, 0, allowSubprocess_, overwrite_, true), append_);
	}
	
	public TextStreamWriter(FileFormat ff, boolean append_){
		FASTQ=ff.fastq() || ff.text();
		FASTA=ff.fasta();
		BREAD=ff.bread();
		SAM=ff.samOrBam();
		BAM=ff.bam();
		SITES=ff.sites();
		INFO=ff.attachment();
		OTHER=(!FASTQ && !FASTA && !BREAD && !SAM && !BAM && !SITES && !INFO);
		
		
		fname=ff.name();
		overwrite=ff.overwrite();
		allowSubprocess=ff.allowSubprocess();
		assert(!(overwrite&append_));
		assert(append_ || ff.canWrite()) : "File "+fname+" exists and overwrite=="+overwrite;
		if(append_ && !ff.raw()){throw new RuntimeException("Can't append to compressed files.");}
		
		if(!BAM || !Data.SAMTOOLS() || !Data.SH()){
			myOutstream=ReadWrite.getOutputStream(fname, append_, true, allowSubprocess);
		}else{
			myOutstream=ReadWrite.getOutputStreamFromProcess(fname, "samtools view -S -b -h - ", true);
		}
		myWriter=new PrintWriter(myOutstream);
		
		queue=new ArrayBlockingQueue<ArrayList<CharSequence>>(5);
		buffer=new ArrayList<CharSequence>(buffersize);
	}
	
	@Override
	public void run() {
		assert(open) : fname;
		synchronized(this){
			started=true;
			this.notify();
		}
		
		ArrayList<CharSequence> job=null;

		while(job==null){
			try {
				job=queue.take();
//				job.list=queue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		while(job!=null && job!=POISON2){
			if(!job.isEmpty()){
				for(final CharSequence cs : job){
					assert(cs!=POISON);
					myWriter.print(cs);
				}
			}
			
			job=null;
			while(job==null){
				try {
					job=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
//		assert(false);
		open=false;
		ReadWrite.finishWriting(myWriter, myOutstream, fname, allowSubprocess);
		synchronized(this){notifyAll();}
	}
	
	public void print(CharSequence cs){
//		System.err.println("Added line '"+cs+"'");
		assert(open) : cs;
		buffer.add(cs);
		bufferLen+=cs.length();
		if(buffer.size()>=buffersize || bufferLen>=maxBufferLen){
			addJob(buffer);
			buffer=new ArrayList<CharSequence>(buffersize);
			bufferLen=0;
		}
	}
	
	public void println(CharSequence cs){
		print(cs);
		print("\n");
	}
	
	public void println(Read r){
		assert(!OTHER);
		StringBuilder sb=(FASTQ ? r.toFastq() : FASTA ? r.toFasta() : SAM ? r.toSam() : 
			SITES ? r.toSites() : INFO ? r.toInfo() : r.toText(true)).append('\n');
		print(sb);
	}
	
	public void print(Read r){
		assert(!OTHER);
		StringBuilder sb=(FASTQ ? r.toFastq() : FASTA ? r.toFasta() : SAM ? r.toSam() : 
			SITES ? r.toSites() : INFO ? r.toInfo() : r.toText(true));
		print(sb);
	}

	public synchronized void poison(){
		if(!open){return;}
		addJob(buffer);
		buffer=null;
//		System.err.println("Poisoned!");
//		assert(false);
		
		//Don't allow thread to shut down before it has started
		while(this.getState()==Thread.State.NEW || !started){
			try {
				this.wait(20);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
//		assert(false) : open+", "+this.getState()+", "+started;
		open=false;
		addJob(POISON2);
	}
	
	public void waitForFinish(){
		while(this.getState()!=Thread.State.TERMINATED){
			try {
				this.join(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	public void poisonAndWait(){
		poison();
		waitForFinish();
	}
	
	//TODO Why is this synchronized?
	public synchronized void addJob(ArrayList<CharSequence> j){
//		System.err.println("Got job "+(j.list==null ? "null" : j.list.size()));
		boolean success=false;
		while(!success){
			try {
				queue.put(j);
				success=true;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				assert(!queue.contains(j)); //Hopefully it was not added.
			}
		}
	}
	
	private ArrayList<CharSequence> buffer;
	
	public int buffersize=100;
	public int maxBufferLen=60000;
	private int bufferLen=0;
	public final boolean overwrite;
	public final boolean allowSubprocess;
	public final String fname;
	private final OutputStream myOutstream;
	private final PrintWriter myWriter;
	private final ArrayBlockingQueue<ArrayList<CharSequence>> queue;
	private boolean open=true;
	private volatile boolean started=false;
	
	/** TODO */
	public boolean errorState=false;
	
	private final boolean BAM;
	private final boolean SAM;
	private final boolean FASTQ;
	private final boolean FASTA;
	private final boolean BREAD;
	private final boolean SITES;
	private final boolean INFO;
	private final boolean OTHER;

	private static final String POISON=new String("POISON_TextStreamWriter");
	private static final ArrayList<CharSequence> POISON2=new ArrayList<CharSequence>(1);
	
}
