package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import dna.Data;

import fileIO.ReadWrite;
import fileIO.FileFormat;

public abstract class ReadStreamWriter extends Thread {

//	public ReadStreamWriter(String fname_, String qfname_, boolean read1_, int bufferSize, 
//			boolean outputSamFile, boolean outputBamFile, boolean fastq, boolean fasta, boolean sitesOnly, boolean attachment, boolean stdout,
//			boolean useSharedHeader, boolean makeWriter, boolean buffered){
//		this(fname_, qfname_, read1_, bufferSize, 
//				outputSamFile, outputBamFile, fastq, fasta, sitesOnly, attachment, stdout,
//				useSharedHeader, makeWriter, buffered, true);
//		
//	}

	protected ReadStreamWriter(String fname_, String qfname_, boolean read1_, int bufferSize, 
			boolean outputSamFile, boolean outputBamFile, boolean fastq, boolean fasta, boolean sitesOnly, boolean attachment, boolean stdout,
			boolean useSharedHeader, boolean makeWriter, boolean buffered, boolean allowSubprocess_){
		fname=fname_;
		qfname=qfname_;
		read1=read1_;
		OUTPUT_SAM=(outputSamFile || outputBamFile);
		OUTPUT_BAM=outputBamFile;
		OUTPUT_FASTQ=fastq;
		OUTPUT_FASTA=fasta;
		OUTPUT_ATTACHMENT=attachment;
		OUTPUT_STANDARD_OUT=stdout;
		SITES_ONLY=sitesOnly;
		allowSubprocess=allowSubprocess_;
		assert(((OUTPUT_SAM ? 1 : 0)+(OUTPUT_FASTQ ? 1 : 0)+(OUTPUT_FASTA ? 1 : 0)+(OUTPUT_ATTACHMENT ? 1 : 0)+(SITES_ONLY ? 1 : 0))<=1) : 
			OUTPUT_SAM+", "+SITES_ONLY+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT;
//		assert(fname==null || (fname.contains(".sam") || fname.contains(".bam"))==OUTPUT_SAM) : "Outfile name and sam output mode flag disagree: "+fname;
		assert(read1 || !OUTPUT_SAM) : "Attempting to output paired reads to different sam files.";
//		System.err.println("Called ReadStreamWriter with fname="+fname+", qfname="+qfname);
		
		if(qfname==null){
			myQOutstream=null;
			myQWriter=null;
		}else{
			myQOutstream=ReadWrite.getOutputStream(qfname, false, buffered, allowSubprocess);
			myQWriter=(makeWriter ? new PrintWriter(myQOutstream) : null);
		}
//		System.err.println("myQWriter="+myQWriter);
		
		if(fname==null && !OUTPUT_STANDARD_OUT){
			myOutstream=null;
			myWriter=null;
		}else{
//			if(OUTPUT_STANDARD_OUT){myOutstream=System.out;}
//			else 
			if(!OUTPUT_BAM || !Data.SAMTOOLS() || !Data.SH()){
				myOutstream=ReadWrite.getOutputStream(fname, false, buffered, allowSubprocess);
			}else{
				myOutstream=ReadWrite.getOutputStreamFromProcess(fname, "samtools view -S -b -h - ", true);
			}
			myWriter=(makeWriter ? new PrintWriter(myOutstream) : null);
			
			if(HEADER!=null){
				if(makeWriter){
					myWriter.println(HEADER);
				}else{
					ByteBuilder bb=new ByteBuilder(HEADER);
					try {
						if(bb.length>0){myOutstream.write(bb.array, 0, bb.length);}
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}else if(OUTPUT_SAM){
				if(useSharedHeader){
					ArrayList<byte[]> list=SamReadInputStream.getSharedHeader(true);
					if(list==null){
						System.err.println("Header was null.");
					}else{
						try {
							for(byte[] line : list){
								myOutstream.write(line);
								myOutstream.write('\n');
								//myWriter.println(new String(line));
							}
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}else{
					if(makeWriter){
						myWriter.println(SamLine.header0());
						int a=(MINCHROM==-1 ? 1 : MINCHROM);
						int b=(MAXCHROM==-1 ? Data.numChroms : MAXCHROM);
						for(int chrom=a; chrom<=b; chrom++){
							//					myWriter.print(SamLine.header1(chrom, chrom));
							SamLine.printHeader1(chrom, chrom, myWriter);
						}
						myWriter.println(SamLine.header2());
					}else{
						ByteBuilder bb=new ByteBuilder(4096);
						SamLine.header0B(bb);
						bb.append('\n');
						int a=(MINCHROM==-1 ? 1 : MINCHROM);
						int b=(MAXCHROM==-1 ? Data.numChroms : MAXCHROM);
						for(int chrom=a; chrom<=b; chrom++){
							SamLine.printHeader1B(chrom, chrom, bb, myOutstream);
						}
						SamLine.header2B(bb);
						bb.append('\n');
						
						
						try {
							if(bb.length>0){myOutstream.write(bb.array, 0, bb.length);}
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
			}else if(!OUTPUT_SAM && !OUTPUT_FASTQ && !OUTPUT_FASTA && !SITES_ONLY && !OUTPUT_ATTACHMENT){
//				myWriter.println("#"+Read.header());
			}
		}
		
		assert(bufferSize>=1);
		queue=new ArrayBlockingQueue<Job>(bufferSize);
	}
	
	protected ReadStreamWriter(FileFormat ff, String qfname_, boolean read1_, int bufferSize, CharSequence header, boolean makeWriter, boolean buffered, boolean useSharedHeader){
//		assert(false) : useSharedHeader+", "+header;
		assert(ff!=null);
		assert(ff.write()) : "FileFormat is not in read mode for "+ff.name();
		
		assert(!ff.text() && !ff.unknownFormat()) : "Unknown format for "+ff;
		OUTPUT_FASTQ=ff.fastq();
		OUTPUT_FASTA=ff.fasta();
//		boolean bread=(ext==TestFormat.txt);
		OUTPUT_SAM=ff.samOrBam();
		OUTPUT_BAM=ff.bam();
		OUTPUT_ATTACHMENT=ff.attachment();
		SITES_ONLY=ff.sites();
		OUTPUT_STANDARD_OUT=ff.stdio();
		assert(((OUTPUT_SAM ? 1 : 0)+(OUTPUT_FASTQ ? 1 : 0)+(OUTPUT_FASTA ? 1 : 0)+(OUTPUT_ATTACHMENT ? 1 : 0)+(SITES_ONLY ? 1 : 0))<=1) : 
			OUTPUT_SAM+", "+SITES_ONLY+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT;
		
		fname=ff.name();
		qfname=qfname_;
		read1=read1_;
		allowSubprocess=ff.allowSubprocess();
//		assert(fname==null || (fname.contains(".sam") || fname.contains(".bam"))==OUTPUT_SAM) : "Outfile name and sam output mode flag disagree: "+fname;
		assert(read1 || !OUTPUT_SAM) : "Attempting to output paired reads to different sam files.";
		
		if(qfname==null){
			myQOutstream=null;
			myQWriter=null;
		}else{
			myQOutstream=ReadWrite.getOutputStream(fname, false, buffered, allowSubprocess);
			myQWriter=(makeWriter ? new PrintWriter(myQOutstream) : null);
		}
		
		if(header==null){header=HEADER;} //new line; test.
		
		
		if(fname==null && !OUTPUT_STANDARD_OUT){
			myOutstream=null;
			myWriter=null;
		}else{
			if(OUTPUT_STANDARD_OUT){myOutstream=System.out;}
			else if(!OUTPUT_BAM || !Data.SAMTOOLS() || !Data.SH()){
				myOutstream=ReadWrite.getOutputStream(fname, false, buffered, allowSubprocess);
			}else{
				if(!allowSubprocess){System.err.println("Warning! Spawning a samtools process when allowSubprocess="+allowSubprocess);}
				myOutstream=ReadWrite.getOutputStreamFromProcess(fname, "samtools view -S -b -h - ", true);
			}
			
			
			
			myWriter=(makeWriter ? new PrintWriter(myOutstream) : null);

			if(header!=null){
				if(myWriter!=null){
					myWriter.println(header);
				}else{
					byte[] temp=new byte[header.length()];
					for(int i=0; i<temp.length; i++){temp[i]=(byte)header.charAt(i);}
					try {
						myOutstream.write(temp);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}else if(OUTPUT_SAM){
				if(useSharedHeader){
					ArrayList<byte[]> list=SamReadInputStream.getSharedHeader(true);
					if(list==null){
						System.err.println("Header was null.");
					}else{
						try {
							for(byte[] line : list){
								myOutstream.write(line);
								myOutstream.write('\n');
								//myWriter.println(new String(line));
							}
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}else{
					if(myWriter!=null){
						myWriter.println(SamLine.header0());
						int a=(MINCHROM==-1 ? 1 : MINCHROM);
						int b=(MAXCHROM==-1 ? Data.numChroms : MAXCHROM);
						for(int chrom=a; chrom<=b; chrom++){
							//					myWriter.print(SamLine.header1(chrom, chrom));
							SamLine.printHeader1(chrom, chrom, myWriter);
						}
						myWriter.println(SamLine.header2());
					}else{
						ByteBuilder bb=new ByteBuilder(4096);
						SamLine.header0B(bb);
						bb.append('\n');
						int a=(MINCHROM==-1 ? 1 : MINCHROM);
						int b=(MAXCHROM==-1 ? Data.numChroms : MAXCHROM);
						for(int chrom=a; chrom<=b; chrom++){
							SamLine.printHeader1B(chrom, chrom, bb, myOutstream);
						}
						SamLine.header2B(bb);
						bb.append('\n');


						try {
							if(bb.length>0){myOutstream.write(bb.array, 0, bb.length);}
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
			}else if(ff.bread()){
				if(myWriter!=null){
					myWriter.println("#"+Read.header());
				}else{
					try {
						myOutstream.write(("#"+Read.header()).getBytes());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		assert(bufferSize>=1);
		queue=new ArrayBlockingQueue<Job>(bufferSize);
	}

	@Override
	public abstract void run();

	/** Uses this thread to transform reads to text, and the ReadStreamWriter thread to write text to disk */
	public final synchronized void addListAsText(ArrayList<Read> list){
		assert(false) : "TODO";
		addList(list, myWriter, myOutstream, false);
	}

	public final synchronized void poison(){
		addJob(new Job(null, null, null, false, true));
	}

	public final synchronized void addList(ArrayList<Read> list){
		addList(list, myWriter, myOutstream, false);
	}

	public final synchronized void addList(ArrayList<Read> l, PrintWriter w, OutputStream o, boolean c){
		boolean poison=(c && w!=null && w==myWriter);
		Job j=new Job(l, w, o, c, poison);
		addJob(j);
	}
	
	public final synchronized void addJob(Job j){
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
	
	protected static final StringBuilder toQualitySB(byte[] quals, int len){
		if(quals==null){return fakeQualitySB(30, len);}
		assert(quals.length==len);
		StringBuilder sb=new StringBuilder(NUMERIC_QUAL ? len*3+1 : len+1);
		if(NUMERIC_QUAL){
			if(quals.length>0){sb.append(quals[0]);}
			for(int i=1; i<quals.length; i++){
				sb.append(' ');
				sb.append(quals[i]);
			}
		}else{
			final byte b=FASTQ.ASCII_OFFSET_OUT;
			for(int i=0; i<quals.length; i++){
				sb.append((char)(b+quals[i]));
			}
		}
		return sb;
	}
	
	protected static final StringBuilder fakeQualitySB(int q, int len){
		StringBuilder sb=new StringBuilder(NUMERIC_QUAL ? len*3+1 : len+1);
		char c=(char)(q+FASTQ.ASCII_OFFSET_OUT);
		if(NUMERIC_QUAL){
			if(len>0){sb.append(q);}
			for(int i=1; i<len; i++){
				sb.append(' ');
				sb.append(q);
			}
		}else{
			for(int i=0; i<len; i++){sb.append(c);}
		}
		return sb;
	}
	
	protected static final ByteBuilder toQualityB(byte[] quals, int len, ByteBuilder bb){
		if(quals==null){return fakeQualityB(30, len, bb);}
		assert(quals.length==len);
		bb.ensureExtra(NUMERIC_QUAL ? len*3+1 : len+1);
		if(NUMERIC_QUAL){
			if(quals.length>0){bb.append((int)quals[0]);}
			for(int i=1; i<quals.length; i++){
				bb.append(' ');
				bb.append((int)quals[i]);
			}
		}else{
			final byte b=FASTQ.ASCII_OFFSET_OUT;
			for(int i=0; i<quals.length; i++){
				bb.append(b+quals[i]);
			}
		}
		return bb;
	}
	
	protected static final ByteBuilder fakeQualityB(int q, int len, ByteBuilder bb){
		bb.ensureExtra(NUMERIC_QUAL ? len*3+1 : len+1);
		if(NUMERIC_QUAL){
			int c=(q+FASTQ.ASCII_OFFSET_OUT);
			if(len>0){bb.append(q);}
			for(int i=1; i<len; i++){
				bb.append(' ');
				bb.append(q);
			}
		}else{
			byte c=(byte)(q+FASTQ.ASCII_OFFSET_OUT);
			for(int i=0; i<len; i++){bb.append(c);}
		}
		return bb;
	}

	/** Return true if this stream has detected an error */
	public final boolean errorState(){return errorState;}
	/** Return true if this stream has finished */
	public final boolean finishedSuccessfully(){return finishedSuccessfully;}
	/** TODO */
	protected boolean errorState=false;
	protected boolean finishedSuccessfully=false;
	
	public static int MINCHROM=-1; //For generating sam header
	public static int MAXCHROM=-1; //For generating sam header
	public static CharSequence HEADER;
	public static boolean NUMERIC_QUAL=true;
	public static boolean OUTPUT_SAM_SECONDARY_ALIGNMENTS=false;

	public final boolean OUTPUT_SAM;
	public final boolean OUTPUT_BAM;
	public final boolean OUTPUT_FASTQ;
	public final boolean OUTPUT_FASTA;
	public final boolean OUTPUT_ATTACHMENT;
	public final boolean OUTPUT_STANDARD_OUT;
	public final boolean SITES_ONLY;
	public boolean OUTPUT_INTERLEAVED=false;
	
	protected final boolean allowSubprocess;
	
	protected final boolean read1;
	protected final String fname;
	protected final String qfname;
	protected final OutputStream myOutstream;
	protected final PrintWriter myWriter;
	protected final OutputStream myQOutstream;
	protected final PrintWriter myQWriter;
	protected final ArrayBlockingQueue<Job> queue;
	
	protected long readsWritten=0;
	protected long basesWritten=0;
	protected long validReadsWritten=0;
	protected long validBasesWritten=0;
	public String fname(){return fname;}
	public long readsWritten(){return readsWritten;}
	public long basesWritten(){return basesWritten;}
	public long validReadsWritten(){return validReadsWritten;}
	public long validBasesWritten(){return validBasesWritten;}
	
	protected static class Job{
		public Job(ArrayList<Read> list_, PrintWriter writer_, OutputStream outstream_, boolean closeWhenDone_,
				boolean shutdownThread_){
			list=list_;
			writer=writer_;
			outstream=outstream_;
			close=closeWhenDone_;
			poison=shutdownThread_;
		}
		public Job(ArrayList<Read> list_, PrintWriter writer_){
			this(list_, writer_, null, false, false);
		}
		public boolean isEmpty(){return list==null || list.isEmpty();}
		public final ArrayList<Read> list;
		public final PrintWriter writer;
		public final OutputStream outstream;
		public final boolean close;
		public final boolean poison;
	}
	
}
