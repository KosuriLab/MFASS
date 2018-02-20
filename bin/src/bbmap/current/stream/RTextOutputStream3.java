package stream;

import java.io.File;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.FileFormat;

public class RTextOutputStream3 {
	
	public RTextOutputStream3(String fname, String mate_fname, int maxSize, 
			boolean ordered, boolean sam, boolean bam, boolean fastq, boolean fasta, boolean attachment, boolean overwrite, boolean sitesonly){
		this(fname, mate_fname, null, null, maxSize, ordered, sam, bam, fastq, fasta, attachment, overwrite, sitesonly, false);
	}
	
	public RTextOutputStream3(String fname, String mate_fname, String qfname, String mate_qfname, int maxSize, 
			boolean ordered, boolean sam, boolean bam, boolean fastq, boolean fasta, boolean attachment, boolean overwrite, boolean sitesonly, boolean useSharedHeader){
//		System.err.println("Called RTextOutputStream3 with fname="+fname+", mate_fname="+mate_fname+", qfname="+qfname+", mate_qfname="+mate_qfname);
		STANDARD_OUT=("standardout".equalsIgnoreCase(fname) || fname.toLowerCase().startsWith("standardout.")
				|| "stdout".equalsIgnoreCase(fname) || fname.toLowerCase().startsWith("stdout."));
		
		if(verbose){
			System.err.println("RTextOutputStream3("+fname+", "+mate_fname+", "+qfname+", "+mate_qfname+", "+maxSize+", "+ordered+")");
		}
		
		ORDERED=ordered;
		ATTACHMENT=attachment;
		SAM=!ATTACHMENT && (sam || bam);
		BAM=!ATTACHMENT && bam;
		FASTQ=!ATTACHMENT && fastq;
		FASTA=!ATTACHMENT && fasta;
		SITESONLY=!ATTACHMENT && sitesonly;
		
		assert(((SAM ? 1 : 0)+(FASTQ ? 1 : 0)+(FASTA ? 1 : 0)+(ATTACHMENT ? 1 : 0)+(SITESONLY ? 1 : 0))<=1) : 
			SAM+", "+SITESONLY+", "+FASTQ+", "+FASTA+", "+ATTACHMENT;
		
		if(fname!=null && !fname.equals("/dev/null")){
			File f=new File(fname);
			assert(overwrite || !f.exists()) : f.getAbsolutePath()+" already exists; please delete it.";
			if(mate_fname!=null){assert(!fname.equals(mate_fname)) : fname+"=="+mate_fname;}
		}
		
		boolean allowSubprocess=true; //TODO - make parameter
		
		if(BYTE_WRITER){
			readstream1=new ReadStreamByteWriter(fname, qfname, true, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, STANDARD_OUT, useSharedHeader, allowSubprocess);
			readstream2=STANDARD_OUT ? null : ((SAM || SITESONLY || mate_fname==null) ? null : 
				new ReadStreamByteWriter(mate_fname, mate_qfname, false, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, false, useSharedHeader, allowSubprocess));
		}else{
			readstream1=new ReadStreamStringWriter(fname, qfname, true, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, STANDARD_OUT, useSharedHeader, allowSubprocess);
			readstream2=STANDARD_OUT ? null : ((SAM || SITESONLY || mate_fname==null) ? null : 
				new ReadStreamStringWriter(mate_fname, mate_qfname, false, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, false, useSharedHeader, allowSubprocess));
		}
		
		if(readstream2==null && readstream1!=null){
//			System.out.println("RTextOutputStream3 detected interleaved output.");
			readstream1.OUTPUT_INTERLEAVED=true;
		}
		
		table=(ORDERED ? new HashMap<Long, ArrayList<Read>>(MAX_CAPACITY) : null);
		
		assert(readstream1==null || readstream1.read1==true);
		assert(readstream2==null || (readstream2.read1==false));
//		assert(false) : ATTACHMENT;
	}
	
	public RTextOutputStream3(FileFormat ff1, FileFormat ff2, int maxSize, CharSequence header, boolean useSharedHeader){
		this(ff1, ff2, null, null, maxSize, header, useSharedHeader);
	}
	
	public RTextOutputStream3(FileFormat ff1, FileFormat ff2, String qf1, String qf2, int maxSize, CharSequence header, boolean useSharedHeader){
		
		if(verbose){
			System.err.println("RTextOutputStream3("+ff1+", "+ff2+", "+qf1+", "+qf2+", "+maxSize+", "+useSharedHeader+")");
		}

		assert(ff1!=null);
		assert(!ff1.text() && !ff1.unknownFormat()) : "Unknown format for "+ff1;
		
		FASTA=ff1.fasta();
//		boolean bread=(ext==TestFormat.txt);
		SAM=(ff1.sam() || ff1.bam());
		BAM=ff1.bam();
		ATTACHMENT=ff1.attachment();
		SITESONLY=ff1.sites();
		FASTQ=ff1.fastq();
		STANDARD_OUT=ff1.stdio();
		
		assert(((SAM ? 1 : 0)+(FASTQ ? 1 : 0)+(FASTA ? 1 : 0)+(ATTACHMENT ? 1 : 0)+(SITESONLY ? 1 : 0))<=1) : 
			SAM+", "+SITESONLY+", "+FASTQ+", "+FASTA+", "+ATTACHMENT;
		
		ORDERED=ff1.ordered();
		if(ff1.hasName() && ff1.devnull()){
			File f=new File(ff1.name());
			assert(ff1.overwrite() || !f.exists()) : f.getAbsolutePath()+" already exists; please delete it.";
			if(ff2!=null){assert(!ff1.name().equals(ff2.name())) : ff1.name()+"=="+ff2.name();}
		}
		
		if(BYTE_WRITER){
			readstream1=new ReadStreamByteWriter(ff1, qf1, true, maxSize, header, useSharedHeader);
			readstream2=ff1.stdio() || ff2==null ? null : new ReadStreamByteWriter(ff2, qf2, false, maxSize, header, useSharedHeader);
			
//			readstream1=new ReadStreamByteWriter(ff1, qfname, true, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, STANDARD_OUT, useSharedHeader, ff.allowSubprocess());
//			readstream2=STANDARD_OUT ? null : ((SAM || SITESONLY || mate_fname==null) ? null : 
//				new ReadStreamByteWriter(mate_fname, mate_qfname, false, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, false, useSharedHeader, ff.allowSubprocess()));
		}else{
			readstream1=new ReadStreamStringWriter(ff1, qf1, true, maxSize, header, useSharedHeader);
			readstream2=ff1.stdio() || ff2==null ? null : new ReadStreamStringWriter(ff2, qf2, false, maxSize, header, useSharedHeader);
			
//			readstream1=new ReadStreamStringWriter(fname, qfname, true, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, STANDARD_OUT, useSharedHeader, ff.allowSubprocess());
//			readstream2=STANDARD_OUT ? null : ((SAM || SITESONLY || mate_fname==null) ? null : 
//				new ReadStreamStringWriter(mate_fname, mate_qfname, false, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, false, useSharedHeader, ff.allowSubprocess()));
		}
		
		if(readstream2==null && readstream1!=null){
//			System.out.println("RTextOutputStream3 detected interleaved output.");
			readstream1.OUTPUT_INTERLEAVED=true;
		}
		
		table=(ORDERED ? new HashMap<Long, ArrayList<Read>>(MAX_CAPACITY) : null);
		
		assert(readstream1==null || readstream1.read1==true);
		assert(readstream2==null || (readstream2.read1==false));
	}
	
//	public RTextOutputStream3(String fname, String mate_fname, int maxSize, 
//			boolean ordered, int[] exts, CharSequence header, boolean overwrite, boolean useSharedHeader){
//		this(fname, mate_fname, null, null, maxSize, ordered, exts, header, overwrite, useSharedHeader);
//	}
//	
//	public RTextOutputStream3(String fname, String mate_fname, String qfname, String mate_qfname, int maxSize, 
//			boolean ordered, int[] exts, CharSequence header, boolean overwrite, boolean useSharedHeader){
//		
//		if(verbose){
//			System.err.println("RTextOutputStream3("+fname+", "+mate_fname+", "+qfname+", "+mate_qfname+", "+maxSize+", "+ordered+")");
//		}
//		
//		if(exts==null){exts=FileFormat.testFormat(fname, false);}
//		int ext=exts[0], zip=exts[1], itype=exts[2];
//		FASTA=ff.fasta();
////		boolean bread=(ext==TestFormat.txt);
//		SAM=ff.samOrBam();
//		BAM=ff.bam();
//		ATTACHMENT=ff.attachment();
//		SITESONLY=ff.sites();
//		FASTQ=(ff.fastq() || (ext==FileFormat.UNKNOWN && !SAM && !ATTACHMENT && !SITESONLY));
//		STANDARD_OUT=ff.stdio();
//		
//		assert(((SAM ? 1 : 0)+(FASTQ ? 1 : 0)+(FASTA ? 1 : 0)+(ATTACHMENT ? 1 : 0)+(SITESONLY ? 1 : 0))<=1) : 
//			SAM+", "+SITESONLY+", "+FASTQ+", "+FASTA+", "+ATTACHMENT;
//		
//		ORDERED=ordered;
//		if(fname!=null && !fname.equals("/dev/null")){
//			File f=new File(fname);
//			assert(overwrite || !f.exists()) : f.getAbsolutePath()+" already exists; please delete it.";
//			if(mate_fname!=null){assert(!fname.equals(mate_fname)) : fname+"=="+mate_fname;}
//		}
//		
//		boolean allowSubprocess=true; //TODO - make parameter
//		
//		if(BYTE_WRITER){
//			readstream1=new ReadStreamByteWriter(fname, qfname, true, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, STANDARD_OUT, useSharedHeader, allowSubprocess);
//			readstream2=STANDARD_OUT ? null : ((SAM || SITESONLY || mate_fname==null) ? null : 
//				new ReadStreamByteWriter(mate_fname, mate_qfname, false, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, false, useSharedHeader, allowSubprocess));
//		}else{
//			readstream1=new ReadStreamStringWriter(fname, qfname, true, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, STANDARD_OUT, useSharedHeader, allowSubprocess);
//			readstream2=STANDARD_OUT ? null : ((SAM || SITESONLY || mate_fname==null) ? null : 
//				new ReadStreamStringWriter(mate_fname, mate_qfname, false, maxSize, SAM, BAM, FASTQ, FASTA, SITESONLY, ATTACHMENT, false, useSharedHeader, allowSubprocess));
//		}
//		
//		if(readstream2==null && readstream1!=null){
////			System.out.println("RTextOutputStream3 detected interleaved output.");
//			readstream1.OUTPUT_INTERLEAVED=true;
//		}
//		
//		table=(ORDERED ? new HashMap<Long, ArrayList<Read>>(MAX_CAPACITY) : null);
//		
//		assert(readstream1==null || readstream1.read1==true);
//		assert(readstream2==null || (readstream2.read1==false));
//	}
	
	public synchronized void add(ArrayList<Read> list, long listnum){
		
		if(ORDERED){
			int size=table.size();
//			System.err.print(size+", ");
			final boolean flag=(size>=HALF_LIMIT);
			if(listnum>nextListID && size>=ADD_LIMIT){
				System.err.println("Output buffer became full; key "+listnum+" waiting on "+nextListID+".");
				while(listnum>nextListID && size>=HALF_LIMIT){
					try {
						this.wait(20000);
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
					size=table.size();
				}
				System.err.println("Output buffer became clear for key "+listnum+"; next="+nextListID+", size="+size);
			}
			addOrdered(list, listnum);
			assert(listnum!=nextListID);
			if(flag && listnum<nextListID){this.notifyAll();}
		}else{
			addDisordered(list, listnum);
		}
	}
	
	private synchronized void addOrdered(ArrayList<Read> list, long listnum){
//		System.err.println("RTOS got "+listnum+" of size "+(list==null ? "null" : list.size())+
//				" with first read id "+(list==null || list.isEmpty() || list.get(0)==null ? "null" : ""+list.get(0).numericID));
		assert(list!=null) : listnum;
		assert(listnum>=nextListID) : listnum+", "+nextListID;
//		assert(list.isEmpty() || list.get(0)==null || list.get(0).numericID>=nextReadID) : list.get(0).numericID+", "+nextReadID;
		assert(!table.containsKey(listnum));
		
		table.put(listnum, new ArrayList<Read>(list));
		
		while(table.containsKey(nextListID)){
//			System.err.println("Writing list "+first.get(0).numericID);
			ArrayList<Read> value=table.remove(nextListID);
			write(value);
			nextListID++;
		}
		if(table.isEmpty()){notifyAll();}
	}
	
	private synchronized void addDisordered(ArrayList<Read> list, long listnum){
		assert(list!=null);
		assert(table==null);
		write(new ArrayList<Read>(list));
	}
	
	private synchronized void write(ArrayList<Read> list){
		if(readstream1!=null){
			if(readstream1.getState()==State.TERMINATED){throw new RuntimeException("Writing to a terminated thread.");}
			readstream1.addList(list);
		}
		if(readstream2!=null){
			if(readstream1.getState()==State.TERMINATED){throw new RuntimeException("Writing to a terminated thread.");}
			readstream2.addList(list);
		}
	}
	
	public synchronized void close(){
		
		assert(table==null || table.isEmpty()); //Seems like a race condition.  Probably, I should wait at this point until the condition is true before proceeding.
		
//		readstream1.addList(null);
//		if(readstream2!=null){readstream2.addList(null);}
		readstream1.poison();
		if(readstream2!=null){readstream2.poison();}
	}
	
	public synchronized void start(){
		if(started){
			System.err.println("Resetting output stream.");
			nextListID=0;
			throw new RuntimeException();
		}else{
			started=true;
			if(readstream1!=null){readstream1.start();}
			if(readstream2!=null){readstream2.start();}
		}
	}
	
	public void join(){
		while(readstream1!=null && readstream1.getState()!=Thread.State.TERMINATED){
			try {
				readstream1.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		while(readstream2!=null && readstream2.getState()!=Thread.State.TERMINATED){
			try {
				if(readstream2!=null){readstream2.join();}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		assert(table==null || table.isEmpty());
		finishedSuccessfully=true;
	}
	
	public synchronized void resetNextListID(){
		for(int i=0; i<2000 && !table.isEmpty(); i++){
			try {this.wait(2000);} 
			catch (InterruptedException e) {e.printStackTrace();}
		}
		if(!table.isEmpty()){
			System.err.println("WARNING! resetNextListID() waited a long time and the table never cleared.  Process may have stalled.");
		}
		while(!table.isEmpty()){
			try {this.wait(2000);} 
			catch (InterruptedException e) {e.printStackTrace();}
		}
		nextListID=0;
	}
	
	public final String fname(){
//		if(STANDARD_OUT){return "stdout";}
		return readstream1.fname();
	}
	
	/** Return true if this stream has detected an error */
	public boolean errorState(){
		return errorState || (readstream1!=null && readstream1.errorState()) || (readstream2!=null && readstream2.errorState());
	}

	public boolean finishedSuccessfully(){
		return finishedSuccessfully && (readstream1==null || readstream1.finishedSuccessfully()) && (readstream2==null || readstream2.finishedSuccessfully());
	}

	private boolean errorState=false;
	private boolean finishedSuccessfully=false;

	public final boolean ORDERED;
	public final boolean SAM;
	public final boolean BAM;
	public final boolean FASTQ;
	public final boolean FASTA;
	public final boolean ATTACHMENT;
	public final boolean SITESONLY;

	public final boolean STANDARD_OUT;

	public final ReadStreamWriter getRS1(){return readstream1;}
	public final ReadStreamWriter getRS2(){return readstream2;}
	
	private final ReadStreamWriter readstream1;
	private final ReadStreamWriter readstream2;
	private long nextListID=0;
	private boolean started=false;
	
	/** Number of lists held before the stream blocks */
	private final int MAX_CAPACITY=256;
	private final int ADD_LIMIT=MAX_CAPACITY-2;
	private final int HALF_LIMIT=ADD_LIMIT/2;
	
	/** For ordered output */
	private final HashMap<Long, ArrayList<Read>> table;
	
	{if(HALF_LIMIT<1){throw new RuntimeException("Capacity too low.");}}
	
	public static boolean BYTE_WRITER=true;
	public static boolean verbose=false;
	
}
