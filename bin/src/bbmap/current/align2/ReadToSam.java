package align2;

import java.io.File;
import java.util.ArrayList;

import stream.ConcurrentReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.ReadStreamStringWriter;
import stream.ReadStreamWriter;

import dna.Data;
import dna.Timer;
import fileIO.ReadWrite;

public class ReadToSam {
	
	
	public static void main(String[] args){
		
		for(String arg : args){
			if(arg.startsWith("b=") || arg.startsWith("build=")){
				String[] split=arg.split("=");
				Data.setGenome(Integer.parseInt(split[1]));
			}
		}
		
		System.err.println("Using header for build "+Data.GENOME_BUILD);
		String reads1=args[0];
		String reads2=args[1].equalsIgnoreCase("null") ?  null : args[1];
		String outname=args[2].equalsIgnoreCase("null") ?  "" : args[2];
		
		ReadToSam smr=new ReadToSam(reads1, reads2, outname);
		smr.process();
		
	}
	
	public ReadToSam(String fname1, String fname2, String outname_){
		this(new RTextInputStream(fname1, fname2, -1), outname_);
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
	}
	
	public ReadToSam(RTextInputStream stream_, String outname_){
		stream=stream_;
		outname=outname_;
		paired=stream.paired();
//		assert(outname.contains("#")) : "Output file name must contain the character '#' to be used for chromosome number.";
		
		cris=(USE_CRIS ? new ConcurrentReadInputStream(stream, -1) : null);
	}
	
	public void process(){

		Timer t=new Timer();
		t.start();
		
		final String fname1=outname.replaceFirst("#", "1");
		if(fname1!=null && new File(fname1).exists()){throw new RuntimeException("Destination file "+fname1+" already exists.");}
		
		ReadStreamWriter wt1=new ReadStreamStringWriter(fname1, true, 4, true, fname1.endsWith(".bam"), false, false, false, false, false, false, true);

		Thread wtt1=new Thread(wt1);
		
//		while(t.hashCode()!=0){}

		if(wtt1!=null){wtt1.start();}
		
		assert(USE_CRIS);
		new Thread(cris).start();
		System.err.println("Started cris");
		
//		System.err.println
		
		long count=0;
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(reads!=null && reads.size()>0){
				
				ArrayList<Read> reads2=new ArrayList<Read>(reads.size());
				
				if(paired){
					for(Read r : reads){
//						r.setPaired(false);
//						if(r.mapped() && !r.discarded() && r.valid()){
//							reads2.add(r);
//						}
						assert(r!=null);
						assert(r.mate!=null);
						reads2.add(r);
						reads2.add(r.mate);
					}
				}else{
					for(Read r : reads){
						r.setPaired(false);
						if(r.mapped() && !r.discarded() && r.valid()){
							reads2.add(r);
						}
					}
				}
				
//				ArrayList<Read> reads2=(ArrayList<Read>) reads.clone();
				
				if(wt1!=null){wt1.addList(reads2);}
				
				//					System.err.println("Added list of length "+reads.size());

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
		}

		//Add poison
		//			if(wt1!=null){wt1.addList(null);}
		//			if(wt2!=null){wt2.addList(null);}
		wt1.poison();

		while(wtt1.isAlive()){
			try {
				wtt1.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		t.stop();
		Data.sysout.println("Time:\t"+t);
		
//		if(cris!=null){
//			new Thread(cris).start();
//			ListNum<Read> ln=cris.nextList();
//			ArrayList<Read> reads=(ln!=null ? ln.list : null);
//			
//			while(reads!=null && reads.size()>0){
//				processReads(reads);
//				cris.returnList(ln, ln.list.isEmpty());
//				ln=cris.nextList();
//				reads=(ln!=null ? ln.list : null);
//			}
//			cris.returnList(ln, ln.list.isEmpty());
//		}else{
//			ArrayList<Read> reads=stream.nextList();
//			while(reads!=null && reads.size()>0){
//				processReads(reads);
//				reads=stream.nextList();
//			}
//		}
//		
//		synchronized(this){this.notifyAll();}
//		
//		finish();
	}
	
	
	public final String outname;
	private final RTextInputStream stream;
	private final ConcurrentReadInputStream cris;
	
	public final boolean paired;
	
	public static boolean USE_CRIS=true; //Similar speed either way.  "true" may be better with many threads.
	
	public static final int WRITE_BUFFER=400; //Bigger number uses more memory, for less frequent writes.
	
	
}
