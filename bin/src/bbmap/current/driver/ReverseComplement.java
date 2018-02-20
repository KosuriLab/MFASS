package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;

import align2.ListNum;
import align2.Shared;
import align2.Tools;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Jul 19, 2013
 *
 */
public class ReverseComplement {
	
public static void main(String[] args){
		
		if(args==null || args.length==0 || (args.length==1 && 
				(args[0].equalsIgnoreCase("-h") || args[0].equals("-help") || args[0].equals("--help") || args[0].equals("-?") || args[0].equals("?")))){
			printOptions();
			System.exit(0);
		}
		ReverseComplement rc=new ReverseComplement(args);
		rc.process();
	}
	
	private static void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("\njava -ea -Xmx100m -cp <path> jgi.ReverseComplement <input file> <output file>");
		outstream.println("\nOptional flags:");
		outstream.println("in=<file>    \tThe 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.");
		outstream.println("out=<file>   \tThe 'out=' flag is needed if the output file is not the second parameter.  'out=stdout' will pipe to standard out.");
		outstream.println("showspeed=t  \tSet to 'f' to suppress display of processing speed.");
	}
	
	public ReverseComplement(String[] args){
		for(String s : args){if(s.contains("standardout") || s.contains("stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		boolean setOut=false;

		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("in")){
				in=b;
			}else if(a.equals("out")){
				out=b;
				setOut=true;
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("bf1")){
				ByteFile.FORCE_MODE_BF1=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF2=!ByteFile.FORCE_MODE_BF1;
			}else if(a.equals("bf2")){
				ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF1=!ByteFile.FORCE_MODE_BF2;
			}else if(a.equals("usegzip") || a.equals("gzip")){
				ReadWrite.USE_GZIP=Tools.parseBoolean(b);
			}else if(a.equals("usepigz") || a.equals("pigz")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					int zt=Integer.parseInt(b);
					if(zt<1){ReadWrite.USE_PIGZ=false;}
					else{
						ReadWrite.USE_PIGZ=true;
						if(zt>1){
							ReadWrite.MAX_ZIP_THREADS=zt;
							ReadWrite.ZIP_THREAD_DIVISOR=1;
						}
					}
				}else{ReadWrite.USE_PIGZ=Tools.parseBoolean(b);}
			}else if(a.equals("usegunzip") || a.equals("gunzip")){
				ReadWrite.USE_GUNZIP=Tools.parseBoolean(b);
			}else if(a.equals("useunpigz") || a.equals("unpigz")){
				ReadWrite.USE_UNPIGZ=Tools.parseBoolean(b);
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.equals("showspeed")){
				showspeed=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(i==0 && in==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				in=args[i];
			}else if(i==1 && out==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				out=args[i];
				setOut=true;
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(!setOut){out="stdout.fa";}
		if("stdout".equalsIgnoreCase(out) || "standarddout".equalsIgnoreCase(out)){
			out="stdout.fa";
			outstream=System.err;
		}
		if(!overwrite){
			if(out!=null && new File(out).exists()){throw new RuntimeException("Output file "+out+" already exists, and overwrite="+overwrite);}
		}
		assert(!in.equalsIgnoreCase(out));
	}
	
	public void process(){
		
		Timer t=new Timer();
		t.start();
		
		boolean dq0=FASTQ.DETECT_QUALITY;
		boolean ti0=FASTQ.TEST_INTERLEAVED;
		int rbl0=Shared.READ_BUFFER_LENGTH;
		FASTQ.DETECT_QUALITY=false;
		FASTQ.TEST_INTERLEAVED=false;
		Shared.READ_BUFFER_LENGTH=8;
		
		process2();
		
		FASTQ.DETECT_QUALITY=dq0;
		FASTQ.TEST_INTERLEAVED=ti0;
		Shared.READ_BUFFER_LENGTH=rbl0;
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                            \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException(this.getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	public void process2(){
		
		final TextStreamWriter tsw=(out==null ? null : new TextStreamWriter(out, overwrite, false, true));
		if(tsw!=null){tsw.start();}
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			FileFormat ff1=FileFormat.testInput(in, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, null);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		boolean paired=cris.paired();
		assert(paired);
		if(verbose){System.err.println("Paired: "+paired);}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		while(reads!=null && reads.size()>0){

			for(Read r : reads){
				assert(r.mate==null);
				readsProcessed++;
				basesProcessed+=r.bases==null ? 0 : r.bases.length;
				r.reverseComplement();
				if(tsw!=null){tsw.println(r);}
			}
			
			cris.returnList(ln, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln, ln.list.isEmpty());
		
		ReadWrite.closeStream(cris);
		
		if(tsw!=null){tsw.poisonAndWait();}
		
		errorState|=(cris.errorState() /*|| (tsw!=null && tsw.errorState())*/);
	}
	
	private String in, out;
	private long maxReads=-1;
	private long readsProcessed=0;
	private long basesProcessed=0;
	public boolean errorState=false;
	
	private static PrintStream outstream=System.err;
	public static boolean overwrite=false;
	public static boolean showspeed=true;
	public static boolean verbose=false;
	
}
