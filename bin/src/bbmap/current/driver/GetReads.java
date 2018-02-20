package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.RTextOutputStream3;
import stream.Read;

import align2.ListNum;
import align2.Tools;
import dna.Data;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextStreamWriter;

/**
 * Grab reads with specified numbers from a file.
 * TODO Note that much of this is ripped directly from ReformatReads, but is incorrect, because this class does not support dual output files.
 * @author Brian Bushnell
 * @date Jul 10, 2013
 *
 */
public class GetReads {
	
	public static void main(String[] args){
		new GetReads(args);
	}
	
	public GetReads(String[] args){
		if(args==null || args.length==0){
			throw new RuntimeException("No arguments.");
		}

		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		Timer t=new Timer();
		t.start();

		String in1=null;
		String in2=null;
		
		String qfin1=null;
		String qfin2=null;

		String out1=null;
		String out2=null;

		String qfout1=null;
		String qfout2=null;

		boolean parsecustom=false;
		boolean errorState=false;
		long maxReads=-1;
		int passes=1;
		boolean testsize=false;
		boolean overwrite=false;
		float samplerate=1f;
		long sampleseed=1;
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		byte qin=-1;
		byte qout=-1;

		FastaReadInputStream.SPLIT_READS=false;
		stream.FastaReadInputStream.MIN_READ_LEN=1;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		
		HashSet<Long> table=new HashSet<Long>();
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=(split.length>1 ? split[1] : "true");
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("null") || a.equals(in2)){
				// do nothing
			}else if(a.equals("id") || a.equals("number")){
				String[] b2=b.split(",");
				for(String c : b2){
					final long x, y; 
					if(c.indexOf('-')>=0){
						String[] c2=c.split("-");
						assert(c2.length==2) : c;
						x=Long.parseLong(c2[0]);
						y=Long.parseLong(c2[1]);
					}else{
						x=y=Long.parseLong(c);
					}
					for(long z=x; z<=y; z++){
						table.add(z);
					}
				}
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
//				align2.FastqReadInputStream.verbose=verbose;
			}else if(a.equals("bf1")){
				ByteFile.FORCE_MODE_BF1=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF2=!ByteFile.FORCE_MODE_BF1;
			}else if(a.equals("bf2")){
				ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF1=!ByteFile.FORCE_MODE_BF2;
			}else if(a.equals("usegzip") || a.equals("gzip")){
				ReadWrite.USE_GZIP=Tools.parseBoolean(b);
			}else if(a.equals("usepigz") || a.equals("pigz")){
				ReadWrite.USE_PIGZ=Tools.parseBoolean(b);
			}else if(a.equals("usegunzip") || a.equals("gunzip")){
				ReadWrite.USE_GUNZIP=Tools.parseBoolean(b);
			}else if(a.equals("useunpigz") || a.equals("unpigz")){
				ReadWrite.USE_UNPIGZ=Tools.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
				if(b.indexOf('#')>-1 && !new File(b).exists()){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				out1=b;
				if(b.indexOf('#')>-1){
					out1=b.replace("#", "1");
					out2=b.replace("#", "2");
				}
			}else if(a.equals("out2") || a.equals("output2")){
				out2=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfout1=b;
			}else if(a.equals("qfin2")){
				qfin2=b;
			}else if(a.equals("qfout2")){
				qfout2=b;
			}else if(a.equals("parsecustom")){
				parsecustom=Tools.parseBoolean(b);
			}else if(a.equals("testsize")){
				testsize=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.startsWith("fastareadlen")){
				FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
				FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
			}else if(a.startsWith("fastaminread") || a.startsWith("fastaminlen")){
				FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
			}else if(a.equals("ascii") || a.equals("quality") || a.equals("qual")){
				byte x;
				if(b.equalsIgnoreCase("sanger")){x=33;}
				else if(b.equalsIgnoreCase("illumina")){x=64;}
				else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;}
				else{x=(byte)Integer.parseInt(b);}
				qin=qout=x;
			}else if(a.equals("asciiin") || a.equals("qualityin") || a.equals("qualin") || a.equals("qin")){
				byte x;
				if(b.equalsIgnoreCase("sanger")){x=33;}
				else if(b.equalsIgnoreCase("illumina")){x=64;}
				else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=true;}
				else{x=(byte)Integer.parseInt(b);}
				qin=x;
			}else if(a.equals("asciiout") || a.equals("qualityout") || a.equals("qualout") || a.equals("qout")){
				byte x;
				if(b.equalsIgnoreCase("sanger")){x=33;}
				else if(b.equalsIgnoreCase("illumina")){x=64;}
				else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY_OUT=true;}
				else{x=(byte)Integer.parseInt(b);}
				qout=x;
			}else if(a.equals("qauto")){
				FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;
			}else if(a.equals("testinterleaved")){
				FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
				outstream.println("Set TEST_INTERLEAVED to "+FASTQ.TEST_INTERLEAVED);
				setInterleaved=true;
			}else if(a.equals("forceinterleaved")){
				FASTQ.FORCE_INTERLEAVED=Tools.parseBoolean(b);
				outstream.println("Set FORCE_INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				setInterleaved=true;
			}else if(a.equals("interleaved") || a.equals("int")){
				if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
				else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
					setInterleaved=true;
				}
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
				assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
			}else if(a.equals("sampleseed")){
				sampleseed=Long.parseLong(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				int x=Integer.parseInt(b);
				stream.FastaReadInputStream.MIN_READ_LEN=(x>0 ? x : Integer.MAX_VALUE);
			}else if(in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(out1==null && i==1 && !arg.contains("=")){
				out1=arg;
				if(arg.indexOf('#')>-1){
					out1=b.replace("#", "1");
					out2=b.replace("#", "2");
				}
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		assert(FastaReadInputStream.settingsOK());
//		if(maxReads!=-1){ReadWrite.USE_GUNZIP=ReadWrite.USE_UNPIGZ=false;}
		
		if(in1==null){
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(out1==null){
			if(out2!=null){
				throw new RuntimeException("Error - cannot define out2 without defining out1.");
			}
			out1="stdout";
		}
		
		if(!setInterleaved){
			assert(in1!=null && out1!=null) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, false, out1, out2)){
			throw new RuntimeException("\n\nOVERWRITE="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		FASTQ.PARSE_CUSTOM=parsecustom;
		
		
		if(qin!=-1 && qout!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY=false;
		}else if(qin!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.DETECT_QUALITY=false;
		}else if(qout!=-1){
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY_OUT=false;
		}
		

		FileFormat ffin=FileFormat.testInput(in1, 0, null, true, true);
		FileFormat ffout=FileFormat.testOutput(out1, 0, null, true, overwrite, false);
		
		
		final boolean useSharedHeader=(ffin!=null && ffout!=null && ffin.samOrBam() && ffout.samOrBam());
		
		if(ffin!=null && ffout!=null && ffin.samOrBam() && (ffout.samOrBam() || ffout.bread())){
			throw new RuntimeException("\nDirect conversion of sam to sam or bread are not currently supported.\nAll other conversions are possible.");
		}
		
		
		ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, useSharedHeader, ff1, ff2);
		}
		
		cris.setSampleRate(samplerate, sampleseed);
		outstream.println("Input is "+(cris.paired() ? "paired" : "unpaired"));
		Thread cristhread=new Thread(cris);
		cristhread.start();

		TextStreamWriter tsw=new TextStreamWriter(out1, overwrite, false, false);
		tsw.start();
		
		
		long readsProcessed=0;
		long basesProcessed=0;

		for(int pass=1; pass<=passes; pass++){
//			outstream.println("pass="+pass);
			if(pass>1){
				FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
				FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
				cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, useSharedHeader, ff1, ff2);
				cris.setSampleRate(samplerate, sampleseed);
				cristhread=new Thread(cris);
				cristhread.start();
			}
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(ffin.samOrBam() || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0 && !table.isEmpty()){

				for(Read r1 : reads){
					{
						readsProcessed++;
						basesProcessed+=r1.bases==null ? 0 : r1.bases.length;
					}
					Read r2=r1.mate;
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=r2.bases==null ? 0 : r2.bases.length;
					}
					
					if(table.remove(r1.numericID)){
						tsw.println(r1);
						if(r2!=null){tsw.println(r2);}
						if(table.isEmpty()){break;}
					}
				}

				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln, ln.list.isEmpty());
			errorState|=ReadWrite.closeStream(cris);
		}

		if(tsw!=null){
			tsw.poisonAndWait();
		}
		
		errorState|=(cris.errorState());
		
		t.stop();

		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}

		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		if(testsize){
			long bytesProcessed=(new File(in1).length()+(in2==null ? 0 : new File(in2).length()))*passes;
			double xpnano=bytesProcessed/(double)(t.elapsed);
			String xpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
			while(xpstring.length()<8){xpstring=" "+xpstring;}
			outstream.println("Bytes Processed:    "+xpstring+" \t"+String.format("%.2fm bytes/sec", xpnano*1000));
		}
		
		if(errorState){
			throw new RuntimeException("GetReads terminated in an error state; the output may be corrupt.");
		}
		
	}

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
