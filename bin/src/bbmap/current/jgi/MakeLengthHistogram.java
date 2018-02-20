package jgi;

import java.io.File;
import java.util.ArrayList;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;

import dna.Data;
import dna.Timer;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;

import align2.ListNum;
import align2.Tools;

/**
 * @author Brian Bushnell
 * @date Jul 16, 2012
 *
 */
public class MakeLengthHistogram {
	
	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		
		String in1=null, in2=null;
		String out=null;
		
		Data.GENOME_BUILD=-1;
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("max") || a.equals("maxlength")){
				MAX_LENGTH=Integer.parseInt(b);
			}else if(a.startsWith("mult") || a.startsWith("div") || a.startsWith("bin")){
				MULT=Integer.parseInt(b);
			}else if(i==0 && !arg.contains("=")){
				in1=arg;
			}else if(i==1 && !arg.contains("=")){
				in2=arg;
			}else if(i==3 && !arg.contains("=")){
				out=arg;
			}else{
				throw new RuntimeException("Unknown argument: "+arg);
			}
		}
		
		MAX_LENGTH/=MULT;
		
		calc(in1, in2, out);
		t.stop();
		System.err.println("Time: \t"+t);
	}
	
	public static void calc(String in1, String in2, String out){
		if(fileIO.FileFormat.hasFastaExtension(in1)){
			FastaReadInputStream.SPLIT_READS=false;
			FastaReadInputStream.MIN_READ_LEN=1;
		}else{
			FASTQ.PARSE_CUSTOM=false;
		}
		long maxReads=-1;
		
		final ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, ff2);
//			if(verbose){System.err.println("Started cris");}
			Thread th=new Thread(cris);
			th.start();
		}
		boolean paired=cris.paired();
//		if(verbose){System.err.println("Paired: "+paired);}
		
		
		final int max=MAX_LENGTH;
		long[] hist=new long[max+1];
		long[] bhist=new long[max+1];
		
		int maxFound=0;
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(paired==(r.mate!=null));
			}
			
			while(reads!=null && reads.size()>0){
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){
					
//					System.out.println("Processing read "+r.numericID);
					
					if(r!=null && r.bases!=null){
						readsProcessed++;
						int x=r.bases.length;
						int y=Tools.min(max, (x+MULT/2)/MULT);
						hist[y]++;
						bhist[y]+=x;
						maxFound=Tools.max(maxFound, x);
					}
					
					if(r.mate!=null){
						readsProcessed++;
						Read r2=r.mate;
						int x=r2.bases.length;
						int y=Tools.min(max, (x+MULT/2)/MULT);
						hist[y]++;
						bhist[y]+=x;
						maxFound=Tools.max(maxFound, x);
					}
					
				}
				//System.err.println("returning list");
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
			System.err.println("Processed "+readsProcessed+" reads.");
		}

		double[] histF=new double[max+1];
		long[] histC=new long[max+1];
		double[] histCF=new double[max+1];
		
		double[] bhistF=new double[max+1];
		long[] bhistC=new long[max+1];
		double[] bhistCF=new double[max+1];
		

		histC[max]=hist[max];
		bhistC[max]=bhist[max];
		for(int i=max; i>0; i--){
			histC[i-1]=histC[i]+hist[i-1];
			bhistC[i-1]=bhistC[i]+bhist[i-1];
		}
		for(int i=0; i<=max; i++){
			histCF[i]=histC[i]*100d/histC[0];
			histF[i]=hist[i]*100d/histC[0];
			bhistCF[i]=bhistC[i]*100d/bhistC[0];
			bhistF[i]=bhist[i]*100d/bhistC[0];
		}

		TextStreamWriter tsw=new TextStreamWriter(out==null ? "stdout" : out, true, false, false);
		tsw.start();
		tsw.println("Reads:      \t"+readsProcessed);
		tsw.println("Bases:      \t"+bhistC[0]);
		tsw.println("Avg Length: \t"+String.format("%.1f",(bhistC[0]*1d/readsProcessed)));
		tsw.println("Read Length Histogram:\n");
		tsw.println("Length\treads\t%reads\tcum reads\tcum %reads\tbases\t%bases\tcum bases\tcum %bases");
		for(int i=0; i<=max; i++){
			tsw.println((i*MULT)+"\t"+hist[i]+String.format("\t%.3f%%", histF[i])+"\t"+histC[i]+String.format("\t%.3f%%", histCF[i])+
					"\t"+bhist[i]+String.format("\t%.3f%%", bhistF[i])+"\t"+bhistC[i]+String.format("\t%.3f%%", bhistCF[i]));
			if(i*MULT>=maxFound){break;}
		}
		tsw.poisonAndWait();
	}
	
	public static long readsProcessed=0;
	public static int MAX_LENGTH=4000;
	public static int MULT=10;
	
}
