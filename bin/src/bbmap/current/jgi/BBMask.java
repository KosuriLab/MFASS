package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.LinkedHashMap;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.RTextOutputStream3;
import stream.Read;
import stream.SamLine;
import align2.ListNum;
import align2.Shared;
import align2.Tools;
import dna.AminoAcid;
import dna.Data;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;

/**
 * Masks a fasta file by inserting 'N' in place of low-complexity short repeats, 
 * and anything covered by mapped reads in a sam file.
 * 
 * @author Brian Bushnell
 * @date Feb 18, 2014
 *
 */
public class BBMask{

	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		BBMask masker=new BBMask(args);
		masker.process(t);
	}

	public BBMask(String[] args){
		
		if(args==null || args.length==0){
			printOptions();
			System.exit(0);
		}

		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");

		FastaReadInputStream.SPLIT_READS=false;
		stream.FastaReadInputStream.MIN_READ_LEN=1;
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.READ_BUFFER_NUM_BUFFERS=Tools.min(8, Shared.READ_BUFFER_NUM_BUFFERS);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=16;
		ReadWrite.ZIP_THREAD_DIVISOR=1;
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				//			align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
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
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("t") || a.equals("threads")){
				Shared.THREADS=Tools.max(Integer.parseInt(b), 1);
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("bf1")){
				ByteFile.FORCE_MODE_BF1=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF2=!ByteFile.FORCE_MODE_BF1;
			}else if(a.equals("bf2")){
				ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF1=!ByteFile.FORCE_MODE_BF2;
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1") || a.equals("ref")){
				inRef=b;
			}else if(a.equals("insam") || a.equals("samin") || a.equals("sam")){
				inSam=b.split(",");
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1") || a.equals("output1")){
				outRef=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfinRef=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfoutRef=b;
			}else if(a.equals("extin")){
				extinRef=b;
			}else if(a.equals("extout")){
				extoutRef=b;
			}else if(a.equals("mink") || a.equals("kmin")){	
				mink=Integer.parseInt(b);
			}else if(a.equals("maxk") || a.equals("kmax")){	
				maxk=Integer.parseInt(b);
			}else if(a.equals("k")){	
				mink=maxk=Integer.parseInt(b);
			}else if(a.equals("minlen")){	
				minlen=Integer.parseInt(b);
			}else if(a.equals("mincount")){	
				mincount=Integer.parseInt(b);
			}else if(a.equals("trd") || a.equals("trc") || a.equals("trimreaddescription")){
				Shared.TRIM_READ_COMMENTS=Tools.parseBoolean(b);
			}else if(a.equals("parsecustom")){
				parsecustom=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("fastareadlen") || a.equals("fastareadlength")){
				FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
				FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
			}else if(a.equals("fastaminread") || a.equals("fastaminlen") || a.equals("fastaminlength")){
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
			}else if(a.equals("tuc") || a.equals("touppercase")){
				Read.TO_UPPER_CASE=Tools.parseBoolean(b);
			}else if(a.equals("tossbrokenreads") || a.equals("tbr")){
				boolean x=Tools.parseBoolean(b);
				Read.NULLIFY_BROKEN_QUALITY=x;
				ConcurrentGenericReadInputStream.REMOVE_DISCARDED_READS=x;
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				int x=Integer.parseInt(b);
				stream.FastaReadInputStream.MIN_READ_LEN=(x>0 ? x : Integer.MAX_VALUE);
			}else if(inRef==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				inRef=arg;
			}else if(outRef==null && i==1 && !arg.contains("=")){
				outRef=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		assert(FastaReadInputStream.settingsOK());

		if(inRef==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.THREADS>2){
			//		if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}

		if(outRef==null){
			outRef="stdout";
		}

		if(outRef!=null && outRef.equalsIgnoreCase("null")){outRef=null;}

		if(!Tools.testOutputFiles(overwrite, false, outRef)){
			throw new RuntimeException("\n\nOVERWRITE="+overwrite+"; Can't write to output file "+outRef+"\n");
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

		ffoutRef=FileFormat.testOutput(outRef, FileFormat.FASTA, extoutRef, true, overwrite, false);

		ffinRef=FileFormat.testInput(inRef, FileFormat.FASTA, extinRef, true, true);
		
		if(inSam!=null && inSam.length>0){
			ffinSam=new FileFormat[inSam.length];
			for(int i=0; i<inSam.length; i++){
				ffinSam[i]=FileFormat.testInput(inSam[i], FileFormat.SAM, null, true, false);
			}
		}else{
			ffinSam=null;
		}
		
		SamLine.CONVERT_CIGAR_TO_MATCH=false;
	}


	public void process(Timer t0){
		
		Timer t=new Timer();
		{
			t.start();
			map=hashRef();
			t.stop();
//			outstream.println("Loaded "+repeats+" repeat bases in "+t);
		}
		
		long repeats=0, mapping=0;

		if(maxk>0){
			t.start();
			repeats=maskRepeats();
			t.stop();
			
			double rpnano=refReads/(double)(t.elapsed);
			double bpnano=refBases/(double)(t.elapsed);
			
			String rpstring=""+refReads;
			String bpstring=""+refBases;
			String bmstring=""+repeats;
	
			while(rpstring.length()<12){rpstring=" "+rpstring;}
			while(bpstring.length()<12){bpstring=" "+bpstring;}
			while(bmstring.length()<12){bmstring=" "+bmstring;}
	
			outstream.println("Time:                         \t"+t);
			outstream.println("Ref Scaffolds:          "+rpstring+" \t"+String.format("%.2fk scafs/sec", rpnano*1000000));
			outstream.println("Ref Bases:              "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
			outstream.println("Repeat Bases Masked:    "+bmstring);
		}
		
		if(ffinSam!=null){
			t.start();
			mapping=maskSam();
			t.stop();
			
			double rpnano=samReads/(double)(t.elapsed);
			double bpnano=samBases/(double)(t.elapsed);
			
			String rpstring=""+samReads;
			String bpstring=""+samBases;
			String bmstring=""+mapping;
	
			while(rpstring.length()<12){rpstring=" "+rpstring;}
			while(bpstring.length()<12){bpstring=" "+bpstring;}
			while(bmstring.length()<12){bmstring=" "+bmstring;}
	
			outstream.println("Time:                         \t"+t);
			outstream.println("Sam Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Sam Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
			outstream.println("Sam Bases Masked:       "+bmstring);
		}
		long total=repeats+mapping, masked=0;
		
		if(total>0 || true){
			masked=maskFromBitsets();
		}
		
		assert(total==masked) : repeats+", "+mapping+", "+total+", "+masked;
		
		{
			writeOutput();
			t0.stop();
			String tstring=""+total;
			while(tstring.length()<12){tstring=" "+tstring;}
			outstream.println("\nTotal Bases Masked:     "+tstring+"/"+refBases);
			outstream.println("Total Time:                   \t"+t);
//			outstream.println("Total Time:             "+t0);
		}
		
		
		
		if(errorState){
			throw new RuntimeException("\nBBMask terminated in an error state; the output may be corrupt.");
		}
	}
	
	private long maskFromBitsets(){
		long sum=0;
		for(Read r : map.values()){
			BitSet bs=((BitSet)r.obj);
			byte[] bases=r.bases;
			for(int i=0; i<bases.length; i++){
				if(bs.get(i)){
					if(bases[i]!='N'){sum++;}
					bases[i]='N';
				}else if(CONVERT_NON_ACGTN && !AminoAcid.isACGTN(bases[i])){
					bases[i]='N';
				}
			}
		}
		return sum;
	}
	
	private void writeOutput(){
		
		RTextOutputStream3 ros=null;
		if(ffoutRef!=null){
			final int buff=16;
			ros=new RTextOutputStream3(ffoutRef, null, qfoutRef, null, buff, null, false);
			ros.start();
		}
		
		long i=0;
		for(String name : map.keySet()){
			Read r=map.get(name);
			ArrayList<Read> list=new ArrayList<Read>(1);
			list.add(r);
			ros.add(list, i);
			i++;
		}
		errorState|=ReadWrite.closeStream(ros);
	}
	
	private long maskRepeats(){
		long sum=0;
		for(Read r : map.values()){
			sum+=maskRepeats(r, mink, maxk, mincount, minlen);
		}
		return sum;
	}
	
	private long maskSam(){
		long before=0, after=0;
		for(Read r : map.values()){
			before+=((BitSet)r.obj).cardinality();
		}
		for(FileFormat ff : ffinSam){
			maskSam(ff);
		}
		for(Read r : map.values()){
			after+=((BitSet)r.obj).cardinality();
		}
		return after-before;
	}
	
	private void maskSam(FileFormat ff){
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, colorspace, false, ff, null, null, null);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		
		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					assert(r.mate==null);

					final int initialLength1=(r.bases==null ? 0 : r.bases.length);

					{
						samReads++;
						samBases+=initialLength1;
						
						if(r.mapped()){
							SamLine sl=(SamLine)r.obj;
							assert(sl!=null) : "No sam line for read "+r;
							byte[] rname=sl.rname();
							assert(rname!=null) : "No rname for sam line "+sl;
							Read ref=map.get(new String(rname));
							assert(ref!=null) : "Could not find reference scaffold '"+rname+"' for samline \n"+sl+"\n in set \n"+map.keySet()+"\n";
							byte[] bases=ref.bases;
							BitSet bs=(BitSet)ref.obj;
							int start=Tools.max(0, sl.start());
							int stop=Tools.min(sl.stop()+1, ref.bases.length);
							bs.set(start, stop);
						}
						
					}
				}

				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln, ln.list==null || ln.list.isEmpty());
			}
		}

		errorState|=ReadWrite.closeStreams(cris);

		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}

	private LinkedHashMap<String, Read> hashRef(){


		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, colorspace, false, ffinRef, null, qfinRef, null);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		
		final LinkedHashMap<String, Read> hmr=new LinkedHashMap<String, Read>();
		
		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);

					final int initialLength1=(r.bases==null ? 0 : r.bases.length);

					r.obj=new BitSet(initialLength1);
					refReads++;
					refBases+=initialLength1;
					Read old=hmr.put(r.id, r);
					assert(old==null) : "Duplicate reference scaffold name "+r.id;
				}

				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln, ln.list==null || ln.list.isEmpty());
			}
		}

		errorState|=ReadWrite.closeStreams(cris);

		if(errorState){
			throw new RuntimeException("BBMask terminated in an error state; the output may be corrupt.");
		}
		
		return hmr;
	}
	
	/*--------------------------------------------------------------*/
	
	private static int maskRepeats(Read r, int mink, int maxk, int mincount, int minlen){
		
		final byte[] bases=r.bases;
		final BitSet bs;
		if(r.obj==null){
			r.obj=bs=new BitSet(bases.length);
		}else{
			bs=(BitSet)r.obj;
			for(int i=0; i<bases.length; i++){
				if(bases[i]=='N'){bs.set(i);}
			}
		}
		
		int before=bs.cardinality();
//		System.err.println("\nbefore="+before+"\n"+new String(bases)+"\n"+bs);
		
		for(int k=mink; k<=maxk; k++){
			maskRepeats(bases, bs, k, Tools.max(minlen, k*mincount));
		}

		int after=bs.cardinality();
		
//		System.err.println("before="+before+", after="+after+"\n"+new String(bases)+"\n"+bs);
		
		return after-before;
	}
	
	
	private static void maskRepeats(final byte[] bases, final BitSet bs, final int k, final int minlen){
		
		final int lim=bases.length-k;
		final int mask=(k>15 ? -1 : ~((-1)<<(2*k)));
		for(int loc=0; loc<lim; loc++){
			int len=repeatLength(bases, k, mask, loc);
			if(len>=minlen){
				int a=loc-k, b=loc-k+len;
				bs.set(a, b);
//				System.err.println("len="+len+", minlen="+minlen+", set "+(loc-k)+"-"+(loc-k+len));
				loc=Tools.max(loc, b-minlen);
//				System.err.println("Reset loc to "+loc);
			}else{
//				System.err.println("len="+len+" < minlen="+minlen);
			}
		}
		
	}
	
	
	private static int repeatLength(final byte[] bases, final int k, final int mask, final int loc){
		
		final int lim=bases.length;
		final int key=getInitialKey(bases, loc, k);
		if(key<0){return 0;}
		int kmer=key;
		int gap=0, last=-1;
		for(int i=loc; i<lim && gap<k; i++){
			final byte b=bases[i];
			final int n=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|n)&mask;
			if(kmer==key){
				last=i;
				gap=0;
			}else{
				gap++;
			}
//			System.err.println("i="+i+", lim="+lim+", gap="+gap+", last="+last+", b="+(char)b+", n="+n+", key="+key+", kmer="+kmer);
		}
		
//		System.err.println("k="+k+", mask="+mask+", loc="+loc+", last="+last);
		
		return (last<0 ? 0 : last-loc+k+1);
	}
	
	private static int getInitialKey(byte[] bases, int loc, int k){
		assert(k<16);
		int start=loc-k;
		int key=0;
		if(start<0){return -1;}
		for(int i=start; i<loc; i++){
			final byte b=bases[i];
			final int n=Dedupe.baseToNumber[b];
			key=(key<<2)|n;
		}
		assert(key>=0);
		return key;
	}

	/*--------------------------------------------------------------*/

	private void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("java -ea -Xmx15g -cp <path> jgi.BBMask ref=<file> sam=<file,file,...file> out=<file>");
		outstream.println("sam and out are optional.\n");
		outstream.println("Other parameters and their defaults:\n");
		outstream.println("overwrite=false  \tOverwrites files that already exist");
		outstream.println("ziplevel=2       \tSet compression level, 1 (low) to 9 (max)");
		outstream.println("fastawrap=100    \tLength of lines in fasta output");
		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
	}


	/*--------------------------------------------------------------*/

	private LinkedHashMap<String, Read> map=null;

	private long refReads=0;
	private long refBases=0;
//	private long repeatsMasked=0;
	
	private long samReads=0;
	private long samBases=0;
//	private long samMasked=0;
	
	public boolean errorState=false;

	private String inRef=null;
	private String inSam[]=null;

	private String qfinRef=null;

	private String outRef=null;

	private String qfoutRef=null;

	private String extinRef=null;
	private String extoutRef=null;

	private boolean parsecustom=false;
	private boolean overwrite=false;
	private boolean colorspace=false;

	private long maxReads=-1;

	private byte qin=-1;
	private byte qout=-1;

	private int mink=0;
	private int maxk=0;
	private int minlen=30;
	private int mincount=3;

	private final FileFormat ffinRef;
	private final FileFormat[] ffinSam;

	private final FileFormat ffoutRef;


	/*--------------------------------------------------------------*/

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public static boolean CONVERT_NON_ACGTN=true;

}
