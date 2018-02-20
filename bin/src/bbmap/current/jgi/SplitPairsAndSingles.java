package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.RTextOutputStream3;
import stream.Read;

import align2.ListNum;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.FileFormat;

/**
 * @author Brian Bushnell
 * @date Sep 4, 2013
 *
 */
public final class SplitPairsAndSingles {
	

	
	public static void main(String[] args){
		
		if(args==null || args.length==0 || (args.length==1 && 
				(args[0].equalsIgnoreCase("-h") || args[0].equals("-help") || args[0].equals("--help") || args[0].equals("-?") || args[0].equals("?")))){
			printOptions();
			System.exit(0);
		}
		SplitPairsAndSingles dd=new SplitPairsAndSingles(args);
		dd.process();
	}
	
	private static void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("\njava -ea -Xmx100m -cp <path> jgi.SplitPairsAndSingles in=<input file> out=<pair output file> outs=<singleton output file> minlen=20");
		outstream.println("\nOptional flags:");
		outstream.println("in=<file>          \tThe 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.");
		outstream.println("in2=<file>         \tUse this if 2nd read of pairs are in a different file.");
		outstream.println("out=<file>         \tThe 'out=' flag is needed if the output file is not the second parameter.  'out=stdout' will pipe to standard out.");
		outstream.println("out2=<file>        \tUse this to write 2nd read of pairs to a different file.");
		outstream.println("outsingle=<file>   \t(outs) Write singleton reads here.");
		outstream.println("");
		outstream.println("overwrite=t        \t(ow) Set to false to force the program to abort rather than overwrite an existing file.");
		outstream.println("showspeed=t        \t(ss) Set to 'f' to suppress display of processing speed.");
		outstream.println("interleaved=auto   \t(int) If true, forces fastq input to be paired and interleaved.");
		outstream.println("qtrim=f            \tTrim read ends to remove bases with quality below minq.");
		outstream.println("                   \tValues: t (trim both ends), f (neither end), r (right end only), l (left end only).");
		outstream.println("trimq=4            \tTrim quality threshold.");
		outstream.println("minlen=2           \t(ml) Reads shorter than this after trimming will be discarded.");
		outstream.println("ziplevel=20        \t(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.");
		outstream.println("fixpairs=f         \t(fp, fint) Fixes corrupted interleaved files by examining paired read names.");
		
	}
	
	public SplitPairsAndSingles(String[] args){
		for(String s : args){if(s.contains("standardout") || s.contains("stdout")){outstream=System.err;}}
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		FastaReadInputStream.SPLIT_READS=false;
		ByteFile.FORCE_MODE_BF2=Shared.THREADS>2;
		boolean setOut=false, setOuts=false, trimRight_=false, trimLeft_=false, setInterleaved=false, fixPairs_=false;
		
		{
			boolean b=false;
			assert(b=true);
			EA=b;
		}

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
			}else if(a.equals("out") || a.equals("out1") || a.equals("outp") || a.equals("outp1") || a.equals("outpair") || a.equals("outpair1")){
				out1=b;
				setOut=true;
			}else if(a.equals("out2") || a.equals("outp2") || a.equals("outpair2")){
				out2=b;
			}else if(a.equals("outs") || a.equals("outsingle") || a.equals("outb") || a.equals("outbad")){
				outsingle=b;
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
			}else if(a.equals("interleaved") || a.equals("int")){
				if("auto".equalsIgnoreCase(b)){
					FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);
				}else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
				setInterleaved=true;
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
			}else if(a.equals("trim") || a.equals("qtrim")){
				if(b==null){trimRight_=trimLeft_=true;}
				else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){trimLeft_=true;trimRight_=false;}
				else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){trimLeft_=false;trimRight_=true;}
				else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){trimLeft_=trimRight_=true;}
				else{trimRight_=trimLeft_=Tools.parseBoolean(b);}
			}else if(a.equals("optitrim") || a.equals("otf") || a.equals("otm")){
				if(b!=null && (b.charAt(0)=='.' || Character.isDigit(b.charAt(0)))){
					TrimRead.optimalMode=true;
					TrimRead.optimalBias=Float.parseFloat(b);
					assert(TrimRead.optimalBias>=0 && TrimRead.optimalBias<1);
				}else{
					TrimRead.optimalMode=Tools.parseBoolean(b);
				}
			}else if(a.equals("trimright")){
				trimRight_=Tools.parseBoolean(b);
			}else if(a.equals("trimleft")){
				trimLeft_=Tools.parseBoolean(b);
			}else if(a.equals("trimq") || a.equals("trimquality") || a.equals("minq")){
				trimq=Byte.parseByte(b);
			}else if(a.equals("fixpairs") || a.equals("fp") || a.equals("fint")){
				fixPairs_=Tools.parseBoolean(b);
			}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength") || a.equals("minreadlength")){
				minReadLength=Integer.parseInt(b);
				assert(minReadLength>=0) : "minReadLength must be at least 0";
			}else if(a.equals("asciiin") || a.equals("qualityin") || a.equals("qualin") || a.equals("qin")){
				if(b.equalsIgnoreCase("auto")){
					FASTQ.DETECT_QUALITY=true;
				}else{
					byte ascii_offset=Byte.parseByte(b);
					FASTQ.ASCII_OFFSET=ascii_offset;
					System.err.println("Set fastq input ASCII offset to "+FASTQ.ASCII_OFFSET);
					FASTQ.DETECT_QUALITY=false;
				}
			}else if(a.equals("asciiout") || a.equals("qualityout") || a.equals("qualout") || a.equals("qout")){
				if(b.equalsIgnoreCase("auto")){
					FASTQ.DETECT_QUALITY_OUT=true;
				}else{
					byte ascii_offset=Byte.parseByte(b);
					FASTQ.ASCII_OFFSET_OUT=ascii_offset;
					System.err.println("Set fastq output ASCII offset to "+FASTQ.ASCII_OFFSET_OUT);
					FASTQ.DETECT_QUALITY_OUT=false;
				}
			}else if(a.equals("qauto")){
				FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;
			}else if(i==0 && in1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				in1=args[i];
			}else if(i==1 && out1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				out1=args[i];
				setOut=true;
			}else if(i==2 && outsingle==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				outsingle=args[i];
				setOuts=true;
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		fixPairs=fixPairs_;
		trimRight=trimRight_;
		trimLeft=trimLeft_;
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(in1!=null && in1.contains("#") && !new File(in1).exists()){
			int pound=in1.lastIndexOf('#');
			String a=in1.substring(0, pound);
			String b=in1.substring(pound+1);
			in1=a+1+b;
			in2=a+2+b;
		}
		if(in2!=null && (in2.contains("=") || in2.equalsIgnoreCase("null"))){in2=null;}
		
		if(fixPairs){
			if(in2!=null){
				System.err.println("ERROR: 'FixPairs' mode only works with a single interleaved input file, not paired input files.");
				System.err.println("Aborting.");
				System.exit(1);
			}
			setInterleaved=true;
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
			outstream.println("Paired input disabled; running in FixPairs mode");
		}
		
		if(!setInterleaved && in2==null){
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true;
			outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(out1!=null && out1.contains("#")){
			int pound=out1.lastIndexOf('#');
			String a=out1.substring(0, pound);
			String b=out1.substring(pound+1);
			out1=a+1+b;
			out2=a+2+b;
		}

		if(!setOut){
			out1="stdout.fq";
			outstream=System.err;
			out2=null;
		}else if("stdout".equalsIgnoreCase(out1) || "standarddout".equalsIgnoreCase(out1)){
			out1="stdout.fq";
			outstream=System.err;
			out2=null;
		}
		if(out1!=null && !Tools.canWrite(out1, overwrite)){throw new RuntimeException("Output file "+out1+" already exists, and overwrite="+overwrite);}

		assert(!in1.equalsIgnoreCase(out1));
		assert(!in1.equalsIgnoreCase(outsingle));
		assert(!in1.equalsIgnoreCase(in2));
		assert(out1==null || !out1.equalsIgnoreCase(out2));
		assert(out1==null || !out1.equalsIgnoreCase(outsingle));
	}

	public void process(){
		
		Timer t=new Timer();
		t.start();
		
		process2();
		
		t.stop();
		
		outstream.println("\nInput:                  \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		
		if(trimLeft || trimRight){
			outstream.println("Trimmed:                \t"+readsTrimmed+" reads ("+String.format("%.2f",readsTrimmed*100.0/readsIn)+"%) \t"+
					basesTrimmed+" bases ("+String.format("%.2f",basesTrimmed*100.0/basesIn)+"%)");
		}
		outstream.println("Result:                 \t"+readsOut+" reads ("+String.format("%.2f",readsOut*100.0/readsIn)+"%) \t"+
				basesOut+" bases ("+String.format("%.2f",basesOut*100.0/basesIn)+"%)");
		outstream.println("Pairs:                  \t"+pairsOut+" reads ("+String.format("%.2f",pairsOut*100.0/readsIn)+"%) \t"+
				pairBasesOut+" bases ("+String.format("%.2f",pairBasesOut*100.0/basesIn)+"%)");
		outstream.println("Singletons:             \t"+singlesOut+" reads ("+String.format("%.2f",singlesOut*100.0/readsIn)+"%) \t"+
				singleBasesOut+" bases ("+String.format("%.2f",singleBasesOut*100.0/basesIn)+"%)");
		
		double rpnano=readsIn/(double)(t.elapsed);
		double bpnano=basesIn/(double)(t.elapsed);

		String rpstring=(readsIn<100000 ? ""+readsIn : readsIn<100000000 ? (readsIn/1000)+"k" : (readsIn/1000000)+"m");
		String bpstring=(basesIn<100000 ? ""+basesIn : basesIn<100000000 ? (basesIn/1000)+"k" : (basesIn/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		if(showSpeed){
			outstream.println("\nTime:   \t\t\t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		}
		
		if(errorState){
			throw new RuntimeException("BBDuk terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void process2(){
		final Thread cristhread;
		final ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		final RTextOutputStream3 ros, rosb;
		final int buff=4;
		if(out1!=null){
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, false);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, false);
			ros=new RTextOutputStream3(ff1, ff2, buff, null, true);
			ros.start();
		}else{ros=null;}
		if(outsingle!=null){
			FileFormat ff1=FileFormat.testOutput(outsingle, FileFormat.FASTQ, null, true, overwrite, false);
			rosb=new RTextOutputStream3(ff1, null, buff, null, true);
			rosb.start();
		}else{rosb=null;}
		if(ros!=null || rosb!=null){
			outstream.println("Started output stream.");
		}
		
//		assert(false) : out1+", "+out2+", "+outsingle;
		if(fixPairs){
			process3_fixPairs(cris, ros, rosb);
		}else{
			process3(cris, ros, rosb);
		}

		
		ReadWrite.closeStreams(cris, ros, rosb);
	}
//	
//	private void process3_old(final ConcurrentReadStreamInterface cris, final RTextOutputStream3 ros, final RTextOutputStream3 rosb){
//
//		ListNum<Read> ln=cris.nextList();
//		ArrayList<Read> reads0=(ln!=null ? ln.list : null);
//		ArrayList<Read> single=(rosb==null ? null : new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
//		
//		while(reads0!=null && reads0.size()>0){
//			ArrayList<Read> reads=(ArrayList<Read>) reads0.clone();
//			int removed=0;
//			for(int i=0; i<reads.size(); i++){
//				Read r1=reads.get(i);
//				Read r2=r1.mate;
//				
//				readsIn++;
//				basesIn+=(r1.bases==null ? 0 : r1.bases.length);
//				if(r2!=null){
//					readsIn++;
//					basesIn+=(r2.bases==null ? 0 : r2.bases.length);
//				}
//				
//				{
//					if(trimLeft || trimRight){
//						if(r1!=null){
//							int x=TrimRead.trimFast(r1, trimLeft, trimRight, trimq, 1);
//							basesTrimmed+=x;
//							readsTrimmed+=(x>0 ? 1 : 0);
//						}
//						if(r2!=null){
//							int x=TrimRead.trimFast(r2, trimLeft, trimRight, trimq, 1);
//							basesTrimmed+=x;
//							readsTrimmed+=(x>0 ? 1 : 0);
//						}
//					}
//
//					final int rlen1=(r1==null ? -1 : r1.bases==null ? 0 : r1.bases.length);
//					final int rlen2=(r2==null ? -1 : r2.bases==null ? 0 : r2.bases.length);
//					
//					if(verbose){System.err.println("rlen1="+rlen1+", rlen2="+rlen2);}
//					
//					if(rlen1<minReadLength || rlen2<minReadLength){
//						reads.set(i, null);
//						removed++;
//						r1.mate=null;
//						if(r2!=null){
//							r2.mate=null;
//						}
//						
//						if(rlen1>=minReadLength){
//							single.add(r1);
//							singlesOut++;
//							singleBasesOut+=rlen1;
//						}
//						if(rlen2>=minReadLength){
//							single.add(r2);	
//							singlesOut++;
//							singleBasesOut+=rlen2;
//						}
//					}else{
//						if(r1!=null){
//							pairsOut++;
//							pairBasesOut+=rlen2;
//						}
//						if(r2!=null){
//							pairsOut++;
//							pairBasesOut+=rlen2;
//						}
//					}
//				}
//			}
//			
//			if(rosb!=null){
//				if(verbose){System.err.println("Adding "+single.size()+" to single out.");}
//				rosb.add(new ArrayList<Read>(single), ln.id);
//				single.clear();
//			}
//			
//			if(ros!=null){
//				if(removed>0){Tools.condenseStrict(reads);}
//				ArrayList<Read> x=new ArrayList<Read>(reads.size());
//				x.addAll(reads);
//				if(verbose){System.err.println("Adding "+x.size()+" to pair out.");}
//				ros.add(x, ln.id);
//			}
//			
//			cris.returnList(ln, ln.list.isEmpty());
//			ln=cris.nextList();
//			reads0=(ln!=null ? ln.list : null);
//		}
//		cris.returnList(ln, ln.list.isEmpty());
//		
//		readsOut+=singlesOut+pairsOut;
//		basesOut+=singleBasesOut+pairBasesOut;
//	}
	
	private void process3(final ConcurrentReadStreamInterface cris, final RTextOutputStream3 ros, final RTextOutputStream3 rosb){

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=ln.list;
		
		final ArrayList<Read> pairs=(ros==null ? null : new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
		final ArrayList<Read> singles=(rosb==null ? null : new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
		
		while(reads!=null && reads.size()>0){
			for(int i=0; i<reads.size(); i++){
				Read r1=reads.get(i);
				Read r2=r1.mate;
				processPair(r1, r2, pairs, singles);
			}
			
			cris.returnList(ln, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			
			if(rosb!=null){
				if(verbose){System.err.println("Adding "+singles.size()+" to single out.");}
				rosb.add(new ArrayList<Read>(singles), ln.id);
				singles.clear();
			}
			
			if(ros!=null){
				if(verbose){System.err.println("Adding "+pairs.size()+" to pair out.");}
				ros.add(new ArrayList<Read>(pairs), ln.id);
				pairs.clear();
			}
		}
		cris.returnList(ln, ln.list.isEmpty());
		
		readsOut+=singlesOut+pairsOut;
		basesOut+=singleBasesOut+pairBasesOut;
	}
	
	private void process3_fixPairs(final ConcurrentReadStreamInterface cris, final RTextOutputStream3 ros, final RTextOutputStream3 rosb){

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=ln.list;
		
		final ArrayList<Read> pairs=(ros==null ? null : new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
		final ArrayList<Read> singles=(rosb==null ? null : new ArrayList<Read>(Shared.READ_BUFFER_LENGTH));
		
		Read current=null, prev=null;
		
		while(reads!=null && reads.size()>0){
			for(int i=0; i<reads.size(); i++){
				
				current=reads.get(i);
//				if(verbose){System.err.println("Fetched "+current);}
				
				if(prev!=null){
					boolean b=FASTQ.testPairNames(prev, current);
					if(b){
						if(verbose){System.err.println("A");}
						processPair(prev, current, pairs, singles);
						prev=null;
						current=null;
					}else{
						if(verbose){System.err.println("B");}
						processPair(prev, null, null, singles);
						prev=null;
					}
				}
				prev=current;
				current=null;
			}
			
//			if(verbose){System.err.println("X\n"+current+"\n"+prev+"\n");}
			
			cris.returnList(ln, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			
			if((ln==null || reads==null || reads.isEmpty()) && prev!=null){ //Process last read
				boolean b=FASTQ.testPairNames(prev, current);
				if(b){
					if(verbose){System.err.println("C");}
					processPair(prev, current, pairs, singles);
					prev=null;
					current=null;
				}else{
					if(verbose){System.err.println("D");}
					processPair(prev, null, null, singles);
					prev=null;
				}
			}
			
			if(rosb!=null){
				if(verbose){System.err.println("Adding "+singles.size()+" to single out.");}
				rosb.add(new ArrayList<Read>(singles), ln.id);
				singles.clear();
			}
			
			if(ros!=null){
				if(verbose){System.err.println("Adding "+pairs.size()+" to pair out.");}
				ros.add(new ArrayList<Read>(pairs), ln.id);
				pairs.clear();
			}
		}
		cris.returnList(ln, ln.list.isEmpty());
		
		readsOut+=singlesOut+pairsOut;
		basesOut+=singleBasesOut+pairBasesOut;
	}
	
	
	private int processPair(Read r1, Read r2, ArrayList<Read> pairs, ArrayList<Read> singles){
		int removed=0;
		readsIn++;
		basesIn+=(r1.bases==null ? 0 : r1.bases.length);
		if(r2!=null){
			readsIn++;
			basesIn+=(r2.bases==null ? 0 : r2.bases.length);
		}
		
		if(trimLeft || trimRight){
			if(r1!=null){
				int x=TrimRead.trimFast(r1, trimLeft, trimRight, trimq, 1);
				basesTrimmed+=x;
				readsTrimmed+=(x>0 ? 1 : 0);
			}
			if(r2!=null){
				int x=TrimRead.trimFast(r2, trimLeft, trimRight, trimq, 1);
				basesTrimmed+=x;
				readsTrimmed+=(x>0 ? 1 : 0);
			}
		}
		final int rlen1=(r1==null ? -1 : r1.bases==null ? 0 : r1.bases.length);
		final int rlen2=(r2==null ? -1 : r2.bases==null ? 0 : r2.bases.length);
		if(verbose){System.err.println("rlen="+rlen1+", rlen2="+rlen2);}
		
		if(rlen1>=minReadLength && rlen2>=minReadLength){
			if(verbose){System.err.println("Sending to pair out:\t"+r1.id+"\t"+r2.id);}
			r1.mate=r2;
			r2.mate=r1;
			r1.setPairnum(0);
			r2.setPairnum(1);
			if(pairs!=null){pairs.add(r1);}					
			pairsOut+=2;
			pairBasesOut+=(rlen1+rlen2);
		}else if(rlen1>=minReadLength){
			if(verbose){System.err.println("Sending r1 to single out:\t"+r1.id+"\t"+(r2==null ? "*" : r2.id));}
			r1.mate=null;
			r1.setPairnum(0);
			if(singles!=null){singles.add(r1);}
			singlesOut++;
			singleBasesOut+=rlen1;
			if(r2!=null){removed++;}
		}else if(rlen2>=minReadLength){
			if(verbose){System.err.println("Sending r2 to single out:\t"+(r1==null ? "*" : r1.id)+"\t"+r2.id);}
			r2.mate=null;
			r2.setPairnum(0);
			if(singles!=null){singles.add(r2);}
			singlesOut++;
			singleBasesOut+=rlen2;
			if(r1!=null){removed++;}
		}else{
			if(verbose){System.err.println("Removed both reads:\t"+(r1==null ? "*" : r1.id)+"\t"+(r2==null ? "*" : r2.id));}
			if(r1!=null){removed++;}
			if(r2!=null){removed++;}
		}
		return removed;
	}
	
	
	private String in1=null, in2=null;
	private String out1=null, out2=null;
	private String outsingle=null;
	private long maxReads=-1;
	public boolean errorState=false;

	long readsIn=0;
	long basesIn=0;
	long readsOut=0;
	long basesOut=0;
	long pairsOut=0;
	long pairBasesOut=0;
	long singlesOut=0;
	long singleBasesOut=0;
	long readsTrimmed=0;
	long basesTrimmed=0;

	private byte trimq=4;
	private int minReadLength=20;
	private final boolean trimLeft, trimRight;

	private final boolean EA;
	private final boolean fixPairs;
	
	private static PrintStream outstream=System.err;
	public static boolean overwrite=false;
	public static boolean showSpeed=true;
	public static boolean verbose=false;
	
}
