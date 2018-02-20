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
import stream.SamLine;

import dna.Data;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import fileIO.FileFormat;

import align2.ListNum;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;

/**
 * @author Brian Bushnell
 * @date Sep 11, 2012
 *
 */
public class ReformatReads {

	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		ReformatReads rr=new ReformatReads(args);
		rr.process(t);
	}
	
	public ReformatReads(String[] args){
		if(args==null || args.length==0){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");

		boolean setxs=false, setintron=false;
		boolean setInterleaved=false; //Whether it was explicitly set.

		FastaReadInputStream.SPLIT_READS=false;
		stream.FastaReadInputStream.MIN_READ_LEN=1;
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.READ_BUFFER_NUM_BUFFERS=Tools.min(8, Shared.READ_BUFFER_NUM_BUFFERS);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("null") || a.equals(in2)){
				// do nothing
			}else if(a.equals("passes")){
				assert(false) : "'passes' is disabled.";
//				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
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
			}else if(a.equals("sample") || a.equals("samplereads") || a.equals("samplereadstarget") || a.equals("srt")){
				sampleReadsTarget=Long.parseLong(b);
				sampleReadsExact=(sampleReadsTarget>0);
			}else if(a.equals("samplebases") || a.equals("samplebasestarget") || a.equals("sbt")){
				sampleBasesTarget=Long.parseLong(b);
				sampleBasesExact=(sampleBasesTarget>0);
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("bf1")){
				ByteFile.FORCE_MODE_BF1=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF2=!ByteFile.FORCE_MODE_BF1;
			}else if(a.equals("bf2")){
				ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b);
				ByteFile.FORCE_MODE_BF1=!ByteFile.FORCE_MODE_BF2;
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				out1=b;
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
			}else if(a.equals("extin")){
				extin=b;
			}else if(a.equals("extout")){
				extout=b;
			}else if(a.equals("trd") || a.equals("trc") || a.equals("trimreaddescription")){
				Shared.TRIM_READ_COMMENTS=Tools.parseBoolean(b);
			}else if(a.equals("parsecustom")){
				parsecustom=Tools.parseBoolean(b);
			}else if(a.equals("testsize")){
				testsize=Tools.parseBoolean(b);
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
			}else if(a.equals("verifyinterleaved") || a.equals("verifyinterleaving") || a.equals("vint")){
				verifyinterleaving=Tools.parseBoolean(b);
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
			}else if(a.equals("rcompmate") || a.equals("rcm")){
				RCOMPMATE=Tools.parseBoolean(b);
				outstream.println("Set RCOMPMATE to "+RCOMPMATE);
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
				assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
			}else if(a.equals("sampleseed")){
				sampleseed=Long.parseLong(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				int x=Integer.parseInt(b);
				stream.FastaReadInputStream.MIN_READ_LEN=(x>0 ? x : Integer.MAX_VALUE);
			}else if(a.equals("samversion") || a.equals("samv") || a.equals("sam")){
				SamLine.VERSION=Float.parseFloat(b);
			}else if(a.equals("mdtag") || a.equals("md")){
				SamLine.MAKE_MD_TAG=Tools.parseBoolean(b);
			}else if(a.equals("xstag") || a.equals("xs")){
				SamLine.MAKE_XS_TAG=true;
				if(b!=null){
					b=b.toLowerCase();
					if(b.startsWith("fr-")){b=b.substring(3);}
					if(b.equals("ss") || b.equals("secondstrand")){
						SamLine.XS_SECONDSTRAND=true;
					}else if(b.equals("fs") || b.equals("firststrand")){
						SamLine.XS_SECONDSTRAND=false;
					}else if(b.equals("us") || b.equals("unstranded")){
						SamLine.XS_SECONDSTRAND=false;
					}else{
						SamLine.MAKE_XS_TAG=Tools.parseBoolean(b);
					}
				}
				setxs=true;
			}else if(a.equals("intronlen") || a.equals("intronlength")){
				SamLine.INTRON_LIMIT=Integer.parseInt(b);
				setintron=true;
			}else if(a.equals("idtag")){
				SamLine.MAKE_IDENTITY_TAG=Tools.parseBoolean(b);
			}else if(a.equals("xmtag") || a.equals("xm")){
				SamLine.MAKE_XM_TAG=Tools.parseBoolean(b);
			}else if(a.equals("stoptag")){
				SamLine.MAKE_STOP_TAG=Tools.parseBoolean(b);
			}else if(a.equals("ftl") || a.equals("forcetrimleft")){
				forceTrimLeft=Integer.parseInt(b);
			}else if(a.equals("ftr") || a.equals("forcetrimright")){
				forceTrimRight=Integer.parseInt(b);
//			}else if(a.equals("ftlb") || a.equals("forcetrimleftby")){
//				forceTrimLeftBy=Integer.parseInt(b);
//			}else if(a.equals("ftrb") || a.equals("forcetrimrightby")){
//				forceTrimRightBy=Integer.parseInt(b);
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
			}else if(a.equals("trimright") || a.equals("qtrimright")){
				trimRight_=Tools.parseBoolean(b);
			}else if(a.equals("trimleft") || a.equals("qtrimleft")){
				trimLeft_=Tools.parseBoolean(b);
			}else if(a.equals("trimq") || a.equals("trimquality")){
				trimq=Byte.parseByte(b);
			}else if(a.equals("q102matrix") || a.equals("q102m")){
				CalcTrueQuality.q102matrix=b;
			}else if(a.equals("qbpmatrix") || a.equals("bqpm")){
				CalcTrueQuality.qbpmatrix=b;
			}else if(a.equals("adjustquality") || a.equals("adjq")){
				TrimRead.ADJUST_QUALITY=Tools.parseBoolean(b);
			}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
				minReadLength=Integer.parseInt(b);
			}else if(a.equals("mlf") || a.equals("minlenfrac") || a.equals("minlenfraction") || a.equals("minlengthfraction")){
				minLenFraction=Float.parseFloat(b);
			}else if(a.equals("minavgquality") || a.equals("maq")){
				minAvgQuality=Byte.parseByte(b);
			}else if(in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(out1==null && i==1 && !arg.contains("=")){
				out1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(TrimRead.ADJUST_QUALITY){CalcTrueQuality.initializeMatrices();}
		
		if(setxs && !setintron){SamLine.INTRON_LIMIT=10;}
		qtrim=trimLeft_||trimRight_;

		if(verifyinterleaving){
			setInterleaved=true;
//			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true;
		}
		if(fixinterleaving){
			setInterleaved=true;
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true; //? Not sure which way to go
			throw new RuntimeException("fixinterleaving: TODO");
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
//		if(maxReads!=-1){ReadWrite.USE_GUNZIP=ReadWrite.USE_UNPIGZ=false;}
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.THREADS>2){
//			if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null){
			if(out2!=null){
				printOptions();
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
		
//		assert(false) : "qin="+qin+", qout="+qout+", FASTQ.ASCII_OFFSET="+FASTQ.ASCII_OFFSET+"\n"+
//			", FASTQ.ASCII_OFFSET_OUT="+FASTQ.ASCII_OFFSET_OUT+", FASTQ.DETECT_QUALITY="+FASTQ.DETECT_QUALITY;
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, false);
//		extsOut=(out1==null ? null : FileFormat.testFormat(out1, false));
//		outsam=(extsOut!=null && (extsOut[0]==FileFormat.SAM || extsOut[0]==FileFormat.BAM));
//		outbread=(extsOut!=null && extsOut[0]==FileFormat.BREAD);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

//		System.err.println("\n"+ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+Data.PIGZ()+", "+Data.UNPIGZ()+", "+ffin1+"\n");
//		assert(false) : ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+Data.PIGZ()+", "+Data.UNPIGZ()+", "+ffin1;
		
//		extsIn=(in1==null ? null : FileFormat.testFormat(in1, false));
//		insam=(extsIn!=null && (extsIn[0]==FileFormat.SAM || extsIn[0]==FileFormat.BAM));

		if(ffin1!=null && ffout1!=null && ffin1.samOrBam()){
			if(ffout1.samOrBam()){
				useSharedHeader=true;
				SamLine.CONVERT_CIGAR_TO_MATCH=true;
			}else if(ffout1.bread()){
				SamLine.CONVERT_CIGAR_TO_MATCH=true;
			}
		}
	}
	
	void process(Timer t){
		
		long readsRemaining=0;
		long basesRemaining=0;
		
		if(sampleReadsExact || sampleBasesExact){
			long[] counts=countReads(in1, in2, maxReads);
			readsRemaining=counts[0];
			basesRemaining=counts[2];
			setSampleSeed(sampleseed);
		}
		
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
//			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
//			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, colorspace, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){System.err.println("Started cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}

		RTextOutputStream3 ros=null;
		if(out1!=null){
			final int buff=4;
			
//			if(!fq && !fa && !bread && !sam){
//				outstream.println("Unspecified output format; defaulting to uncompressed fastq.");
//				fq=true;
//			}
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
//			System.err.println("Calling RTextOutputStream3 with out1="+out1+", out2="+out2+", qfout1="+qfout1+", qfout2="+qfout2);
			ros=new RTextOutputStream3(ffout1, ffout2, qfout1, qfout2, buff, null, useSharedHeader);
			ros.start();
		}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		long basesTrimmedT=0;
		long readsTrimmedT=0;
		
		long lowqBasesT=0;
		long lowqReadsT=0;
		
		long readShortDiscardsT=0;
		long baseShortDiscardsT=0;
		
//		for(int pass=1; pass<=passes; pass++){
////			outstream.println("pass="+pass);
//			if(pass>1){
//				cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, colorspace, useSharedHeader, true, in1, in2);
//				cris.setSampleRate(samplerate, sampleseed);
//				cristhread=new Thread(cris);
//				cristhread.start();
//			}
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=(r1.bases==null ? 0 : r1.bases.length);
					final int initialLength2=(r2==null ? 0 : r2.bases==null ? 0 : r2.bases.length);

					final int minlen1=(int)Tools.max(initialLength1*minLenFraction, minReadLength);
					final int minlen2=(int)Tools.max(initialLength2*minLenFraction, minReadLength);
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
						if(RCOMPMATE){r2.reverseComplement();}
					}
					
					if(verifyinterleaving){
						String s1=r1==null ? null : r1.id;
						String s2=r2==null ? null : r2.id;
						boolean b=FASTQ.testPairNames(s1, s2);
						if(!b){
							outstream.println("Names do not appear to be correctly paired.\n"+s1+"\n"+s2+"\n");
							ReadWrite.closeStreams(cris, ros);
							System.exit(1);
						}
					}
					
					boolean remove=false;

					if(minAvgQuality>0){
						final int a=(r1==null ? 50 : r1.avgQuality());
						final int b=(r2==null ? 50 : r2.avgQuality());
						remove=(remove || a<minAvgQuality || b<minAvgQuality);
						if(remove){
							if(r1!=null){
								lowqBasesT+=r1.bases.length;
								lowqReadsT++;
							}
							if(r2!=null){
								lowqBasesT+=r2.bases.length;
								lowqReadsT++;
							}
						}
					}
					
					if(!remove && (forceTrimLeft>0 || forceTrimRight>0)){
						if(r1!=null){
							int x=TrimRead.trimToPosition(r1, forceTrimLeft>0 ? forceTrimLeft : 0, forceTrimRight>0 ? forceTrimRight : r1.bases.length, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2!=null){
							int x=TrimRead.trimToPosition(r2, forceTrimLeft>0 ? forceTrimLeft : 0, forceTrimRight>0 ? forceTrimRight : r2.bases.length, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
					}
					
					if(qtrim && !remove){
						if(r1!=null){
							int x=TrimRead.trimFast(r1, trimLeft_, trimRight_, trimq, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2!=null){
							int x=TrimRead.trimFast(r2, trimLeft_, trimRight_, trimq, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
					}
					
					if((minlen1>0 || minlen2>0) && !remove){
						int rlen=(r1==null || r1.bases==null ? 0 : r1.bases.length);
						int rlen2=(r2==null || r2.bases==null ? 0 : r2.bases.length);
						if(rlen<minlen1 && (r2==null || rlen2<minlen2)){
							remove=true;
							readShortDiscardsT+=(r2==null ? 1 : 2);
							baseShortDiscardsT+=(rlen+rlen2);
						}
					}
					
					if(remove){reads.set(idx, null);}
				}
				
//				assert(false) : sampleReadsExact+", "+sampleBasesExact;
				if(sampleReadsExact || sampleBasesExact){
					ArrayList<Read> reads2=new ArrayList<Read>();
					if(sampleReadsExact){
						for(Read r : reads){
							if(r!=null){
								assert(readsRemaining>0) : readsRemaining;
								double prob=sampleReadsTarget/(double)(readsRemaining);
//								System.err.println("sampleReadsTarget="+sampleReadsTarget+", readsRemaining="+readsRemaining+", prob="+prob);
								if(randy.nextDouble()<prob){
									reads2.add(r);
									sampleReadsTarget--;
								}
							}
							readsRemaining--;
						}
					}else if(sampleBasesExact){
						for(Read r : reads){
							if(r!=null){
								assert(basesRemaining>0) : basesRemaining;
								int bases=r.bases.length+(r.mate==null ? 0 : r.mate.bases.length);
								double prob=sampleBasesTarget/(double)(basesRemaining);
								if(randy.nextDouble()<prob){
									reads2.add(r);
									sampleBasesTarget-=bases;
								}
								basesRemaining-=bases;
							}
						}
					}
					if(ros!=null){ros.add(reads2, ln.id);}
				}else{
					if(ros!=null){ros.add(reads, ln.id);}
				}

				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
//		System.err.println(cris.errorState()+", "+(ros==null ? "null" : (ros.errorState()+", "+ros.finishedSuccessfully())));
//		if(ros!=null){
//			ReadStreamWriter rs1=ros.getRS1();
//			ReadStreamWriter rs2=ros.getRS2();
//			System.err.println(rs1==null ? "null" : rs1.finishedSuccessfully());
//			System.err.println(rs2==null ? "null" : rs2.finishedSuccessfully());
//		}
//		assert(false);
		
		t.stop();

		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}

		if(qtrim || forceTrimLeft>0 || forceTrimRight>0){
			outstream.println("QTrimmed:               \t"+readsTrimmedT+" reads ("+String.format("%.2f",readsTrimmedT*100.0/readsProcessed)+"%) \t"+
					basesTrimmedT+" bases ("+String.format("%.2f",basesTrimmedT*100.0/basesProcessed)+"%)");
		}else if(minReadLength>0){
			outstream.println("Short Read Discards:    \t"+readShortDiscardsT+" reads ("+String.format("%.2f",readShortDiscardsT*100.0/readsProcessed)+"%) \t"+
					baseShortDiscardsT+" bases ("+String.format("%.2f",baseShortDiscardsT*100.0/basesProcessed)+"%)");
		}
		if(minAvgQuality>0){
			outstream.println("Low quality discards:   \t"+lowqReadsT+" reads ("+String.format("%.2f",lowqReadsT*100.0/readsProcessed)+"%) \t"+
					lowqBasesT+" bases ("+String.format("%.2f",lowqBasesT*100.0/basesProcessed)+"%)");
		}
//		if(qtrim || minAvgQuality>0){
//			outstream.println("Result:                 \t"+readsOut+" reads ("+String.format("%.2f",readsOut*100.0/readsProcessed)+"%) \t"+
//					basesOut+" bases ("+String.format("%.2f",basesOut*100.0/basesProcessed)+"%)");
//		}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		if(testsize){
			long bytesProcessed=(new File(in1).length()+(in2==null ? 0 : new File(in2).length())+
					(qfin1==null ? 0 : new File(qfin1).length())+(qfin2==null ? 0 : new File(qfin2).length()));//*passes
			double xpnano=bytesProcessed/(double)(t.elapsed);
			String xpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
			while(xpstring.length()<8){xpstring=" "+xpstring;}
			outstream.println("Bytes Processed:    "+xpstring+" \t"+String.format("%.2fm bytes/sec", xpnano*1000));
		}
		
		if(verifyinterleaving){
			outstream.println("Names appear to be correctly paired.");
		}
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private long[] countReads(String fname1, String fname2, long maxReads){
		{
			String x=fname1.toLowerCase();
			if((x.equals("stdin") || x.startsWith("stdin.")) && !new File(fname1).exists()){
				throw new RuntimeException("Can't precount reads from standard in, only from a file.");
			}
		}
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, colorspace, false, ffin1, ffin2, null, null);
			if(verbose){System.err.println("Counting Reads");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		long count=0, count2=0, bases=0;
		
		while(reads!=null && reads.size()>0){
			count+=reads.size();
			for(Read r : reads){
				bases+=r.bases.length;
				count2++;
				if(r.mate!=null){
					bases+=r.mate.bases.length;
					count2++;
				}
			}
			cris.returnList(ln, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln, ln.list.isEmpty());
		errorState|=ReadWrite.closeStream(cris);
		return new long[] {count, count2, bases};
	}
	
	private void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("java -ea -Xmx512m -cp <path> jgi.ReformatReads in=<infile> in2=<infile2> out=<outfile> out2=<outfile2>");
		outstream.println("\nin2 and out2 are optional.  \nIf input is paired and there is only one output file, it will be written interleaved.\n");
		outstream.println("Other parameters and their defaults:\n");
		outstream.println("overwrite=false  \tOverwrites files that already exist");
		outstream.println("ziplevel=5       \tSet compression level, 1 (low) to 9 (max)");
		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
		outstream.println("fastawrap=100    \tLength of lines in fasta output");
		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
	}
	
	
	public void setSampleSeed(long seed){
		if(seed>-1){
			randy=new java.util.Random(seed);
		}else{
			randy=new java.util.Random();
		}
	}
	
	
	/*--------------------------------------------------------------*/
	
	public boolean errorState=false;
	
	private String in1=null;
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;

	private boolean parsecustom=false;
	private boolean testsize=false;
	private boolean overwrite=false;
	private boolean RCOMPMATE=false;
	private boolean verifyinterleaving=false;
	private boolean fixinterleaving=false;
	private boolean colorspace=false;
	
	private long maxReads=-1;
//	private int passes=1;
	private float samplerate=1f;
	private long sampleseed=-1;
	
	private boolean trimRight_=false;
	private boolean trimLeft_=false;
	private byte trimq=4;
	private byte minAvgQuality=0;
	private int minReadLength=0;
	private float minLenFraction=0;

	private int forceTrimLeft=-1;
	private int forceTrimRight=-1;
	
	private byte qin=-1;
	private byte qout=-1;
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	private boolean useSharedHeader;
	private final boolean qtrim;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	private boolean sampleReadsExact=false;
	private long sampleReadsTarget=0;
	private boolean sampleBasesExact=false;
	private long sampleBasesTarget=0;
	public static boolean verbose=false;
	
	private java.util.Random randy;
	
}
