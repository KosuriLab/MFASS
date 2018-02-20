package jgi;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

import kmer.KCountArray;
import kmer.KmerCount7MT;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.RTextOutputStream3;
import stream.Read;
import stream.SamLine;

import dna.AminoAcid;
import dna.Data;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.FileFormat;

import align2.ListNum;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;

/**
 * @author Brian Bushnell
 * @date Aug 20, 2012
 *
 */
public class ErrorCorrectMT extends Thread{
	
	public static void main(String[] args){
		System.err.println("Executing "+(new Object() { }.getClass().getEnclosingClass().getName())+" "+Arrays.toString(args)+"\n");
		
		if(args.length<1){throw new RuntimeException("No parameters.");}
		
		String reads1=args[0];
		String reads2=(args.length>1 ? args[1] : null);
		if(reads2!=null && "null".equalsIgnoreCase(reads2)){reads2=null;}
		
		{
			{
				File f=new File(reads1);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(reads1+" does not exist.");}
			}
			if(reads2!=null){
				File f=new File(reads2);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(reads2+" does not exist.");}
				if(reads1.equalsIgnoreCase(reads2)){
					throw new RuntimeException("Both input files are the same.");
				}
			}
		}
		
		int k=31;
		int cbits=2;
		int gap=0;
		int hashes=3;
		int thresh1=2;
		int thresh2=2;
//		int matrixbits=-1;
		long cells=-1;
		long maxReads=-1;
		int buildpasses=2;
		long tablereads=-1; //How many reads to process when building the hashtable
		int buildStepsize=4;
		String output=null;
		boolean ordered=true;
		boolean overwrite=true;
		int threads=-1;
		
		boolean auto=true;
		
		FastaReadInputStream.TARGET_READ_LEN=Integer.MAX_VALUE;
		FASTQ.PARSE_CUSTOM=false;
		
		List<String> extra=null;
		
		long memory=Runtime.getRuntime().maxMemory();
		long tmemory=Runtime.getRuntime().totalMemory();
//		assert(false) : memory+", "+tmemory;
		
		for(int i=2; i<args.length; i++){
			if(args[i]==null){args[i]="null";}
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=(split.length>1 ? split[1] : "true");
			if("null".equalsIgnoreCase(b)){b=null;}
			
			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.equals("initialthresh") || a.equals("thresh1")){
				thresh1=Integer.parseInt(b);
			}else if(a.equals("thresh") || a.equals("thresh2") || a.equals("threshgood")){
				thresh2=Integer.parseInt(b);
			}else if(a.equals("threshbad")){
				THRESH_BAD=Integer.parseInt(b);
			}else if(a.startsWith("gap")){
				gap=Integer.parseInt(b);
			}else if(a.startsWith("matrixbits")){
				int matrixbits=Integer.parseInt(b);
				assert(matrixbits<63);
				cells=1L<<matrixbits;
			}else if(a.startsWith("cells")){
				cells=Tools.parseKMG(b);
			}else if(a.startsWith("minq")){
				KmerCount7MT.minQuality=Byte.parseByte(b);
			}else if(a.startsWith("minprob")){
				KmerCount7MT.minProb=Float.parseFloat(b);
			}else if(a.startsWith("hashes") || a.startsWith("multihash")){
				hashes=Integer.parseInt(b);
			}else if(a.equals("ascii") || a.equals("asciioffset")){
				byte q=Byte.parseByte(b);
				Data.sysout.println("Setting quality ASCII offset to "+q);
				if(q!=33 && q!=64){
					Data.sysout.println("Warning! Expected ASCII offsets are 33 (Sanger) or 64 (old Illumina), not "+q);
				}
				FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT=q;
				FASTQ.DETECT_QUALITY=false;
			}else if(a.equals("asciiin") || a.equals("qualityin") || a.equals("qualin") || a.equals("qin")){
				if(args[i].indexOf('=')>0){
					byte ascii_offset=Byte.parseByte(b);
					FASTQ.ASCII_OFFSET=ascii_offset;
					Data.sysout.println("Set fastq input ASCII offset to "+FASTQ.ASCII_OFFSET);
					FASTQ.DETECT_QUALITY=false;
				}else{
					try {
						a=a.replace("-", "");
						a=a.replace("ascii", "");
						byte q=Byte.parseByte(a);
						Data.sysout.println("Setting quality ASCII offset to "+q);
						if(q!=33 && q!=64){
							Data.sysout.println("Warning! Expected ASCII offsets are 33 (Sanger) or 64 (old Illumina), not "+q);
						}
						FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT=q;
						FASTQ.DETECT_QUALITY=false;
					} catch (NumberFormatException e) {
						System.err.println("Warning! Could not parse "+args[i]);
						e.printStackTrace();
					}
				}
			}else if(a.startsWith("maxerrors")){
				ERROR_CORRECTION_LIMIT=Integer.parseInt(b);
			}else if(a.startsWith("maxburst")){
				MAX_ERROR_BURST=Integer.parseInt(b);
			}else if(a.startsWith("passes")){
				buildpasses=Integer.parseInt(b);
			}else if(a.startsWith("stepsize") || a.startsWith("buildstepsize")){
				buildStepsize=Integer.parseInt(b);
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				threads=Integer.parseInt(b);
			}else if(a.startsWith("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.startsWith("tablereads") || a.startsWith("buildreads")){
				tablereads=Long.parseLong(b);
			}else if(a.startsWith("build") || a.startsWith("genome")){
				Data.setGenome(Integer.parseInt(b));
				Data.sysout.println("Set genome to "+Data.GENOME_BUILD);
			}else if(a.equals("outputinfo") || a.startsWith("info")){
				OUTPUT_INFO_ONLY=Tools.parseBoolean(b);
			}else if(a.startsWith("out")){
				output=b;
			}else if(a.startsWith("verbose")){
				KCountArray.verbose=Tools.parseBoolean(b);
//				verbose=KCountArray.verbose=Tools.parseBoolean(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered=Tools.parseBoolean(b);
			}else if(a.startsWith("dontoutputbadreads")){
				DONT_OUTPUT_BAD_READS=Tools.parseBoolean(b);
			}else if(a.startsWith("dontoutputbadpairs")){
				DONT_OUTPUT_BAD_PAIRS=Tools.parseBoolean(b);
			}else if(a.equals("trim") || a.equals("qtrim")){
				TRIM_LEFT=TRIM_RIGHT=Tools.parseBoolean(b);
			}else if(a.equals("optitrim") || a.equals("otf") || a.equals("otm")){
				if(b!=null && (b.charAt(0)=='.' || Character.isDigit(b.charAt(0)))){
					TrimRead.optimalMode=true;
					TrimRead.optimalBias=Float.parseFloat(b);
					assert(TrimRead.optimalBias>=0 && TrimRead.optimalBias<1);
				}else{
					TrimRead.optimalMode=Tools.parseBoolean(b);
				}
			}else if(a.equals("trimleft")){
				TRIM_LEFT=Tools.parseBoolean(b);
			}else if(a.equals("trimright")){
				TRIM_RIGHT=Tools.parseBoolean(b);
			}else if(a.startsWith("trimq")){
				TRIM_QUAL=MAX_TRIM_QUAL=Byte.parseByte(b);
			}else if(a.startsWith("mintrimq")){
				TRIM_QUAL=Byte.parseByte(b);
			}else if(a.startsWith("maxtrimq")){
				MAX_TRIM_QUAL=Byte.parseByte(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("auto") || a.equals("automatic")){
				auto=Tools.parseBoolean(b);
			}else if(a.equals("trybothsides")){
				TRY_BOTH_SIDES=Tools.parseBoolean(b);
			}else if(a.equals("onlycorrectn")){
				ONLY_CORRECT_N=Tools.parseBoolean(b);
			}else if(a.startsWith("parsecustom")){
				FASTQ.PARSE_CUSTOM=Tools.parseBoolean(b);
			}else if(a.equals("testinterleaved")){
				FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
				Data.sysout.println("Set TEST_INTERLEAVED to "+FASTQ.TEST_INTERLEAVED);
			}else if(a.equals("forceinterleaved")){
				FASTQ.FORCE_INTERLEAVED=Tools.parseBoolean(b);
				Data.sysout.println("Set FORCE_INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else if(a.equals("interleaved") || a.equals("int")){
				if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
				else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
					Data.sysout.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}else if(a.startsWith("fastareadlen")){
				FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
				FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
			}else if(a.startsWith("fastaminread") || a.startsWith("fastaminlen")){
				FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(a.startsWith("canonical")){
				CANONICAL=KmerCount7MT.CANONICAL=Tools.parseBoolean(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
			}else if(a.startsWith("extra")){
				if(b!=null && !b.equalsIgnoreCase("null")){
					extra=Arrays.asList(b.split(","));
				}
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		assert(FastaReadInputStream.settingsOK());
		if(k>31){CANONICAL=KmerCount7MT.CANONICAL=false;}
		assert(CANONICAL==KmerCount7MT.CANONICAL);
		
		assert(THRESH_BAD<thresh2);
		
		if(OUTPUT_INFO_ONLY){DONT_OUTPUT_BAD_READS=DONT_OUTPUT_BAD_PAIRS=false;}
		
		if(extra!=null){
			for(String s : extra){
				File f=new File(s);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(s+" does not exist.");}
				assert(!s.equalsIgnoreCase(reads1) && (reads2==null || !s.equalsIgnoreCase(reads2))) : "\nInput file "+s+" should not be included as an extra file.\n";
			}
		}
		
		Data.sysout.println("ForceInterleaved = "+FASTQ.FORCE_INTERLEAVED);
		
//		assert(false) : reads1+", "+reads2+", "+output;
//		if(FASTQ.FORCE_INTERLEAVED && in2==null){
//			Data.sysout.println()
//		}
		
		if(threads>0){
			THREADS=threads;
		}else{
			THREADS=Data.LOGICAL_PROCESSORS;
		}
		assert(THREADS>0 && THREADS<100000);
		
		long cells2=cells;
		if(auto && cells==-1){
			final long usable=(long)Tools.max(((memory-16000000)*.75), memory*0.45);
			long mem=usable;
			if(buildpasses>1){mem/=2;}
			cells=(mem*8)/cbits;
			cells2=cells;
			
//			long tablebytes=((1L<<matrixbits)*cbits)/8;
//			if(tablebytes*3<usable){matrixbits++;}
//			Data.sysout.println(tablebytes/1000000+", "+usable/1000000+", "+(tablebytes*3)/1000000);
		}else if(cells==-1){
			cells=1L<<34;
			cells2=cells;
		}
		
		if(k<32 && cells>(1L<<(2*k))){cells=cells2=(1L<<(2*k));}
		
		Data.sysout.println("\nSettings:");
		Data.sysout.println("threads:          \t"+THREADS);
		Data.sysout.println("k:                \t"+k);
		Data.sysout.println("cbits:            \t"+cbits);
		Data.sysout.println("cells:            \t"+Tools.toKMG(cells));
//		if(buildpasses>1){Data.sysout.println("cells2:           \t"+Tools.toKMG(cells2));}
		Data.sysout.println("hashes:           \t"+hashes);
		Data.sysout.println("passes:           \t"+buildpasses);
		Data.sysout.println("maxerrors:        \t"+ERROR_CORRECTION_LIMIT);
		Data.sysout.println("maxburst:         \t"+MAX_ERROR_BURST);
		if(buildpasses>1){
			Data.sysout.println("thresh1:          \t"+thresh1);
			Data.sysout.println("thresh2:          \t"+thresh2);
		}else{
			Data.sysout.println("thresh:           \t"+thresh2);
		}
		Data.sysout.println("output bad reads: \t"+(!DONT_OUTPUT_BAD_READS));
		Data.sysout.println();
		
//		KmerCount7MT.THREADS=Tools.max(THREADS/2, KmerCount7MT.THREADS);  //Seems like 4 is actually optimal...
		
		FastaReadInputStream.MIN_READ_LEN=k;
		
		if(DONT_OUTPUT_BAD_PAIRS){DONT_OUTPUT_BAD_READS=true;}
		
		Timer t=new Timer();
		t.start();
//		assert(false) : cells+", "+cells2;
		KCountArray kca=makeTable(reads1, reads2, extra, k, cbits, gap, hashes, buildpasses, cells, cells2, tablereads, buildStepsize, thresh1, thresh2);
		
		long bases=detect(reads1, reads2, kca, k, thresh2, maxReads, output, ordered, overwrite);
		t.stop();
		Data.sysout.println("Total time:      \t"+t+"   \t"+String.format("%.2f", bases*1000000.0/(t.elapsed))+" kb/sec");
		
	}
	
	public static KCountArray makeTable(String reads1, String reads2, Iterable<String> extra, int k, int cbits, int gap, int hashes, int buildpasses, long cells1,
			long cells2, long maxreads, int stepsize, int thresh1, int thresh2){
		
		long extraMaxreads=-1;
		
		Timer thash=new Timer();
		
		KmerCount7MT.maxReads=maxreads;
		int kbits=Tools.max(2, Tools.min(2*k, 62));

		long mcells=(buildpasses&1)==1 ? cells1 : cells2;
		
		thash.start();
//		Data.sysout.println("kbits="+(kbits)+" -> "+(1L<<kbits)+", matrixbits="+(matrixbits)+" -> "+(1L<<matrixbits)+", cbits="+cbits+", gap="+gap+", hashes="+hashes);
		KCountArray kca=KCountArray.makeNew(1L<<kbits, mcells, cbits, gap, hashes);
		
		assert(gap==0) : "TODO";
		if(buildpasses==1){

			Data.sysout.println("Running pass "+1+" with k="+k+", cbits="+cbits+", cells="+Tools.toKMG(kca.cells)+", thresh="+thresh2+", step="+stepsize+", conservative="+true);
//			KmerCount7MT.countFastq(reads1, reads2, k, cbits, gap, true, kca);
			try {
				KmerCount7MT.countFastq(reads1, reads2, k, cbits, true, kca, null, maxreads, thresh2, stepsize, true);
				if(extra!=null){
					for(String s : extra){
						KmerCount7MT.countFastq(s, null, k, cbits, true, kca, null, extraMaxreads, thresh2, stepsize, true);
					}
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			kca.shutdown();
			
		}else{
			assert(buildpasses>1);
			KCountArray trusted=null;
			for(int i=1; i<buildpasses; i++){
				boolean conservative=i>2;// /*or, alternately, (trusted==null || trusted.capacity()>0.3)
				int step=(stepsize==1 ? 1 : stepsize+i%2);
//				if(!conservative){step=(step+3)/4;}
				if(!conservative){step=Tools.min(3, (step+3)/4);}
				
				Data.sysout.println("Running pass "+i+" with k="+k+", cbits="+cbits+", cells="+Tools.toKMG(kca.cells)+", thresh="+thresh1+", step="+step+", conservative="+conservative);
				try {
					KmerCount7MT.countFastq(reads1, reads2, k, cbits, true, kca, trusted, maxreads, thresh1, step, conservative);
					if(extra!=null){
						for(String s : extra){
							KmerCount7MT.countFastq(s, null, k, cbits, true, kca, trusted, extraMaxreads, thresh1, step, conservative);
						}
					}
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				kca.shutdown();
				Data.sysout.println("Trusted:   \t"+kca.toShortString());
				trusted=kca;
//				assert(false) : cells2+", "+cells1;
				kca=KCountArray.makeNew(1L<<kbits, ((i&1)==(buildpasses&1) ? cells2 : cells1), cbits, gap, hashes);
				
			}
			
			Data.sysout.println("Running pass "+buildpasses+" with k="+k+", cbits="+cbits+", cells="+Tools.toKMG(kca.cells)+
					", thresh="+thresh2+", step="+stepsize+", conservative="+true);
			try {
				KmerCount7MT.countFastq(reads1, reads2, k, cbits, true, kca, trusted, maxreads, thresh2, stepsize, true);
				if(extra!=null){
					for(String s : extra){
						KmerCount7MT.countFastq(s, null, k, cbits, true, kca, trusted, extraMaxreads, thresh2, stepsize, true);
					}
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			kca.shutdown();
		}
		
		
		thash.stop();
//		Data.sysout.println(kca.toString());
		Data.sysout.println("Table:     \t"+kca.toShortString());
		Data.sysout.println("Build time:      \t"+thash);
		return kca;
	}
	
	public static long detect(String in1, String in2, KCountArray kca, int k, int thresh, long maxReads, String output, boolean ordered, boolean overwrite) {
		final ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
			Thread th=new Thread(cris);
			th.start();
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		RTextOutputStream3 ros=null;
		RTextOutputStream3 rosbad=null;
		if(output!=null){
			String out1=output.replaceFirst("#", "1"), out2=null;
			
			if(cris.paired()){
				if(output.contains("#")){
					out2=output.replaceFirst("#", "2");
				}else{
					System.err.println("Writing interleaved.");
				}
			}

			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, OUTPUT_INFO_ONLY ? ".info" : null, true, overwrite, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, OUTPUT_INFO_ONLY ? ".info" : null, true, overwrite, ordered);
			ros=new RTextOutputStream3(ff1, ff2, buff, null, true);
			
			assert(!ff1.sam()) : "Sam files need reference info for the header.";

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));
			
			if(DONT_OUTPUT_BAD_READS){
				out1=out1.replace('\\', '/');
				out2=out2==null ? null : out2.replace('\\', '/');
				

				String outbad1=null;
				if(out1!=null && out1.lastIndexOf('/')>0){
					int slash=out1.lastIndexOf('/');
					String a=out1.substring(0, slash+1);
					String b=out1.substring(slash+1);
					outbad1=a+(b.replaceFirst("\\.", "_BAD."));
//					assert(false) : "\n"+a+"\n"+b+"\n"+outbad1+"\n";
				}else{
					outbad1=out1.replaceFirst("\\.", "_BAD.");
				}
//				assert(false) : outbad1;

				String outbad2=null;
				if(out2!=null && out2.lastIndexOf('/')>0){
					int slash=out2.lastIndexOf('/');
					String a=out2.substring(0, slash+1);
					String b=out2.substring(slash+1);
					outbad2=a+(b.replaceFirst("\\.", "_BAD."));
//					assert(false) : "\n"+a+"\n"+b+"\n";
				}else if(out2!=null){
					outbad2=out2.replaceFirst("\\.", "_BAD.");
				}
				
				FileFormat ffb1=FileFormat.testOutput(outbad1, FileFormat.FASTQ, OUTPUT_INFO_ONLY ? ".info" : null, true, overwrite, ordered);
				FileFormat ffb2=FileFormat.testOutput(outbad2, FileFormat.FASTQ, OUTPUT_INFO_ONLY ? ".info" : null, true, overwrite, ordered);
				ros=new RTextOutputStream3(ffb1, ffb2, buff, null, true);
			}
		}
		
		
		if(ros!=null){
			ros.start();
			Data.sysout.println("Started output threads.");
		}
		if(rosbad!=null){
			rosbad.start();
		}
		
		long bases=detect(cris, kca, k, thresh, maxReads, ros, rosbad);
		
		ReadWrite.closeStreams(cris, ros, rosbad);
		if(verbose){System.err.println("Closed stream");}
		return bases;
	}
	
	public static long detect(ConcurrentReadStreamInterface cris, KCountArray kca, int k, int thresh, long maxReads, RTextOutputStream3 ros, RTextOutputStream3 rosbad) {
		Timer tdetect=new Timer();
		tdetect.start();

		long covered=0;
		long uncovered=0;
		
		long coveredFinal=0;
		long uncoveredFinal=0;
		
		long fullyCorrected=0;
		long failed=0;
		
		long errorsCorrected=0;
		
		long totalBases=0;
		long totalReads=0;
		long readsOut=0;
		long basesOut=0;
		long readsTrimmed=0;
		long basesTrimmed=0;
		
		ProcessThread[] pta=new ProcessThread[THREADS];
		for(int i=0; i<pta.length; i++){
			pta[i]=new ProcessThread(cris, kca, k, thresh, ros, rosbad);
			pta[i].start();
		}
		
		for(int i=0; i<pta.length; i++){
			ProcessThread ct=pta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				covered+=ct.covered;
				uncovered+=ct.uncovered;
				coveredFinal+=ct.coveredFinal;
				uncoveredFinal+=ct.uncoveredFinal;
				fullyCorrected+=ct.fullyCorrected;
				failed+=ct.failed;
				errorsCorrected+=ct.errorsCorrected;
				totalBases+=ct.totalBases;
				totalReads+=ct.totalReads;
				readsOut+=ct.readsOut;
				basesOut+=ct.basesOut;
				readsTrimmed+=ct.readsTrimmed;
				basesTrimmed+=ct.basesTrimmed;
			}
		}
		
//		ros.
		
		long perfect=totalReads-fullyCorrected-failed;
		
		Data.sysout.println();
		tdetect.stop();
		Data.sysout.println("Detect time:     \t"+tdetect+"   \t"+String.format("%.2f", totalBases*1000000.0/(tdetect.elapsed))+" kb/sec");
		Data.sysout.println("Total reads:     \t"+totalReads);
		Data.sysout.println("Total bases:     \t"+totalBases);
		Data.sysout.println("Reads Perfect:   \t"+perfect+"   \t"+String.format("%.2f%%", perfect*100.0/totalReads));
		Data.sysout.println("Reads Corrected: \t"+fullyCorrected+"   \t"+String.format("%.2f%%", fullyCorrected*100.0/totalReads));
		Data.sysout.println("Errors Corrected:\t"+errorsCorrected+"   \t"+String.format("%.3f avg", errorsCorrected*1.0/fullyCorrected));
		Data.sysout.println("Reads Failed:    \t"+failed+"   \t"+String.format("%.2f%%", failed*100.0/totalReads));
		if(TRIM_LEFT || TRIM_RIGHT){
			Data.sysout.println("Reads Trimmed:   \t"+readsTrimmed+"   \t"+String.format("%.2f%%", readsTrimmed*100.0/totalReads));
			Data.sysout.println("Bases Trimmed:   \t"+basesTrimmed+"   \t"+String.format("%.2f%%", basesTrimmed*100.0/totalBases));
		}
		Data.sysout.println("Reads Out:       \t"+readsOut+"   \t"+String.format("%.2f%%", readsOut*100.0/totalReads));
		Data.sysout.println("Bases Out:       \t"+basesOut+"   \t"+String.format("%.2f%%", basesOut*100.0/totalBases));
		
		Data.sysout.println("\n - bases before correction - ");
		Data.sysout.println("Covered:         \t"+covered);
		Data.sysout.println("Uncovered:       \t"+uncovered);
		
		Data.sysout.println("\n -  bases after correction - ");
		Data.sysout.println("Covered:         \t"+coveredFinal);
		Data.sysout.println("Uncovered:       \t"+uncoveredFinal);
		
		return totalBases;
	}
	
	
	public static BitSet detectErrors(Read r, KCountArray kca, int k, int thresh){
		if(kca.gap>0){return detectErrorsSplit(r, kca, k, thresh);}
		
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		final int gap=kca.gap;
		
		int bslen=r.bases.length-k-gap+1;
		if(bslen<1){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(bslen);
		
		int len=0;
		long kmer=0;
		byte[] bases=r.bases;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				
				if(len>=k){
					int count=kca.read(CANONICAL ? KCountArray.makeCanonical2(kmer, k) : kmer);
					if(count>=thresh){
						bs.set(i+1-k);
					}
				}
			}
		}
		return bs;
	}
	
	public static BitSet detectErrorsBulk(final Read r, final KCountArray kca, final int k, final int thresh, final int stepsize/*, final int offset*/){
		if(ONLY_CORRECT_N){return detectNBulk(r);}
		if(kca.gap>0){return detectErrorsSplit(r, kca, k, thresh);}
		
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		final int gap=kca.gap;
		
		if(r.bases==null || r.bases.length<k+gap){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.bases.length);
		final int setlen=k+gap;
		
		int len=0;
		long kmer=0;
		byte[] bases=r.bases;
		
//		final int sub=
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				
				if(len>=k && ((len-k)%stepsize==0 || i==bases.length-1)){
					int count=kca.read(CANONICAL ? KCountArray.makeCanonical2(kmer, k) : kmer);
					if(count>=thresh){
						bs.set(i+1-setlen, i+1);
					}
				}
			}
		}
		
		r.errors=r.bases.length-bs.cardinality();
		
//		assert(bases.length==r.bases.length);
		return bs;
	}
	
	public static BitSet detectNBulk(final Read r){
		if(r.bases==null){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.bases.length);
		
		for(int i=0; i<r.bases.length; i++){
			if(AminoAcid.isFullyDefined(r.bases[i])){bs.set(i);}
		}
		r.errors=r.bases.length-bs.cardinality();
		return bs;
	}
	
	public static BitSet detectTrusted(final Read r, final KCountArray kca, final int k, final int thresh, final int detectStepsize){
		if(kca.gap>0){throw new RuntimeException("TODO");}
		
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		final int gap=kca.gap;
		
		if(r.bases==null || r.bases.length<k+gap){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.bases.length);
		bs.set(0, r.bases.length);
		final int setlen=k+gap;
		
		int len=0;
		long kmer=0;
		byte[] bases=r.bases;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				
				if(len>=k && (i%detectStepsize==0 || i==bases.length-1)){
					int count=kca.read(CANONICAL ? KCountArray.makeCanonical2(kmer, k) : kmer);
					if(count<thresh){
						bs.clear(i+1-setlen, i+1);
//						bs.clear(i+1-setlen+detectStepsize, i+1-detectStepsize);
//						bs.clear(i+k/2-detectStepsize, i+k/2+detectStepsize);
//						bs.clear(i+k/2);
					}
				}
			}
		}
//		assert(bases.length==r.bases.length);
		return bs;
	}
	
	public static BitSet detectErrorsTips(final Read r, final KCountArray kca, final int k, final int thresh){
		if(kca.gap>0){return detectErrorsSplit(r, kca, k, thresh);}
		
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		final int gap=kca.gap;
		
		if(r.bases==null || r.bases.length<k+gap){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.bases.length);
		final int setlen=k+gap;
		
		int len=0;
		long kmer=0;
		byte[] bases=r.bases;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				
				if(len>=k){
					int count=kca.read(CANONICAL ? KCountArray.makeCanonical2(kmer, k) : kmer);
					if(count>=thresh){
						bs.set(i+1-setlen);
						bs.set(i);
					}
				}
			}
		}
		return bs;
	}

	/**
	 * @param r
	 * @param kca
	 * @param k
	 * @param thresh
	 * @return
	 */
	private static BitSet detectErrorsSplit(Read r, KCountArray kca, int k,
			int thresh) {
		assert(false) : "TODO";
		return null;
	}
	

	/** Assumes bulk mode was used; e.g., any '0' bit is covered by no correct kmers */
	public static BitSet correctErrors(final Read r, final KCountArray kca, final int k, final int thresh, BitSet bs, final int maxCorrections, final int maxBurst){
		assert(!TRY_BOTH_SIDES);
		if(kca.gap>0){assert(false) : "TODO";}
		assert(!OUTPUT_INFO_ONLY) : "TODO: Outputting correction data is not yet supported.";
		
		int corrections=0; //Alternately, corrections=r.errorsCorrected
		r.errors=0;

		int initialErrors=r.bases.length-bs.cardinality();
		if(initialErrors==0){return bs;}
		
		if(initialErrors>r.bases.length-k){//Cannot be corrected
			r.errors=r.bases.length;
			return bs;
		}
		
		
		byte[] bases0=Arrays.copyOf(r.bases, r.bases.length);
		byte[] qual0=(r.quality==null ? null : Arrays.copyOf(r.quality, r.quality.length));

//		verbose=!bs.get(0);
		if(verbose){
			Data.sysout.println();
			Data.sysout.println(toString(bs, r.bases.length));
			Data.sysout.println(toString(detectErrorsTips(r, kca, k, thresh), r.bases.length));
			Data.sysout.println(toString(detectErrors(r, kca, k, thresh), r.bases.length-k+1));
		}
		
		int lastloc=-99;
		int burst=1;
		while(!bs.get(0) && corrections<maxCorrections){//While the read starts with a '0', correct from the right.
			assert(initialErrors>0);
			assert(initialErrors>corrections);
//			Data.sysout.println("Could not correct.");
//			return bs;
			int errorLoc=bs.nextSetBit(0)-1;//Location to left of first '1'
			if(Tools.absdif(errorLoc,lastloc)<=BURST_THRESH){burst++;}
			else{burst=1;}
			lastloc=errorLoc;
			boolean success=(burst<=maxBurst) && correctFromRight(r, kca, k, thresh, bs, errorLoc);
			if(success){
				corrections++;
				bs=detectErrorsBulk(r, kca, k, thresh, 1);
				if(verbose){System.err.println(">\n"+toString(bs, r.bases.length));}
			}else{
				if(verbose){System.err.println("Could not correct.");}
				
				r.bases=bases0;
				r.quality=qual0;
				r.errors=initialErrors;
				
//				r.errors=r.bases.length-bs.cardinality();
//				r.errorsCorrected+=corrections;
//				r.bases[errorLoc]='N';
//				r.quality[errorLoc]=0;
				
				return bs;
			}
		}
		
		burst=1;
		while(bs.cardinality()<r.bases.length && corrections<maxCorrections){
			assert(initialErrors>0);
			assert(initialErrors>corrections);
			if(bs.get(0)){//First bit is a "1", can correct from the left
				int errorLoc=bs.nextClearBit(0);//Location to left of first '0'
				if(Tools.absdif(errorLoc,lastloc)<=BURST_THRESH){burst++;}
				else{burst=1;}
				lastloc=errorLoc;
				boolean success=(burst<=maxBurst) && correctFromLeft(r, kca, k, thresh, bs, errorLoc);
				if(success){
					corrections++;
					bs=detectErrorsBulk(r, kca, k, thresh, 1);
					if(verbose){System.err.println(">\n"+toString(bs, r.bases.length));}
				}else{
					if(verbose){System.err.println("Could not correct.");}
					
					r.bases=bases0;
					r.quality=qual0;
					r.errors=initialErrors;
					
//					r.errors=r.bases.length-bs.cardinality();
//					r.errorsCorrected+=corrections;
//					r.bases[errorLoc]='N';
//					r.quality[errorLoc]=0;
					
					return bs;
				}
			}
		}
		
		if(corrections>=maxCorrections && bs.cardinality()<r.bases.length){
			r.bases=bases0;
			r.quality=qual0;
			r.errors=initialErrors;
			
//			r.errors=r.bases.length-bs.cardinality();
//			r.errorsCorrected+=corrections;
//			r.bases[errorLoc]='N';
//			r.quality[errorLoc]=0;
			
			return bs;
		}
		
		r.errors=r.bases.length-bs.cardinality();
		r.errorsCorrected+=corrections;
		assert(corrections<=initialErrors);
		assert(corrections<=maxCorrections);
		assert(corrections>0 || maxCorrections<1) : "\ncorrections="+corrections+", maxCorrections="+maxCorrections+",\n" +
			"r.bases.length="+r.bases.length+", initialErrors="+initialErrors+", r.errors="+r.errors;
		return bs;
	}
	
	
	
	
	/** Assumes bulk mode was used; e.g., any '0' bit is covered by no correct kmers.
	 * This function  */
	public static BitSet correctErrorsBothSides(final Read r, final KCountArray kca, final int k, final int goodThresh, final int badThresh, BitSet bs, final int maxCorrections){
		
//		verbose=r.numericID==405093;
		
		assert(goodThresh>badThresh) : goodThresh+", "+badThresh;
		assert(goodThresh<=kca.maxValue) : goodThresh+", "+kca.maxValue+", "+kca.cellBits;
		
//		assert(false) : "TODO";
		if(kca.gap>0){assert(false) : "TODO";}
		assert(!OUTPUT_INFO_ONLY) : "TODO: Outputting correction data is not yet supported.";
		
		int corrections=0; //Alternately, corrections=r.errorsCorrected
		r.errors=0;
		
		byte[] bclone=r.bases.clone();
		
		int initialErrors=r.bases.length-bs.cardinality();
		if(initialErrors>r.bases.length-k){//Cannot be corrected
			r.errors=r.bases.length;
			return bs;
		}
		if(initialErrors==0){return bs;} //Nothing to correct.
		

		int prevBlock=-1;
		int blockStart=bs.nextClearBit(0);
		int blockStop=bs.nextSetBit(blockStart)-1;
		if(blockStop<0){blockStop=r.bases.length-1;}

//		verbose=!bs.get(0);
		if(verbose){
			Data.sysout.println();
			Data.sysout.println(new String(r.bases));
			Data.sysout.println(toString(bs, r.bases.length));
			Data.sysout.println(toString(detectErrorsTips(r, kca, k, goodThresh), r.bases.length));
			Data.sysout.println(toString(detectErrors(r, kca, k, goodThresh), r.bases.length-k+1));
			Data.sysout.println("prevBlock="+prevBlock+", blockStart="+blockStart+", blockStop="+blockStop);
			Data.sysout.println();
		}
		
		while(corrections<initialErrors && blockStart<r.bases.length){
			
			if(verbose){
				Data.sysout.println("Entering loop.");
				Data.sysout.println(new String(r.bases));
				Data.sysout.println(toString(bs, r.bases.length));
				Data.sysout.println("prevBlock="+prevBlock+", blockStart="+blockStart+", blockStop="+blockStop);
			}
			
			assert(blockStop>=0) : "\n"+blockStart+"\n"+blockStop+"\n"+prevBlock+"\n"+toString(bs, r.bases.length)+"\n"+r+"\n";
			int nextBlock=bs.nextClearBit(blockStop+1);
			if(verbose){
				Data.sysout.println("nextBlock="+nextBlock);
			}
			
			assert(prevBlock<blockStart);
			assert(blockStart<=blockStop);
			assert(blockStop<nextBlock);
			int leftSpace=blockStart-prevBlock-1;
			
			int x=0;
			if(leftSpace+1>=k+kca.gap){//correct from left
				x=correctFullyFromLeft(r, kca, k, goodThresh, badThresh, bs, blockStart);
				if(verbose){System.err.println("Left: "+x);}
			}
			if(!bs.get(blockStop)){
				int rightSpace=nextBlock-blockStop-1;
				if(rightSpace+1>=k+kca.gap){
					x+=correctFullyFromRight(r, kca, k, goodThresh, badThresh, bs, blockStop);
					if(verbose){System.err.println("Right: "+x);}
				}
			}
			
			corrections+=x;
			if(corrections>=blockStop-blockStart+1){//then the whole block was cleared
				
			}else{
				prevBlock=bs.nextClearBit(blockStart);
				prevBlock=bs.nextSetBit(prevBlock)-1;
			}
			
			blockStart=nextBlock;
			blockStop=bs.nextSetBit(blockStart)-1;
			if(blockStop<0){blockStop=r.bases.length-1;}
		}
		
		
		r.errors=0;//r.bases.length-bs.cardinality();
		r.errorsCorrected+=corrections;
		assert(corrections<=initialErrors) : corrections+", "+initialErrors+"\n"+new String(r.bases)+"\n"+new String(bclone)+"\n";
//		assert(corrections<=maxCorrections);
//		assert(corrections>0 || maxCorrections<1);
		return bs;
	}
	
	
	
	/**
	 * @param r
	 * @param kca
	 * @param k
	 * @param thresh
	 * @param bs
	 * @param errorLoc
	 * @return
	 */
	private static boolean correctFromLeft(Read r, KCountArray kca, int k, int thresh, BitSet bs, int error) {
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		final int gap=kca.gap;
		final int setlen=k+gap;
		final int startLoc=error-(setlen)+1;
		final byte oldBase=r.bases[error];
		final byte[] bases=r.bases;
		
		final int minAdvance=Tools.min(MIN_ADVANCE, bases.length-error);
		
		long kmer=0;
		int len=0;
		for(int i=startLoc; i<error; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				throw new RuntimeException("Can't correct from left!\nerror="+error+"\n"+toString(bs, bases.length)+"\n"+new String(bases)+"\nreadID: "+r.numericID);
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
			}
		}
		assert(len==setlen-1) : setlen+", "+len+", "+error+", "+startLoc;
		
		int[] counts=new int[4];
		int[] dists=new int[4];
		int maxLoc=Tools.min(bases.length-1, error+setlen-1);
		if(!bs.get(error+1)){maxLoc=Tools.min(maxLoc, error+9);}
		else{
			for(int i=error+2; i<=maxLoc; i++){
				if(!bs.get(i)){
					maxLoc=i-1;
					break;
				}
			}
		}
		
		if(verbose){System.err.println("correctFromLeft.  Error = "+error+", maxloc="+maxLoc);}
		for(int bnum=0; bnum<4; bnum++){
			byte c=AminoAcid.numberToBase[bnum];
			bases[error]=c;
			if(verbose){System.err.println("Considering "+(char)c);}
			long key=kmer;
			for(int loc=error; loc<=maxLoc; loc++){
				c=bases[loc];
				int x=AminoAcid.baseToNumber[c];
				if(x<0){
					if(verbose){System.err.println("break: N");}
					break;
				}
				key=((key<<2)|x)&mask;
				int count=kca.read(CANONICAL ? KCountArray.makeCanonical2(key, k) : key);
				if(count<thresh){
					if(verbose){System.err.println("break: count="+count);}
					break;
				}
				dists[bnum]++;
				counts[bnum]+=count;
			}
		}
		bases[error]=oldBase;
		
		//Note:  I could require both to be the same, to decrease false-positives
		
		final int muid=maxUniqueIndex(dists);
		Arrays.sort(dists);
		final int advance=dists[3];
		final int delta=dists[3]-dists[2];
//		if(advance<minAdvance){return false;}
		if(delta<minAdvance){return false;}
		
		int best=(muid<0 ? maxUniqueIndex(counts) : muid);
		
		if(verbose){System.err.println("Best="+best+": "+Arrays.toString(dists)+"  \t"+Arrays.toString(counts));}
		if(best<0){return false;}
		byte bestC=AminoAcid.numberToBase[best];
		if(bestC==oldBase){return false;}
		bases[error]=bestC;
		
		if(r.quality!=null){r.quality[error]=(byte)Tools.min(25, 18+delta);}
		
		return true;
	}
	
	
	
	/**
	 * @param r
	 * @param kca
	 * @param k
	 * @param thresh
	 * @param bs
	 * @param errorLoc
	 * @return
	 */
	private static boolean correctFromRight(Read r, KCountArray kca, int k, int thresh, BitSet bs, int error) {
		final int kbits=2*k;
		final int shift=kbits-2;
		final long mask=~((-1L)<<(kbits));
		final int gap=kca.gap;
		final int setlen=k+gap;
		final int stopLoc=error+(setlen)-1;
		final byte oldBase=r.bases[error];
		final byte[] bases=r.bases;
		
		final int minAdvance=Tools.min(MIN_ADVANCE, error+1);
		
		long kmer=0;
		int len=0;
		
		for(int i=error+1; i<=stopLoc; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				throw new RuntimeException("Can't correct from right!\nerror="+error+"\n"+toString(bs, bases.length)+"\n"+new String(bases));
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
			}
//			Data.sysout.print((char)b);
		}
		kmer<<=2;
		
		if(verbose){
			Data.sysout.println();
			String s=Long.toBinaryString(kmer);
			while(s.length()<kbits){s="0"+s;}
			Data.sysout.println("kmer = \t"+s);
		}
		assert(len==setlen-1) : setlen+", "+len+", "+error+", "+stopLoc;
		
		int[] counts=new int[4];
		int[] dists=new int[4];
		int minLoc=Tools.max(0, error-setlen+1);
		if(error==0 || !bs.get(error-1)){minLoc=Tools.max(minLoc, error-9);}
		else{
			for(int i=error-2; i>=minLoc; i--){
				if(!bs.get(i)){
					minLoc=i+1;
					break;
				}
			}
		}
		
		if(verbose){
			Data.sysout.println("correctFromRight.  Error = "+error+", minloc="+minLoc);
			Data.sysout.println(new String(r.bases));
		}
		for(int bnum=0; bnum<4; bnum++){
			byte c=AminoAcid.numberToBase[bnum];
			bases[error]=c;
			if(verbose){System.err.println("Considering "+(char)c);}
			long key=kmer;
			for(int loc=error; loc>=minLoc; loc--){
				c=bases[loc];
				int x=AminoAcid.baseToNumber[c];
				if(x<0){
					if(verbose){System.err.println("break: N");}
					break;
				}
				key=((key>>2)|(((long)x)<<shift))&mask;
//				{
//					String s=Long.toBinaryString(key);
//					while(s.length()<kbits){s="0"+s;}
//					Data.sysout.println("mask="+Long.toBinaryString(mask)+", shift="+shift+", c="+c+", x="+x+", key  = \t"+s);
//				}
				int count=kca.read(CANONICAL ? KCountArray.makeCanonical2(key, k) : key);
				if(count<thresh){
					if(verbose){System.err.println("break: count="+count);}
					break;
				}
				dists[bnum]++;
				counts[bnum]+=count;
			}
		}
		bases[error]=oldBase;
		
		//Note:  I could require both to be the same, to decrease false-positives
		
		final int muid=maxUniqueIndex(dists);
		Arrays.sort(dists);
		final int advance=dists[3];
		final int delta=dists[3]-dists[2];
//		if(advance<minAdvance){return false;}
		if(delta<minAdvance){return false;}
		
		int best=(muid<0 ? maxUniqueIndex(counts) : muid);
		
		if(verbose){System.err.println("Best="+best+": "+Arrays.toString(dists)+"  \t"+Arrays.toString(counts));}
		if(best<0){return false;}
		byte bestC=AminoAcid.numberToBase[best];
		if(bestC==oldBase){return false;}
		bases[error]=bestC;
		
		if(r.quality!=null){r.quality[error]=(byte)Tools.min(25, 18+delta);}
		
		return true;
	}
	
	
	private static int correctFullyFromRight(Read r, KCountArray kca, int k, int threshGood, int threshBad, BitSet bs, int error) {
//		assert(false) : "TODO";
		final int kbits=2*k;
		final int shift=kbits-2;
		final long mask=~((-1L)<<(kbits));
		final int gap=kca.gap;
		final int setlen=k+gap;
		final int stopLoc=error+(setlen)-1;
		final byte[] bases=r.bases;
		
		long kmer=0;
		int len=0;
		
		for(int i=error+1; i<=stopLoc; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				throw new RuntimeException("Can't correct from right!\nerror="+error+"\n"+toString(bs, bases.length)+"\n"+new String(bases));
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
			}
//			Data.sysout.print((char)b);
		}
		kmer<<=2;
		
		if(verbose){
			Data.sysout.println();
			String s=Long.toBinaryString(kmer);
			while(s.length()<kbits){s="0"+s;}
			Data.sysout.println("kmer = \t"+s);
		}
		assert(len==setlen-1) : setlen+", "+len+", "+error+", "+stopLoc;
		
		int[] counts=new int[4];
		int minLoc=error;
		while(minLoc>0 && !bs.get(minLoc)){minLoc--;}
		if(bs.get(minLoc)){minLoc++;}
		
		if(verbose){
			Data.sysout.println("correctFullyFromRight.  Error = "+error+", minloc="+minLoc);
			Data.sysout.println(new String(bases));
		}
		
		boolean success=true;
		int corrected=0;
		for(int loc=error; loc>=minLoc; loc--){
			int lows=0, highs=0, zeros=0;
			for(int bnum=0; bnum<4; bnum++){
				if(verbose){System.err.println("Considering "+(char)AminoAcid.numberToBase[bnum]);}
				long key=kmer;
				
				if(bnum<0){
					if(verbose){System.err.println("break: N");}
					break;
				}
				key=((key>>2)|(((long)bnum)<<shift))&mask;
//				{
//					String s=Long.toBinaryString(key);
//					while(s.length()<kbits){s="0"+s;}
//					Data.sysout.println("mask="+Long.toBinaryString(mask)+", shift="+shift+", c="+c+", x="+x+", key  = \t"+s);
//				}
				int count=kca.read(CANONICAL ? KCountArray.makeCanonical2(key, k) : key);
				if(count<=threshBad){lows++;}
				else if(count>=threshGood){highs++;}
				if(count==0){zeros++;}
				counts[bnum]=count;
			}
			assert(zeros<=lows);
			
			if((highs==1 && lows==3) || (bases[loc]=='N' && zeros==3)){
				int x=Tools.maxIndex(counts);
				if(verbose){System.err.println("Best="+x+":  \t"+Arrays.toString(counts));}
				bases[loc]=AminoAcid.numberToBase[x];
				if(r.quality!=null){r.quality[loc]=(byte)(20+zeros+(highs==1 ? 3 : 0));}
				corrected++;
				kmer=((kmer>>2)|(((long)x)<<shift))&mask;
				bs.set(loc);
			}else{
				success=false;
				break;
			}
		}
		
		return corrected;
	}
	
	
	
	private static int correctFullyFromLeft(Read r, KCountArray kca, int k, int threshGood, int threshBad, BitSet bs, int error) {
//		assert(false);
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		final int gap=kca.gap;
		final int setlen=k+gap;
		final int startLoc=error-(setlen)+1;
		final byte[] bases=r.bases;
		
		long kmer=0;
		int len=0;
		for(int i=startLoc; i<error; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				throw new RuntimeException("Can't correct from left!\nerror="+error+"\n"+toString(bs, bases.length)+"\n"+new String(bases)+"\nreadID: "+r.numericID);
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
			}
		}
		assert(len==setlen-1) : setlen+", "+len+", "+error+", "+startLoc;
		
		int[] counts=new int[4];
		int maxLoc=error;
		while(maxLoc<bases.length && !bs.get(maxLoc)){maxLoc++;}
		if(maxLoc>=bases.length){maxLoc=bases.length-1;}
		if(bs.get(maxLoc)){maxLoc--;}
		
		if(verbose){
			Data.sysout.println("correctFullyFromLeft.  Error = "+error+", minloc="+maxLoc);
			Data.sysout.println(new String(bases));
		}
		
		boolean success=true;
		int corrected=0;
		for(int loc=error; loc<=maxLoc; loc++){

//			Data.sysout.println(bs);
//			Data.sysout.println(loc);
//			Data.sysout.println(bases[loc]);
//			Data.sysout.println();
			
			int lows=0, highs=0, zeros=0;
			for(int bnum=0; bnum<4; bnum++){
				if(verbose){System.err.println("Considering "+(char)AminoAcid.numberToBase[bnum]);}
				long key=kmer;
				
				if(bnum<0){
					if(verbose){System.err.println("break: N");}
					break;
				}
				key=((key<<2)|bnum)&mask;
//				{
//					String s=Long.toBinaryString(key);
//					while(s.length()<kbits){s="0"+s;}
//					Data.sysout.println("mask="+Long.toBinaryString(mask)+", shift="+shift+", c="+c+", x="+x+", key  = \t"+s);
//				}
				int count=kca.read(CANONICAL ? KCountArray.makeCanonical2(key, k) : key);
				if(count<=threshBad){lows++;}
				else if(count>=threshGood){highs++;}
				if(count==0){zeros++;}
				counts[bnum]=count;
			}
			assert(zeros<=lows);
//			assert(zeros>0);
			
			if((highs==1 && lows==3) || (bases[loc]=='N' && zeros==3)){
				assert(zeros<4);
				int x=Tools.maxIndex(counts);
				if(verbose){System.err.println("Best="+x+":  \t"+Arrays.toString(counts));}
				bases[loc]=AminoAcid.numberToBase[x];
				if(r.quality!=null){r.quality[loc]=(byte)(20+zeros+(highs==1 ? 3 : 0));}
				corrected++;
				kmer=((kmer<<2)|x)&mask;
				bs.set(loc);
//				assert(highs==1) : "\nloc="+loc+", zeros="+zeros+", lows="+lows+", highs="+highs+
//					", counts="+Arrays.toString(counts)+", bases["+loc+"]="+((char)bases[loc])+
//					"\n"+new String(bases)+"\n"+toString(bs, bases.length)+"\n";
			}else{
//				assert(zeros<3 || bases[bases.length-1]=='N') : "\nloc="+loc+", zeros="+zeros+", lows="+lows+", highs="+highs+
//					", counts="+Arrays.toString(counts)+", bases["+loc+"]="+((char)bases[loc])+
//					"\n"+new String(bases)+"\n"+toString(bs, bases.length)+"\n";
//				assert(false) : threshGood+", "+threshBad;
				success=false;
				break;
			}
		}
		
		return corrected;
	}

	/** returns index of highest value, if unique; else a negative number */
	private static int maxUniqueIndex(int[] array){
		int max=array[0];
		int maxIndex=0;
		for(int i=1; i<array.length; i++){
			if(array[i]>max){
				max=array[i];
				maxIndex=i;
			}else if(max==array[i]){
				maxIndex=-1;
			}
		}
		return maxIndex;
	}

	public static final String toString(BitSet bs, int len){
//		assert(verbose);
		StringBuilder sb=new StringBuilder(len);
		for(int i=0; i<len; i++){sb.append(bs.get(i) ? '1' : '0');}
		return sb.toString();
	}
	
	private static ArrayList<Read> removeBad(ArrayList<Read> list){
		
		ArrayList<Read> bad=new ArrayList<Read>();
		
		if(DONT_OUTPUT_BAD_PAIRS){
			for(int i=0; i<list.size(); i++){
				Read r=list.get(i);
				if((r.errors>0 || r.discarded()) || (r.mate!=null && (r.mate.errors>0 || r.mate.discarded()))){
					list.set(i, null);
					bad.add(r);
				}
			}
		}else{
			for(int i=0; i<list.size(); i++){
				Read r=list.get(i);
				if((r.errors>0 || r.discarded()) && (r.mate==null || (r.mate.errors>0 || r.mate.discarded()))){
					list.set(i, null);
					bad.add(r);
				}
			}
		}
		
		return bad;
	}
	
	private static void trim(Read r, BitSet bs, byte minq, int maxTrim) {
		if(bs==null){trim(r, minq, maxTrim);}
		else{
			assert(false) : "TODO";
			trim(r, minq, maxTrim);
		}
	}
	
	private static int trim(Read r, byte minq, int maxTrim) {
//		assert(r.bases.length>=MIN_LEN) : r.bases.length;
		assert(maxTrim>0) : maxTrim;
		if(maxTrim<1){return 0;}

//		Data.sysout.println(Arrays.toString(r.quality));
//		Data.sysout.println("TRIM_LEFT="+TRIM_LEFT+", TRIM_RIGHT="+TRIM_RIGHT+", minq="+minq+", TRIM_N="+TRIM_N+", TRIM_N_ONLY="+TRIM_N_ONLY);
		
		byte[] bases=r.bases;
		byte[] quals=r.quality;
		if(bases.length<MIN_LEN_2){
//			r.bases=null;
//			r.quality=null;
//			r.match=null;
//			assert(false);
			r.setDiscarded(true);
			return 0;
		}
		
		int left=0;
		int right=bases.length-1;
		int safe=0;
		final int minsafe=8;
		
//		assert(r.bases.length>=MIN_LEN) : r.bases.length;
		
		if(TRIM_LEFT){
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				byte q=quals[i];
				if((b=='N' && TRIM_N) || (q<minq && !TRIM_N_ONLY)){
					safe=0;
					left=i+1;
				}else{
					safe++;
					if(safe>=minsafe){break;}
				}
			}
		}
		while(left>maxTrim){left--;}
		
		safe=0;

		if(TRIM_RIGHT && left<maxTrim){
			for(int i=bases.length-1; i>0; i--){
				byte b=bases[i];
				byte q=quals[i];
//				Data.sysout.println("Processing "+(char)b+" of q="+q);
//				Data.sysout.println("safe="+safe+", minsafe="+minsafe+", right="+right);
				if((b=='N' && TRIM_N) || (q<minq && !TRIM_N_ONLY)){
					safe=0;
					right=i-1;
				}else{
					safe++;
					if(safe>=minsafe){break;}
				}
			}
		}
		
		while(left+bases.length-right-1>maxTrim){right++;}
		
		if(left==0 && right==bases.length-1){
			//do nothing
			return 0;
		}else if(right-left+1<MIN_LEN_2){
//			r.match=null;
//			r.bases=null;
//			r.quality=null;
//			r.setMapped(false);
			r.setDiscarded(true);
			while(left>0 && right-left+1<MIN_LEN_2){left--;}
			while(right<bases.length-1 && right-left+1<MIN_LEN_2){right++;}
//			r.clearSite();
//			return;
		}else if(right-left+1<MIN_LEN){
			r.setDiscarded(true);
		}
		r.bases=Arrays.copyOfRange(bases, left, right+1);
		r.quality=Arrays.copyOfRange(quals, left, right+1);
		r.mapLength=r.bases.length;
		assert(r.bases.length>=MIN_LEN || r.discarded());
		return bases.length-r.bases.length;
	}
	
	
	private static class ProcessThread extends Thread{
		
		ProcessThread(ConcurrentReadStreamInterface cris_, KCountArray kca_, int k_, int thresh_, RTextOutputStream3 ros_, RTextOutputStream3 rosbad_){
			cris=cris_;
			kca=kca_;
			k=k_;
			thresh=thresh_;
			ros=ros_;
			rosbad=rosbad_;
		}
		
		public void run(){
			detect();
		}
		
		void detect() {
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			
			while(reads!=null && reads.size()>0){
				for(Read r : reads){
					Read r2=r.mate;
//					assert(ec==0); //*** for testing
//					
//					boolean cor=false;
					{
						
//						if(r.numericID==23){verbose=true;}

						int initialErrorsCorrected=r.errorsCorrected;
						int initialLen=r.bases.length;
						totalReads++;
						if(verbose){System.err.println();}
						totalBases+=r.bases.length;
//						BitSet bs=detectErrors(r, kca, k, thresh);
						BitSet bs=detectErrorsBulk(r, kca, k, thresh, 1);
						if(verbose){System.err.println(ErrorCorrectMT.toString(bs, r.bases.length));}
//						Data.sysout.println(toString(detectErrorsTips(r, kca, k, thresh), r.bases.length));
						if(verbose){System.err.println(ErrorCorrectMT.toString(detectErrors(r, kca, k, thresh), r.bases.length-k+1));}
						if(bs==null){//can't detect errors
//							assert(false);
							r.setDiscarded(true);
						}else{
							int initialCorrect=bs.cardinality();
							covered+=initialCorrect;
							uncovered+=(r.bases.length-initialCorrect);
							if(initialCorrect<r.bases.length){
								if(TRY_BOTH_SIDES){bs=correctErrorsBothSides(r, kca, k, thresh, /*Tools.min(1, thresh-1)*/THRESH_BAD, bs, ERROR_CORRECTION_LIMIT);}
								else{bs=correctErrors(r, kca, k, thresh, bs, ERROR_CORRECTION_LIMIT, MAX_ERROR_BURST, (TRIM_LEFT || TRIM_RIGHT));}
							}
							int correct=bs.cardinality();
							coveredFinal+=correct;
							uncoveredFinal+=(r.bases.length-correct);
							
							int errorsNewlyCorrected=r.errorsCorrected-initialErrorsCorrected;
							errorsCorrected+=errorsNewlyCorrected;
							
							if(initialCorrect<initialLen){
								if(correct==r.bases.length){
									fullyCorrected++;
									assert(errorsNewlyCorrected>0 || initialLen>r.bases.length);
								}else{
									failed++;
									assert(errorsNewlyCorrected==0 || TRY_BOTH_SIDES);
								}
							}else{
								assert(errorsNewlyCorrected==0 || initialLen>r.bases.length);
							}
						}
					}
					if(r2!=null){
//						assert(false); //***
						int initialErrorsCorrected=r2.errorsCorrected;
						int initialLen=r2.bases.length;
						totalReads++;
						totalBases+=r2.bases.length;
//						BitSet bs=detectErrors(r2, kca, k, thresh);
						BitSet bs=detectErrorsBulk(r2, kca, k, thresh, 1);
						if(verbose){System.err.println(ErrorCorrectMT.toString(bs, r2.bases.length));}
//						Data.sysout.println(toString(detectErrorsTips(r2, kca, k, thresh), r2.bases.length));
						if(verbose){System.err.println(ErrorCorrectMT.toString(detectErrors(r2, kca, k, thresh), r2.bases.length-k+1));}
						if(bs==null){//can't detect errors
							r.setDiscarded(true);
						}else{
							int initialCorrect=bs.cardinality();
							covered+=initialCorrect;
							uncovered+=(r2.bases.length-initialCorrect);
							if(initialCorrect<r2.bases.length){
								if(TRY_BOTH_SIDES){bs=correctErrorsBothSides(r2, kca, k, thresh, /*Tools.min(1, thresh-1)*/THRESH_BAD, bs, ERROR_CORRECTION_LIMIT);}
								else{bs=correctErrors(r2, kca, k, thresh, bs, ERROR_CORRECTION_LIMIT, MAX_ERROR_BURST, (TRIM_LEFT || TRIM_RIGHT));}
							}
							int correct=bs.cardinality();
							coveredFinal+=correct;
							uncoveredFinal+=(r2.bases.length-correct);
							
							int errorsNewlyCorrected=r2.errorsCorrected-initialErrorsCorrected;
							errorsCorrected+=errorsNewlyCorrected;
							
							if(initialCorrect<r2.bases.length){
								if(correct==r2.bases.length){
									fullyCorrected++;
									assert(errorsNewlyCorrected>0 || initialLen>r2.bases.length);
								}else{
									failed++;
									assert(errorsNewlyCorrected==0);
								}
							}else{
								assert(errorsNewlyCorrected==0 || initialLen>r2.bases.length);
							}
						}
					}
				}
				

				if(DONT_OUTPUT_BAD_READS){
					ArrayList<Read> bad=removeBad(reads);
					if(rosbad!=null){
						rosbad.add(bad, ln.id);
					}
				}
				for(Read r : reads){
					if(r!=null){
						Read r2=r.mate;
						readsOut++;
						basesOut+=(r.bases==null ? 0 : r.bases.length);
						r.obj=null;
						assert(r.bases!=null);
						if(r.sites!=null && r.sites.isEmpty()){r.sites=null;}
						
						if(r2!=null){
							readsOut++;
							basesOut+=(r2.bases==null ? 0 : r2.bases.length);
							r2.obj=null;
							assert(r2.bases!=null);
							if(r2.numSites()==0){r2.sites=null;}
						}
					}
				}
				//System.err.println("Adding list of length "+readlist.size());
				if(ros!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
					ros.add(reads, ln.id);
				}
				
				
				cris.returnList(ln, ln.list.isEmpty());
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln, ln.list.isEmpty());
			if(verbose){System.err.println("Returned list");}
		}
		

		/** Assumes bulk mode was used; e.g., any '0' bit is covered by no correct kmers */
		private BitSet correctErrors(final Read r, final KCountArray kca, final int k, final int thresh, BitSet bs, final int maxCorrections, final int maxBurst, boolean trim){
			if(r.discarded()){r.setDiscarded(false);}
			
			if(!trim){return ErrorCorrectMT.correctErrors(r, kca, k, thresh, bs, maxCorrections, maxBurst);}
			
			byte q=TRIM_QUAL;
			
			byte[] qual=r.quality;
			byte[] bases=r.bases;
			int initialErrors=r.errors;
			int initialCorrected=r.errorsCorrected;
//			Data.sysout.println("A");
			assert(initialErrors>=0) : initialErrors;
//			Data.sysout.println("initialErrors = "+initialErrors+", errors = "+r.errors+", corrected = "+r.errorsCorrected+", discarded = "+r.discarded());
			bs=ErrorCorrectMT.correctErrors(r, kca, k, thresh, bs, maxCorrections, maxBurst);
			
			int maxtrim=Tools.min(MAX_TRIM_BASES, (r.bases.length*1)/4);
			for(int i=0; r.errors>0 && i<20 && q<=MAX_TRIM_QUAL && !r.discarded() && maxtrim>0; i++){
//				Data.sysout.println("B");
//				Data.sysout.println("Errors = "+r.errors+", corrected = "+r.errorsCorrected+", discarded = "+r.discarded());
				int x=trim(r, q, maxtrim);
				maxtrim-=x;
//				Data.sysout.println("trimmed: "+x);
//				Data.sysout.println("Errors = "+r.errors+", corrected = "+r.errorsCorrected+", discarded = "+r.discarded());
				if(x>0){
//					Data.sysout.println("C");
					bs=detectErrorsBulk(r, kca, k, thresh, 1);
//					Data.sysout.println("Errors = "+r.errors+", corrected = "+r.errorsCorrected+", discarded = "+r.discarded());
					bs=ErrorCorrectMT.correctErrors(r, kca, k, thresh, bs, maxCorrections, maxBurst);
//					Data.sysout.println("Errors = "+r.errors+", corrected = "+r.errorsCorrected+", discarded = "+r.discarded());
				}
				q=(byte) (q+2);
			}
//			Data.sysout.println("D");
//			Data.sysout.println("Errors = "+r.errors+", corrected = "+r.errorsCorrected+", discarded = "+r.discarded());
			
			if(r.errors>0 || r.discarded()){
//				Data.sysout.println("E");
				r.bases=bases;
				r.quality=qual;
				r.errors=initialErrors;
				r.errorsCorrected=initialCorrected;
			}else{
//				Data.sysout.println("F");
				int x=bases.length-r.bases.length;
				if(x>0){
//					Data.sysout.println("G");
					readsTrimmed++;
					basesTrimmed+=x;
					assert(x<=Tools.min(MAX_TRIM_BASES, (bases.length*3)/8)) : x+", "+bases.length+", "+r.bases.length+", "+Tools.min(MAX_TRIM_BASES, (bases.length*3)/8);
//					assert(r.bases.length>38);
				}
			}
//			Data.sysout.println("H");
//			Data.sysout.println("Errors = "+r.errors+", corrected = "+r.errorsCorrected+", discarded = "+r.discarded());
			return bs;
		}
		
		private final ConcurrentReadStreamInterface cris;
		private final KCountArray kca; 
		private final int k; 
		private final int thresh;
		private final RTextOutputStream3 ros;
		private final RTextOutputStream3 rosbad;

		long covered=0;
		long uncovered=0;

		long coveredFinal=0;
		long uncoveredFinal=0;

		long fullyCorrected=0;
		long readsTrimmed=0;
		long basesTrimmed=0;
		long failed=0;
		long readsOut=0;
		long basesOut=0;

		long totalBases=0;
		long totalReads=0;
		
		long errorsCorrected=0;
	}
	
	public static boolean verbose=false;
	/** Bails out if a read still has errors after correcting this many. */
	public static int ERROR_CORRECTION_LIMIT=4;
	/** Max allowed number of nearby corrections.  
	 * A long error burst indicates the read simply has low coverage, and is not being corrected correctly. */
	public static int MAX_ERROR_BURST=2;
	/** Bursts have at most this distance between errors. E.G. '1' means errors are adjacent. */
	public static int BURST_THRESH=2;	
	/** Max frequency of kmer for an alternate (not chosen) base in TRY_BOTH_SIDES mode.*/
	public static int THRESH_BAD=0;
	/** Withhold uncorrectable reads from output. */
	public static boolean DONT_OUTPUT_BAD_READS=false;
	/** Withhold uncorrectable reads from output, if either the read or its mate is bad. */
	public static boolean DONT_OUTPUT_BAD_PAIRS=false;
	/** Do not correct an error if it is at most this far from the next error.  Instead, bail out. */
	public static int MIN_ADVANCE=1;

//	/** Trim bases of this quality or below when trimming the left side of a read */
//	public static byte TRIM_QUAL_LEFT=6;
//	/** Trim bases of this quality or below when trimming the left side of a read */
//	public static byte TRIM_QUAL_RIGHT=6;
	
	/** Trim bases of this quality or below when trimming a read */
	public static byte TRIM_QUAL=5;
	/** Trim bases of up to this quality if a read still has errors after trimming at a lower quality */
	public static byte MAX_TRIM_QUAL=15;

	public static boolean TRIM_LEFT=false;
	public static boolean TRIM_RIGHT=false;
	public static boolean TRIM_N=true;
	public static boolean TRIM_N_ONLY=false;
	public static boolean CANONICAL=false;
	
	/** Number of threads used for error correction.  Does not control number of threads for creating the hash table. 
	 * Additionally, up to 2 threads are used for reading and up to 2 for writing. */ 
	public static int THREADS=8;

	/** Output correction data instead of the corrected read. */
	public static boolean OUTPUT_INFO_ONLY=false;
	/** Only detect/correct N, not called bases.  Mainly for filling read middles. */
	public static boolean ONLY_CORRECT_N=false;
	/** When an uncorrectable error is encountered, don't bail out, but instead try from the other side. This mode is far faster.  */
	public static boolean TRY_BOTH_SIDES=false;
	
	public static int MIN_LEN=45;
	public static int MIN_LEN_2=35;
	
	public static int MAX_TRIM_BASES=60;
	
}
