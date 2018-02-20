package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import align2.ListNum;
import align2.QualityTools;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import dna.AminoAcid;
import dna.Data;
import dna.Gene;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Jan 13, 2014
 *
 */
public class CalcTrueQuality {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		ReadStats.COLLECT_QUALITY_STATS=true;
		CalcTrueQuality ctq=new CalcTrueQuality(args);
		ReadStats.OVERWRITE=ctq.overwrite;
		ctq.process();
		
		ctq.writeMatrices();
	}
	
	public static void printOptions(){
		assert(false) : "No help available.";
	}
	
	public CalcTrueQuality(String[] args){
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
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=true;
		SamLine.CONVERT_CIGAR_TO_MATCH=true;
		
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
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
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
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in=b.split(",");
			}else if(a.equals("out") || a.equals("output") || a.equals("q102") || a.equals("q102out")){
				q102out=b;
			}else if(a.equals("qbp") || a.equals("qbpout")){
				qbpout=b;
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("tuc") || a.equals("touppercase")){
				Read.TO_UPPER_CASE=Tools.parseBoolean(b);
			}else if(in==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				in=arg.split(",");
			}else if(q102out==null && i==1 && !arg.contains("=")){
				q102out=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		assert(FastaReadInputStream.settingsOK());
//		if(maxReads!=-1){ReadWrite.USE_GUNZIP=ReadWrite.USE_UNPIGZ=false;}
		
		if(in==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.THREADS>2){
//			if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(!Tools.testOutputFiles(overwrite, false, q102out, qbpout, q10out, q12out, qb012out, qb234out, qpout, qout, pout)){
			throw new RuntimeException("\n\nOVERWRITE="+overwrite+"; Can't write to output file "+q102out+"\n");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(){
		Timer t=new Timer();
		t.start();
		for(String s : in){
			process(s);
		}
		
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

		rpstring=(readsUsed<100000 ? ""+readsUsed : readsUsed<100000000 ? (readsUsed/1000)+"k" : (readsUsed/1000000)+"m");
		bpstring=(basesUsed<100000 ? ""+basesUsed : basesUsed<100000000 ? (basesUsed/1000)+"k" : (basesUsed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}

		outstream.println("Reads Used:    "+rpstring);
		outstream.println("Bases Used:    "+bpstring);
		
		if(errorState){
			throw new RuntimeException(this.getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	public void process(String fname){
		
		final ConcurrentReadStreamInterface cris;
		final Thread cristhread;
		{
			FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, true, false);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, false, ff, null);
			if(verbose){System.err.println("Starting cris");}
			cristhread=new Thread(cris);
			cristhread.start();
		}
		
		{	
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					Read r1=reads.get(idx);
					Read r2=r1.mate;
					process(r1);
					process(r2);
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
	
	}
	
	public void writeMatrices(){
		if(q102out!=null){writeMatrix(q102out, q102GoodMatrix, q102BadMatrix, overwrite);}
		if(qbpout!=null){writeMatrix(qbpout, qbpGoodMatrix, qbpBadMatrix, overwrite);}
		if(q10out!=null){writeMatrix(q10out, q10GoodMatrix, q10BadMatrix, overwrite);}
		if(q12out!=null){writeMatrix(q12out, q12GoodMatrix, q12BadMatrix, overwrite);}
		if(qb012out!=null){writeMatrix(qb012out, qb012GoodMatrix, qb012BadMatrix, overwrite);}
		if(qb234out!=null){writeMatrix(qb234out, qb234GoodMatrix, qb234BadMatrix, overwrite);}
		if(qpout!=null){writeMatrix(qpout, qpGoodMatrix, qpBadMatrix, overwrite);}
		if(qout!=null){writeMatrix(qout, qGoodMatrix, qBadMatrix, overwrite);}
		if(pout!=null){writeMatrix(pout, pGoodMatrix, pBadMatrix, overwrite);}
		readstats.writeQualityToFile(qhist, false);
	}
	
	public static void writeMatrix(String fname, long[][][][] goodMatrix, long[][][][] badMatrix, boolean overwrite){
		assert(fname!=null) : "No file specified";
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, false, overwrite, false);
		TextStreamWriter tsw=new TextStreamWriter(ff, false);
		System.err.println("Starting tsw for "+fname);
		tsw.start();
		System.err.println("Started tsw for "+fname);
		StringBuilder sb=new StringBuilder();
		
		final int d0=goodMatrix.length, d1=goodMatrix[0].length, d2=goodMatrix[0][0].length, d3=goodMatrix[0][0][0].length;
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					for(int d=0; d<d3; d++){
						long good=goodMatrix[a][b][c][d];
						long bad=badMatrix[a][b][c][d];
						long sum=good+bad;
						if(sum>0){
							sb.append(a);
							sb.append('\t');
							sb.append(b);
							sb.append('\t');
							sb.append(c);
							sb.append('\t');
							sb.append(d);
							sb.append('\t');
							sb.append(sum);
							sb.append('\t');
							sb.append(bad);
							sb.append('\n');
						}
					}
					if(sb.length()>0){
						tsw.print(sb.toString());
						sb.setLength(0);
					}
				}
			}
		}
		System.err.println("Writing "+fname);
		tsw.poisonAndWait();
		System.err.println("Done.");
	}
	
	public static void writeMatrix(String fname, long[][][] goodMatrix, long[][][] badMatrix, boolean overwrite){
		assert(fname!=null) : "No file specified";
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, false, overwrite, false);
		TextStreamWriter tsw=new TextStreamWriter(ff, false);
		System.err.println("Starting tsw for "+fname);
		tsw.start();
		System.err.println("Started tsw for "+fname);
		StringBuilder sb=new StringBuilder();
		
		final int d0=goodMatrix.length, d1=goodMatrix[0].length, d2=goodMatrix[0][0].length;
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					long good=goodMatrix[a][b][c];
					long bad=badMatrix[a][b][c];
					long sum=good+bad;
					if(sum>0){
						sb.append(a);
						sb.append('\t');
						sb.append(b);
						sb.append('\t');
						sb.append(c);
						sb.append('\t');
						sb.append(sum);
						sb.append('\t');
						sb.append(bad);
						sb.append('\n');
					}
				}
				if(sb.length()>0){
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}
		}
		System.err.println("Writing "+fname);
		tsw.poisonAndWait();
		System.err.println("Done.");
	}
	
	public static void writeMatrix(String fname, long[][] goodMatrix, long[][] badMatrix, boolean overwrite){
		assert(fname!=null) : "No file specified";
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, false, overwrite, false);
		TextStreamWriter tsw=new TextStreamWriter(ff, false);
		System.err.println("Starting tsw for "+fname);
		tsw.start();
		System.err.println("Started tsw for "+fname);
		StringBuilder sb=new StringBuilder();
		
		final int d0=goodMatrix.length, d1=goodMatrix[0].length;
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				long good=goodMatrix[a][b];
				long bad=badMatrix[a][b];
				long sum=good+bad;
				if(sum>0){
					sb.append(a);
					sb.append('\t');
					sb.append(b);
					sb.append('\t');
					sb.append(sum);
					sb.append('\t');
					sb.append(bad);
					sb.append('\n');
				}
			}
			if(sb.length()>0){
				tsw.print(sb.toString());
				sb.setLength(0);
			}
		}
		System.err.println("Writing "+fname);
		tsw.poisonAndWait();
		System.err.println("Done.");
	}
	
	public static void writeMatrix(String fname, long[] goodMatrix, long[] badMatrix, boolean overwrite){
		assert(fname!=null) : "No file specified";
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, false, overwrite, false);
		TextStreamWriter tsw=new TextStreamWriter(ff, false);
		System.err.println("Starting tsw for "+fname);
		tsw.start();
		System.err.println("Started tsw for "+fname);
		StringBuilder sb=new StringBuilder();

		final int d0=goodMatrix.length;
		for(int a=0; a<d0; a++){
			long good=goodMatrix[a];
			long bad=badMatrix[a];
			long sum=good+bad;
			if(sum>0){
				sb.append(a);
				sb.append('\t');
				sb.append(sum);
				sb.append('\t');
				sb.append(bad);
				sb.append('\n');
			}
			if(sb.length()>0){
				tsw.print(sb.toString());
				sb.setLength(0);
			}
		}
		System.err.println("Writing "+fname);
		tsw.poisonAndWait();
		System.err.println("Done.");
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void process(Read r){
		if(r==null){return;}
		readsProcessed++;
		basesProcessed+=r.bases==null ? 0 : r.bases.length;
		
		if(verbose){outstream.println(r+"\n");}
		
		if(verbose){outstream.println("A");}
		if(r.match!=null && r.shortmatch()){
			r.match=Read.toLongMatchString(r.match);
			r.setShortMatch(false);
		}
		final byte[] quals=r.quality, bases=r.bases, match=r.match;
		if(quals==null || bases==null || match==null){return;}
		if(verbose){outstream.println("B");}
		if(r.containsNonNMS() || r.containsConsecutiveS(4)){
			if(verbose){System.err.println("*************************************************** "+new String(match));}
			return;
		}
		if(r.strand()==Gene.MINUS){
			Tools.reverseInPlace(match);
		}
		if(verbose){outstream.println("C");}
		
		final byte e='E';
		
		readstats.addToQualityHistogram(r);
		
		readsUsed++;
		for(int i=0, last=quals.length-1; i<quals.length; i++){
			if(verbose){outstream.print("D");}
			final byte q0=(i>0 ? (byte)Tools.mid(QMAX, quals[i-1], 0) : QEND);
			final byte q1=quals[i];
			final byte q2=(i<last ? (byte)Tools.mid(QMAX, quals[i+1], 0) : QEND);
			
			byte b0=i>1 ? bases[i-2] : e;
			byte b1=i>0 ? bases[i-1] : e;
			byte b2=bases[i];
			byte b3=i<last ? bases[i+1] : e;
			byte b4=i<last-1 ? bases[i+2] : e;
			byte n0=baseToNum[b0];
			byte n1=baseToNum[b1];
			byte n2=baseToNum[b2];
			byte n3=baseToNum[b3];
			byte n4=baseToNum[b4];
			byte m=match[i];
			
			if(m=='N' || !AminoAcid.isFullyDefined(b2)){
				if(verbose){outstream.print("E");}
				//do nothing
			}else{

				if(verbose){outstream.print("F");}
				basesUsed++;
				if(m=='m'){
					q102GoodMatrix[q1][q0][q2]++;
					qbpGoodMatrix[q1][n2][i]++;

					q10GoodMatrix[q1][q0]++;
					q12GoodMatrix[q1][q0]++;
					qb012GoodMatrix[q1][n0][n1][n2]++;
					qb234GoodMatrix[q1][n2][n3][n4]++;
					qpGoodMatrix[q1][i]++;
					qGoodMatrix[q1]++;
					pGoodMatrix[i]++;
				}else if(m=='S'){
					q102BadMatrix[q1][q0][q2]++;
					qbpBadMatrix[q1][n2][i]++;

					q10BadMatrix[q1][q0]++;
					q12BadMatrix[q1][q0]++;
					qb012BadMatrix[q1][n0][n1][n2]++;
					qb234BadMatrix[q1][n2][n3][n4]++;
					qpBadMatrix[q1][i]++;
					qBadMatrix[q1]++;
					pBadMatrix[i]++;
				}else{
					throw new RuntimeException("Bad symbol m='"+((char)m)+"'\n"+new String(match)+"\n"+new String(bases)+"\n");
				}
			}
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final void initializeMatrices(){
		initializeMatrices(q102, qbp, q10, q12, qb012, qb234, qp);
	}
	
	public static final void initializeMatrices(boolean q102, boolean qbp, boolean q10, boolean q12, boolean qb012, boolean qb234, boolean qp){
		if(initialized[0]){return;}
		
//		assert(false) : q102+". "+qbp+". "+q10+". "+q12+". "+qb012+". "+qb234+". "+qp;
		
		synchronized(initialized){
			if(initialized[0]){return;}

			if(q102){
				q102CountMatrix=loadMatrix(q102matrix, QMAX2, QMAX2, QMAX2);
				q102ProbMatrix=toProbs(q102CountMatrix[0], q102CountMatrix[1], OBSERVATION_CUTOFF);
			}
			if(qbp){
				qbpCountMatrix=loadMatrix(qbpmatrix, QMAX2, 4, LENMAX);
				qbpProbMatrix=toProbs(qbpCountMatrix[0], qbpCountMatrix[1], OBSERVATION_CUTOFF);
			}
			if(q10){
				q10CountMatrix=loadMatrix(q10matrix, QMAX2, QMAX2);
				q10ProbMatrix=toProbs(q10CountMatrix[0], q10CountMatrix[1], OBSERVATION_CUTOFF);
			}
			if(q12){
				q12CountMatrix=loadMatrix(q12matrix, QMAX2, QMAX2);
				q12ProbMatrix=toProbs(q12CountMatrix[0], q12CountMatrix[1], OBSERVATION_CUTOFF);
			}
			if(qb012){
				qb012CountMatrix=loadMatrix(qb012matrix, QMAX2, BMAX, BMAX, 4);
				qb012ProbMatrix=toProbs(qb012CountMatrix[0], qb012CountMatrix[1], OBSERVATION_CUTOFF);
			}
			if(qb234){
				qb234CountMatrix=loadMatrix(qb234matrix, QMAX2, 4, BMAX, BMAX);
				qb234ProbMatrix=toProbs(qb234CountMatrix[0], qb234CountMatrix[1], OBSERVATION_CUTOFF);
			}
			if(qp){
				qpCountMatrix=loadMatrix(qpmatrix, QMAX2, LENMAX);
				qpProbMatrix=toProbs(qpCountMatrix[0], qpCountMatrix[1], OBSERVATION_CUTOFF);
			}
			
			initialized[0]=true;
		}
	}
	
	public static final float estimateErrorProbAvg(byte[] quals, byte[] bases, int pos){
		if(q102ProbMatrix==null && qbpProbMatrix==null){return PROB_ERROR[quals[pos]];}
		
		final byte e='E';
		final int last=quals.length-1;
		
		final byte q0=(pos>0 ? (byte)Tools.mid(QMAX, quals[pos-1], 0) : QEND);
		final byte q1=quals[pos];
		final byte q2=(pos<last ? (byte)Tools.mid(QMAX, quals[pos+1], 0) : QEND);
		
		byte b0=pos>1 ? bases[pos-2] : e;
		byte b1=pos>0 ? bases[pos-1] : e;
		byte b2=bases[pos];
		byte b3=pos<last ? bases[pos+1] : e;
		byte b4=pos<last-1 ? bases[pos+2] : e;
		byte n0=baseToNum[b0];
		byte n1=baseToNum[b1];
		byte n2=baseToNum[b2];
		byte n3=baseToNum[b3];
		byte n4=baseToNum[b4];
		
		float expected=PROB_ERROR[q1];
		float sum=0;
		int x=0;

//		System.err.println();
//		System.err.println(((char)b0)+"\t"+((char)b1)+"\t"+((char)b2)+"\t"+((char)b3)+"\t"+((char)b4));
//		System.err.println((n0)+"\t"+(n1)+"\t"+(n2)+"\t"+(n3)+"\t"+(n4));
//		System.err.println(" "+"\t"+(q0)+"\t"+(q1)+"\t"+(q2)+"\t"+(" "));
//		System.err.println("Expected: "+expected);
		
		if(q102ProbMatrix!=null){
			float f=q102ProbMatrix[q1][q0][q2];
			sum+=f;
//			System.err.println(f);
			x++;
		}
		if(qbpProbMatrix!=null){
			float f=qbpProbMatrix[q1][n2][pos];
			sum+=f;
//			System.err.println(f);
			x++;
		}
		if(q10ProbMatrix!=null){
			float f=q10ProbMatrix[q1][q0];
			sum+=f;
//			System.err.println(f);
			x++;
		}
		if(q12ProbMatrix!=null){
			float f=q12ProbMatrix[q1][q2];
			sum+=f;
//			System.err.println(f);
			x++;
		}
		if(qb012ProbMatrix!=null){
			float f=qb012ProbMatrix[q1][n0][n1][n2];
			sum+=f;
//			System.err.println(f);
			x++;
		}
		if(qb234ProbMatrix!=null){
			float f=qb234ProbMatrix[q1][n2][n3][n4];
			sum+=f;
//			System.err.println(f);
			x++;
		}
		if(qpProbMatrix!=null){
			float f=qpProbMatrix[q1][pos];
			sum+=f;
//			System.err.println(f);
			x++;
		}
//		System.err.println("result: "+sum+", "+x+", "+sum/(double)x);
//		
//		assert(pos<149) : sum+", "+x+", "+sum/(double)x;
		
		if(x<1){
			assert(false);
			return expected;
		}
		return (float)(sum/(double)x);
	}
	
	public static final float estimateErrorProbGeoAvg(byte[] quals, byte[] bases, int pos){
//		if(q102ProbMatrix==null && qbpProbMatrix==null){return PROB_ERROR[quals[pos]];}
		
		final byte e='E';
		final int last=quals.length-1;
		
		final byte q0=(pos>0 ? (byte)Tools.mid(QMAX, quals[pos-1], 0) : QEND);
		final byte q1=quals[pos];
		final byte q2=(pos<last ? (byte)Tools.mid(QMAX, quals[pos+1], 0) : QEND);
		
		byte b0=pos>1 ? bases[pos-2] : e;
		byte b1=pos>0 ? bases[pos-1] : e;
		byte b2=bases[pos];
		byte b3=pos<last ? bases[pos+1] : e;
		byte b4=pos<last-1 ? bases[pos+2] : e;
		byte n0=baseToNum[b0];
		byte n1=baseToNum[b1];
		byte n2=baseToNum[b2];
		byte n3=baseToNum[b3];
		byte n4=baseToNum[b4];
		
		float expected=PROB_ERROR[q1];
		double product=1;
		int x=0;

//		System.err.println();
//		System.err.println(((char)b0)+"\t"+((char)b1)+"\t"+((char)b2)+"\t"+((char)b3)+"\t"+((char)b4));
//		System.err.println((n0)+"\t"+(n1)+"\t"+(n2)+"\t"+(n3)+"\t"+(n4));
//		System.err.println(" "+"\t"+(q0)+"\t"+(q1)+"\t"+(q2)+"\t"+(" "));
//		System.err.println("Expected: "+expected);
		
		if(q102ProbMatrix!=null){
			float f=q102ProbMatrix[q1][q0][q2];
			product*=f;
//			System.err.println(f);
			x++;
		}
		if(qbpProbMatrix!=null){
			float f=qbpProbMatrix[q1][n2][pos];
			product*=f;
//			System.err.println(f);
			x++;
		}
		if(q10ProbMatrix!=null){
			float f=q10ProbMatrix[q1][q0];
			product*=f;
//			System.err.println(f);
			x++;
		}
		if(q12ProbMatrix!=null){
			float f=q12ProbMatrix[q1][q2];
			product*=f;
//			System.err.println(f);
			x++;
		}
		if(qb012ProbMatrix!=null){
			float f=qb012ProbMatrix[q1][n0][n1][n2];
			product*=f;
//			System.err.println(f);
			x++;
		}
		if(qb234ProbMatrix!=null){
			float f=qb234ProbMatrix[q1][n2][n3][n4];
			product*=f;
//			System.err.println(f);
			x++;
		}
		if(qpProbMatrix!=null){
			float f=qpProbMatrix[q1][pos];
			product*=f;
//			System.err.println(f);
			x++;
		}
//		System.err.println("result: "+product+", "+x+", "+(float)Math.pow(product, 1.0/x));
		
//		assert(pos<149) : product+", "+x+", "+(float)Math.pow(product, 1.0/x);
		
		if(x<1){
			assert(false);
			return expected;
		}
		return (float)Math.pow(product, 1.0/x);
	}
	
	public static final float estimateErrorProb2(byte[] quals, byte[] bases, int pos){
		if(q102ProbMatrix==null && qbpProbMatrix==null){return PROB_ERROR[quals[pos]];}
		
		final byte e='E';
		final int last=quals.length-1;
		
		final byte q0=(pos>0 ? (byte)Tools.mid(QMAX, quals[pos-1], 0) : QEND);
		final byte q1=quals[pos];
		final byte q2=(pos<last ? (byte)Tools.mid(QMAX, quals[pos+1], 0) : QEND);
		
		byte b0=pos>1 ? bases[pos-2] : e;
		byte b1=pos>0 ? bases[pos-1] : e;
		byte b2=bases[pos];
		byte b3=pos<last ? bases[pos+1] : e;
		byte b4=pos<last-1 ? bases[pos+2] : e;
		byte n0=baseToNum[b0];
		byte n1=baseToNum[b1];
		byte n2=baseToNum[b2];
		byte n3=baseToNum[b3];
		byte n4=baseToNum[b4];
		
		long sum=OBSERVATION_CUTOFF, bad=0;
		if(q102CountMatrix!=null){
			sum+=q102CountMatrix[0][q1][q0][q2];
			bad+=q102CountMatrix[1][q1][q0][q2];
		}
		if(qbpCountMatrix!=null){
			sum+=qbpCountMatrix[0][q1][n2][pos];
			bad+=qbpCountMatrix[1][q1][n2][pos];
		}
		if(q10CountMatrix!=null){
			sum+=q10CountMatrix[0][q1][q0];
			bad+=q10CountMatrix[1][q1][q0];
		}
		if(q12CountMatrix!=null){
			sum+=q12CountMatrix[0][q1][q2];
			bad+=q12CountMatrix[1][q1][q2];
		}
		if(qb012CountMatrix!=null){
			sum+=qb012CountMatrix[0][q1][n0][n1][n2];
			bad+=qb012CountMatrix[1][q1][n0][n1][n2];
		}
		if(qb234CountMatrix!=null){
			sum+=qb234CountMatrix[0][q1][n2][n3][n4];
			bad+=qb234CountMatrix[1][q1][n2][n3][n4];
		}
		if(qpCountMatrix!=null){
			sum+=qpCountMatrix[0][q1][pos];
			bad+=qpCountMatrix[1][q1][pos];
		}

		double expected=PROB_ERROR[q1];
		double dbad=bad+expected*OBSERVATION_CUTOFF;
		
		double observed=dbad/sum;
		
		return (float)Math.sqrt(observed*expected);
	}
	
	/*--------------------------------------------------------------*/
	
	private static double modify(final double sum, final double bad, final int phred, final long cutoff){
		double expected=QualityTools.phredToProbError(phred);

		double sum2=sum+cutoff;
		double bad2=bad+expected*cutoff;
		double measured=bad2/sum2;

		double modified=Math.pow(measured*measured*measured*expected, 0.25);
//		double modified=Math.sqrt(measured*expected);
//		double modified=(measured+expected)*.5;
		
		return modified;
	}
	
	public static final float[][][][] toProbs(long[][][][] sumMatrix, long[][][][] badMatrix, final long cutoff){
		final int d0=sumMatrix.length, d1=sumMatrix[0].length, d2=sumMatrix[0][0].length, d3=sumMatrix[0][0][0].length;
		float[][][][] probs=new float[d0][d1][d2][d3];
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					for(int d=0; d<d3; d++){
						double sum=sumMatrix[a][b][c][d];
						double bad=badMatrix[a][b][c][d];
						double modified=modify(sum, bad, a, cutoff);
						probs[a][b][c][d]=(float)modified;
					}
				}
			}
		}
		return probs;
	}
	
	public static final float[][][] toProbs(long[][][] sumMatrix, long[][][] badMatrix, final long cutoff){
		final int d0=sumMatrix.length, d1=sumMatrix[0].length, d2=sumMatrix[0][0].length;
		float[][][] probs=new float[d0][d1][d2];
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					double sum=sumMatrix[a][b][c];
					double bad=badMatrix[a][b][c];
					double modified=modify(sum, bad, a, cutoff);
					probs[a][b][c]=(float)modified;
				}
			}
		}
		return probs;
	}
	
	public static final float[][] toProbs(long[][] sumMatrix, long[][] badMatrix, final long cutoff){
		final int d0=sumMatrix.length, d1=sumMatrix[0].length;
		float[][] probs=new float[d0][d1];
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
					double sum=sumMatrix[a][b];
					double bad=badMatrix[a][b];
					double modified=modify(sum, bad, a, cutoff);
					probs[a][b]=(float)modified;
			}
		}
		return probs;
	}
	
	/*--------------------------------------------------------------*/

	public static final long[][] loadMatrix(String fname, int d0){
		if(fname==null){return null;}
		fname=Data.findPath(fname);
		long[][] matrix=new long[2][d0];

		TextFile tf=new TextFile(fname, false, false);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			String[] split=line.split("\t");
			assert(split.length==3) : Arrays.toString(split);
			int a=Integer.parseInt(split[0]);
			long bases=Integer.parseInt(split[1]);
			long errors=Integer.parseInt(split[2]);
			matrix[0][a]=bases;
			matrix[1][a]=errors;
		}
		return matrix;
	}

	public static final long[][][] loadMatrix(String fname, int d0, int d1){
		if(fname==null){return null;}
		fname=Data.findPath(fname);
		long[][][] matrix=new long[2][d0][d1];

		TextFile tf=new TextFile(fname, false, false);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			String[] split=line.split("\t");
			assert(split.length==4) : Arrays.toString(split);
			int a=Integer.parseInt(split[0]);
			int b=Integer.parseInt(split[1]);
			long bases=Integer.parseInt(split[2]);
			long errors=Integer.parseInt(split[3]);
			matrix[0][a][b]=bases;
			matrix[1][a][b]=errors;
		}
		return matrix;
	}

	public static final long[][][][] loadMatrix(String fname, int d0, int d1, int d2){
		if(fname==null){return null;}
		fname=Data.findPath(fname);
		long[][][][] matrix=new long[2][d0][d1][d2];

		TextFile tf=new TextFile(fname, false, false);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			String[] split=line.split("\t");
			assert(split.length==5) : Arrays.toString(split);
			int a=Integer.parseInt(split[0]);
			int b=Integer.parseInt(split[1]);
			int c=Integer.parseInt(split[2]);
			long bases=Integer.parseInt(split[3]);
			long errors=Integer.parseInt(split[4]);
			matrix[0][a][b][c]=bases;
			matrix[1][a][b][c]=errors;
		}
		return matrix;
	}

	public static final long[][][][][] loadMatrix(String fname, int d0, int d1, int d2, int d3){
		if(fname==null){return null;}
		fname=Data.findPath(fname);
		long[][][][][] matrix=new long[2][d0][d1][d2][d3];

		TextFile tf=new TextFile(fname, false, false);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			String[] split=line.split("\t");
			assert(split.length==6) : Arrays.toString(split);
			int a=Integer.parseInt(split[0]);
			int b=Integer.parseInt(split[1]);
			int c=Integer.parseInt(split[2]);
			int d=Integer.parseInt(split[3]);
			long bases=Integer.parseInt(split[4]);
			long errors=Integer.parseInt(split[5]);
			matrix[0][a][b][c][d]=bases;
			matrix[1][a][b][c][d]=errors;
		}
		return matrix;
	}
	
	private static byte[] fillBaseToNum(){
		byte[] btn=new byte[128];
		Arrays.fill(btn, (byte)5);
		btn['A']=btn['a']=0;
		btn['C']=btn['c']=1;
		btn['G']=btn['g']=2;
		btn['T']=btn['t']=3;
		btn['U']=btn['u']=3;
		btn['E']=4;
		return btn;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ReadStats readstats=new ReadStats();
	
	private long[][][] q102GoodMatrix=new long[QMAX2][QMAX2][QMAX2];
	private long[][][] q102BadMatrix=new long[QMAX2][QMAX2][QMAX2];

	private long[][][] qbpGoodMatrix=new long[QMAX2][BMAX][LENMAX];
	private long[][][] qbpBadMatrix=new long[QMAX2][BMAX][LENMAX];

	private long[][] q10GoodMatrix=new long[QMAX2][QMAX2];
	private long[][] q10BadMatrix=new long[QMAX2][QMAX2];

	private long[][] q12GoodMatrix=new long[QMAX2][QMAX2];
	private long[][] q12BadMatrix=new long[QMAX2][QMAX2];

	private long[][][][] qb012GoodMatrix=new long[QMAX2][BMAX][BMAX][BMAX];
	private long[][][][] qb012BadMatrix=new long[QMAX2][BMAX][BMAX][BMAX];

	private long[][][][] qb234GoodMatrix=new long[QMAX2][BMAX][BMAX][BMAX];
	private long[][][][] qb234BadMatrix=new long[QMAX2][BMAX][BMAX][BMAX];

	private long[][] qpGoodMatrix=new long[QMAX2][LENMAX];
	private long[][] qpBadMatrix=new long[QMAX2][LENMAX];

	private long[] qGoodMatrix=new long[QMAX2];
	private long[] qBadMatrix=new long[QMAX2];

	private long[] pGoodMatrix=new long[LENMAX];
	private long[] pBadMatrix=new long[LENMAX];
	
	private PrintStream outstream=System.err;
	private boolean verbose=false;
	private long maxReads=-1;
	private String[] in;
	
	private String q102out="q102matrix.txt.gz";
	private String qbpout="qbpmatrix.txt.gz";
	private String q10out="q10matrix.txt.gz";
	private String q12out="q12matrix.txt.gz";
	private String qb012out="qb012matrix.txt.gz";
	private String qb234out="qb234matrix.txt.gz";
	private String qpout="qpmatrix.txt.gz";
	private String qout="qmatrix.txt.gz";
	private String pout="pmatrix.txt.gz";
	private String qhist="qhist.txt";
	
	private boolean overwrite=false;
	private long readsProcessed=0;
	private long basesProcessed=0;
	private long readsUsed=0;
	private long basesUsed=0;
	private boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final boolean[] initialized={false};
	
	private static final int QMAX=41;
	private static final int QEND=QMAX+1;
	private static final int QMAX2=QEND+1;
	private static final int BMAX=6;
	private static final int LENMAX=400;
	private static final byte[] baseToNum=fillBaseToNum();
	private static final byte[] numToBase={'A', 'C', 'G', 'T', 'E', 'N'};
	private static final float[] PROB_ERROR=QualityTools.PROB_ERROR;
	
	public static String q102matrix="?q102matrix.txt.gz";
	public static String qbpmatrix="?qbpmatrix.txt.gz";
	public static String q10matrix="?q10matrix.txt.gz";
	public static String q12matrix="?q12matrix.txt.gz";
	public static String qb012matrix="?qb012matrix.txt.gz";
	public static String qb234matrix="?qb234matrix.txt.gz";
	public static String qpmatrix="?qpmatrix.txt.gz";
	
	public static long[][][][] q102CountMatrix;
	public static long[][][][] qbpCountMatrix;
	
	public static long[][][] q10CountMatrix;
	public static long[][][] q12CountMatrix;
	public static long[][][][][] qb012CountMatrix;
	public static long[][][][][] qb234CountMatrix;
	public static long[][][] qpCountMatrix;

	public static float[][][] q102ProbMatrix;
	public static float[][][] qbpProbMatrix;
	
	public static float[][] q10ProbMatrix;
	public static float[][] q12ProbMatrix;
	public static float[][][][] qb012ProbMatrix;
	public static float[][][][] qb234ProbMatrix;
	public static float[][] qpProbMatrix;
	
	public static boolean q102=true;
	public static boolean qbp=true;
	public static boolean q10=true;
	public static boolean q12=true;
	public static boolean qb012=true;
	public static boolean qb234=true;
	public static boolean qp=true;
	
	public static final long OBSERVATION_CUTOFF=10000; //Soft threshold
	
}
