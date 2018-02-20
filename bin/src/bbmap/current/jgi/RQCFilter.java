package jgi;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.TimeZone;

import dna.Data;

import stream.FASTQ;
import stream.Read;

import align2.BBMap;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;
import fileIO.ByteFile1;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;

/**
 * Wrapper for BBDukF to implement Rolling QC's filter stage.
 * @author Brian Bushnell
 * @date Nov 26, 2013
 *
 */
public class RQCFilter {

	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Methods    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Program entrance from command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Create a filter instance
		RQCFilter filter=new RQCFilter(args);
		
		///...and execute it.
		filter.process();
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	RQCFilter(String[] args){
		
		//Optional default parameters to match current pipeline
//		arglist.add("k=22");
//		arglist.add("maxbadkmers=2");
		
		//Symbols to insert in output filename to denote operations performed; may be overriden from command line
		String symbols_=null;//"filtered"
		
		//Parse argument list
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("="); //Expect key=value pairs
			String a=split[0].toLowerCase(); //All keys are converted to lower case
			String b=split.length>1 ? split[1] : null;
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens
			
//			System.out.println("Processing '"+arg+"' a='"+a+"', b='"+b+"'");
			
			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("null") || a.equals(in2)){
				// do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				primaryArgList.add(arg);
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
			}else if(a.equals("ref")){
				if(b!=null){
					if(!b.contains(",") || new File(b).exists()){
						refs.add(b);
					}else{
						String[] split2=b.split(",");
						for(String s2 : split2){
							refs.add(s2);
						}
					}
				}
			}else if(a.equals("artifactdb")){
				mainArtifactFile=b;
			}else if(a.equals("rnadb")){
				artifactFileRna=b;
			}else if(a.equals("dnadb")){
				artifactFileDna=b;
			}else if(a.equals("phixref")){
				phixRef=b;
			}else if(a.equals("fragadapter")){
				fragAdapter=b;
			}else if(a.equals("lfpelinker")){
				lfpeLinker=b;
			}else if(a.equals("cliplinker") || a.equals("jointseq")){
				clipLinker=b;
			}else if(a.equals("clrslinker")){
				clrsLinker=b;
			}else if(a.equals("trimfragadapter")){
				fragAdapterFlag=Tools.parseBoolean(b);
			}else if(a.equals("removehuman")){
				humanFlag=Tools.parseBoolean(b);
			}else if(a.equals("useindex")){
				humanRefIndexedFlag=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("mlf") || a.equals("minlenfrac") || a.equals("minlenfraction") || a.equals("minlengthfraction")){
				minLenFraction=Float.parseFloat(b);
			}else if(a.equals("libtype") || a.equals("library")){
				libType=toLibType(b);
			}else if(a.equals("path") || a.equals("outdir")){
				outDir=b;
			}else if(a.equals("symbols")){
				symbols_=b;
			}else if(a.equals("overallstats") || a.equals("stats")){
				rqcStatsName=b;
			}else if(a.equals("scafstats")){
				scaffoldStatsName=b;
			}else if(a.equals("kmerstats")){
				kmerStatsName=b;
			}else if(a.equals("log")){
				logName=b;
			}else if(a.equals("filelist")){
				fileListName=b;
			}else if(a.equals("compress")){
				compress=Tools.parseBoolean(b);
			}else if(a.equals("rna")){
				rnaFlag=Tools.parseBoolean(b);
			}else if(a.equals("phix")){
				phixFlag=Tools.parseBoolean(b);
			}else if(a.equals("jointseq")){
				jointSeq=b;
			}else if(a.equals("ktrim")){
				ktrim=b;
			}else if(a.equals("mink")){
				mink=Integer.parseInt(b);
			}else if(a.equals("maq")){
				maq=Byte.parseByte(b);
			}else if(a.equals("trimq")){
				trimq=Byte.parseByte(b);
			}else if(a.equals("qtrim")){
				if(b==null){qtrim="rl";}
				else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){qtrim="l";}
				else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){qtrim="r";}
				else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){qtrim="lr";}
				else if(Character.isDigit(b.charAt(0))){
					trimq=Byte.parseByte(b);
					qtrim=(trimq>=0 ? "lr" : "f");
				}else{qtrim=""+Tools.parseBoolean(b);}
			}else if(a.equals("optitrim") || a.equals("otf") || a.equals("otm")){
				if(b!=null && (b.charAt(0)=='.' || Character.isDigit(b.charAt(0)))){
					TrimRead.optimalMode=true;
					TrimRead.optimalBias=Float.parseFloat(b);
					assert(TrimRead.optimalBias>=0 && TrimRead.optimalBias<1);
				}else{
					TrimRead.optimalMode=Tools.parseBoolean(b);
				}
			}else if(a.equals("maxns")){
				maxNs=Integer.parseInt(b);
			}else if(a.equals("usetmpdir")){
				writeTempToTmpdir=Tools.parseBoolean(b);
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
				if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
				else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Tools.parseBoolean(b);
					System.err.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}else if(a.equals("tuc") || a.equals("touppercase")){
				Read.TO_UPPER_CASE=Tools.parseBoolean(b);
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else{
				//Uncaptured arguments are passed to BBDuk
				primaryArgList.add(arg);
			}
		}
		
//		assert(false) : rnaFlag+"\n"+primaryArgList+"\n"+libType+"\n"+outDir;
		
		if(writeTempToTmpdir){tmpDir=Shared.TMPDIR;}
		else{tmpDir=null;}
		
		//Set final field 'symbols'
		symbols=(symbols_==null ? abbreviation() : symbols_);
		
		//Pass overwrite flag to BBDuk
		primaryArgList.add("ow="+overwrite);
		
		if(outDir!=null){
			outDir=outDir.trim().replace('\\', '/');
			if(outDir.length()>0 && !outDir.endsWith("/")){outDir=outDir+"/";}
		}else{outDir="";}
		
		{//Prepend output directory to output files
			if(logName!=null){logName=outDir+logName+".tmp";} //Add '.tmp' to log file
			if(reproduceName!=null){reproduceName=outDir+reproduceName;}
			if(fileListName!=null){fileListName=outDir+fileListName;}
		}

		{//Create unique output file names for second pass
			if(rqcStatsName!=null){
				rqcStatsName_kt=outDir+"ktrim_"+rqcStatsName;
				rqcStatsName=outDir+rqcStatsName;
			}
			if(kmerStatsName!=null){
				kmerStatsName_kt=outDir+"ktrim_"+kmerStatsName;
				kmerStatsName=outDir+kmerStatsName;
			}
			if(scaffoldStatsName!=null){
				scaffoldStatsName_kt=outDir+"ktrim_"+scaffoldStatsName;
				scaffoldStatsName=outDir+scaffoldStatsName;
			}
		}
		
		//Create output filename from input filename if no output filename is specified
		if(out1==null && in1!=null){
			File f=new File(in1);
			String name=f.getName();
			String raw=ReadWrite.rawName(name);
			int x=raw.lastIndexOf('.');
			if(x>-1){
				out1=raw.substring(0, x)+"."+symbols+raw.substring(x)+(compress ? ".gz" : "");
			}else{
				out1=raw+"."+symbols+".fastq"+(compress ? ".gz" : "");
			}
		}
		
		tempSalt=KmerNormalize.getSalt(out1, 1);
		trimPrefix="TEMP_TRIM_"+tempSalt+"_";
		humanPrefix="TEMP_HUMAN_"+tempSalt+"_";
		filterPrefix="TEMP_FILTER_"+tempSalt+"_";
	}

	
	/*--------------------------------------------------------------*/
	/*----------------     Processing Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Primary method to fully execute the program.
	 */
	public void process(){
		
		//Create output directory
		if(outDir!=null && outDir.length()>0){
			File f=new File(outDir);
			if(!f.exists()){
				f.mkdirs();
			}
		}
		
		//Create log file
		if(logName!=null){
			boolean b=Tools.canWrite(logName, overwrite);
			assert(b) : "Can't write to "+logName;
			log("start", false);
		}
		
		//Create file list file
		if(fileListName!=null){
			boolean b=Tools.canWrite(fileListName, overwrite);
			assert(b) : "Can't write to "+fileListName;
			
			StringBuilder sb=new StringBuilder();
			if(out1!=null){sb.append("filtered_fastq="+out1).append('\n');}
			if(qfout1!=null){sb.append("filtered_qual="+qfout1).append('\n');}
			if(out2!=null){sb.append("filtered_fastq_2="+out2).append('\n');}
			if(qfout2!=null){sb.append("filtered_qual_2="+qfout2).append('\n');}
			
			if(sb.length()>0){
				ReadWrite.writeString(sb, fileListName, false);
			}
		}

		final boolean doFilter;
		final boolean doTrim;
		final boolean doHuman=humanFlag;
		
		//Determine execution path
		if(libType==FRAG || ((libType==LFPE && lfpeLinker==null) || (libType==CLIP && clipLinker==null) || (libType==CLRS && clrsLinker==null))){
			doTrim=fragAdapterFlag;
			doFilter=true;
		}else if(libType==LFPE){
			doTrim=true;
			doFilter=true;
		}else if(libType==CLIP){
			doTrim=true;
			doFilter=true;
		}else if(libType==CLRS){
			doTrim=true;
			doFilter=true;
		}else{
			throw new RuntimeException("Unknown library type.");
		}
		
		{
			int step=0;
			final int numSteps=(doFilter ? 1 : 0)+(doTrim ? 1 : 0)+(doHuman ? 1 : 0);
			String inPrefix=null, outPrefix=null;
			if(doTrim){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? trimPrefix : null);
				if(step==1){
					trim(in1, in2, out1, out2, qfin1, qfin2, qfout1, qfout2, inPrefix, outPrefix);
				}else{
					trim(out1, out2, out1, out2, qfout1, qfout2, qfout1, qfout2, inPrefix, outPrefix);
				}
				if(inPrefix!=null){
					delete(inPrefix, out1, out2, qfout1, qfout2);
				}
			}
			
			if(doFilter){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? filterPrefix : null);
				if(step==1){
					filter(in1, in2, out1, out2, qfin1, qfin2, qfout1, qfout2, inPrefix, outPrefix);
				}else{
					filter(out1, out2, out1, out2, qfout1, qfout2, qfout1, qfout2, inPrefix, outPrefix);
				}
				if(step>1){
					delete(inPrefix, out1, out2, qfout1, qfout2);
				}
			}
			
			if(doHuman){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? humanPrefix : null);
				if(step==1){
					dehumanize(in1, in2, out1, out2, qfin1, qfin2, qfout1, qfout2, inPrefix, outPrefix);
				}else{
					dehumanize(out1, out2, out1, out2, qfout1, qfout2, qfout1, qfout2, inPrefix, outPrefix);
				}
				Data.unloadAll();
				if(step>1){
					delete(inPrefix, out1, out2, qfout1, qfout2);
				}
			}
		}
		
		//Write combined stats file (number of reads/bases present/removed in each stage) 
		if(rqcStatsName!=null){
			final TextStreamWriter tsw=new TextStreamWriter(rqcStatsName, overwrite, false, false);
			tsw.start();
			tsw.println(BBDukF.rqcString());
			tsw.poisonAndWait();
		}
		
//		{//Set files to permission 777
//			setPermissions((out1==null ? null : outDir+out1),(out2==null ? null : outDir+out2));
//			setPermissions((qfout1==null ? null : outDir+qfout1),(qfout2==null ? null : outDir+qfout2));
//			setPermissions(reproduceName,fileListName);
//			setPermissions(rqcStatsName,kmerStatsName,scaffoldStatsName);
//			setPermissions(rqcStatsName_kt,kmerStatsName_kt,scaffoldStatsName_kt);
//			setPermissions(outDir);
//		}
		
		//Finish writing log file
		if(logName!=null){
			log("complete", true);
			if(logName.endsWith(".tmp")){ //Remove .tmp extension
				String old=logName;
				logName=logName.substring(0, logName.length()-4);
				new File(old).renameTo(new File(logName));
			}
		}
		
//		//Set log file permission
//		setPermissions(logName);
		
	}
	
	/**
	 * Runs BBMap to perform:
	 * Human contaminant removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param qfin1 Primary input qual file
	 * @param qfin2 Secondary input qual file
	 * @param qfout1 Primary output qual file
	 * @param qfout2 Secondary output qual file
	 * @param inPrefix Append this prefix to input filenames
	 */
	private void dehumanize(String in1, String in2, String out1, String out2, String qfin1, String qfin2, String qfout1, String qfout2, String inPrefix, String outPrefix){
		
		log("dehumanize start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		{
//			argList.add("kfilter="+47);
			argList.add("minratio=.75");
			argList.add("maxindel=20");
			argList.add("bw=20");
			argList.add("bwr=0.18");
			argList.add("minhits=2");
			if(humanRefIndexedFlag){
				argList.add("path="+humanPath);
			}else{
				argList.add("ref="+humanRef);
				argList.add("nodisk");
			}
			argList.add("quickmatch");
			argList.add("overwrite="+overwrite);
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("outu1="+outPre+out1);}
			if(out2!=null){argList.add("outu2="+outPre+out2);}
			if(qfin1!=null){argList.add("qfin1="+inPre+qfin1);}
			if(qfin2!=null){argList.add("qfin2="+inPre+qfin2);}
			if(qfout1!=null){argList.add("qfoutu1="+outPre+qfout1);}
			if(qfout2!=null){argList.add("qfoutu2="+outPre+qfout2);}
			
		}
		
		String[] args=argList.toArray(new String[0]);
		
		{//Run BBMap
			try {
				BBMap.main(args);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbmap.sh", args, true, overwrite);
		}
		
		//Optionally append files to file list here
		
		log("dehumanize finish", true);
	}
	
	/**
	 * Runs BBDuk to perform:
	 * Quality filtering, quality trimming, n removal, short read removal, artifact removal (via kmer filtering), phiX removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param qfin1 Primary input qual file
	 * @param qfin2 Secondary input qual file
	 * @param qfout1 Primary output qual file
	 * @param qfout2 Secondary output qual file
	 * @param inPrefix Append this prefix to input filenames
	 */
	private void filter(String in1, String in2, String out1, String out2, String qfin1, String qfin2, String qfout1, String qfout2, String inPrefix, String outPrefix){
		
		log("filter start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		{//Fill list with BBDuk arguments
			if(maq>-1){argList.add("maq="+maq);}
			if(qtrim!=null){
				argList.add("trimq="+trimq);
				argList.add("qtrim="+qtrim);
			}
			argList.add("overwrite="+overwrite);
			if(maxNs>=0){argList.add("maxns="+maxNs);}
			if(minLen>0){argList.add("minlen="+minLen);}
			if(minLenFraction>0){argList.add("minlenfraction="+minLenFraction);}
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}
			if(qfin1!=null){argList.add("qfin1="+inPre+qfin1);}
			if(qfin2!=null){argList.add("qfin2="+inPre+qfin2);}
			if(qfout1!=null){argList.add("qfout1="+outPre+qfout1);}
			if(qfout2!=null){argList.add("qfout2="+outPre+qfout2);}

//			if(rqcStatsName!=null){al.add("rqc="+rqcStatsName);} //Old style for 2 log files
			if(rqcStatsName!=null){argList.add("rqc=hashmap");}
			if(kmerStatsName!=null){argList.add("outduk="+kmerStatsName);}
			if(scaffoldStatsName!=null){argList.add("stats="+scaffoldStatsName);}
		}
		
		{//Add BBDuk references
			refs.add(mainArtifactFile);
			refs.add(rnaFlag ? artifactFileRna : artifactFileDna);
			
			if(phixFlag){refs.add(phixRef);}

			if(libType==FRAG){

			}else if(libType==LFPE){

			}else if(libType==CLIP){

			}else if(libType==CLRS){

			}else{
				throw new RuntimeException("Unknown library type.");
			}

			StringBuilder refstring=new StringBuilder();
			for(String ref : refs){
				if(ref!=null){
					refstring.append(refstring.length()==0 ? "ref=" : ",");
					refstring.append(ref);
				}
			}

			if(refstring!=null && refstring.length()>0){
				argList.add(refstring.toString());
			}
		}
		
		String[] dukargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbduk.sh", dukargs, true, overwrite);
		}
		
		{//Run BBDuk
			BBDukF duk=new BBDukF(dukargs);
			try {
				duk.process();
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("filter finish", true);
	}
	
	
	/**
	 * Runs BBDuk to perform:
	 * Kmer trimming, short read removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param qfin1 Primary input qual file
	 * @param qfin2 Secondary input qual file
	 * @param qfout1 Primary output qual file
	 * @param qfout2 Secondary output qual file
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void trim(String in1, String in2, String out1, String out2, String qfin1, String qfin2, String qfout1, String qfout2, String inPrefix, String outPrefix){
		
		log("ktrim start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		{//Fill list with BBDuk arguments
			argList.add("mink="+mink);
			argList.add("ktrim="+(ktrim==null ? "f" : ktrim));
			argList.add("overwrite="+overwrite);
			if(minLen>0){argList.add("minlen="+minLen);}
			if(minLenFraction>0){argList.add("minlenfraction="+minLenFraction);}
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}
			if(qfin1!=null){argList.add("qfin1="+inPre+qfin1);}
			if(qfin2!=null){argList.add("qfin2="+inPre+qfin2);}
			if(qfout1!=null){argList.add("qfout1="+outPre+qfout1);}
			if(qfout2!=null){argList.add("qfout2="+outPre+qfout2);}

//			if(rqcStatsName!=null){al.add("rqc="+rqcStatsName_kt);} //Old style for 2 log files
			if(rqcStatsName!=null){argList.add("rqc=hashmap");}
			if(kmerStatsName!=null){argList.add("outduk="+kmerStatsName_kt);}
			if(scaffoldStatsName!=null){argList.add("stats="+scaffoldStatsName_kt);}
		}
		
		{//Add BBDuk references
			ArrayList<String> refs=new ArrayList<String>();

			if(libType==FRAG){
				refs.add(fragAdapter);
			}else if(libType==LFPE){
				refs.add(lfpeLinker);
			}else if(libType==CLIP){
//				refs.add(clipLinker);
				if(clipLinker!=null){
					argList.add("literal="+clipLinker);
					{//Special processing for literal strings of approx 4bp
						String[] split=clipLinker.split(",");
						int min=split[0].length();
						for(String s : split){min=Tools.min(min, s.length());}
						argList.add("k="+min);
						argList.add("mink=-1");
						argList.add("mm=f");
						argList.add("hdist=0");
						argList.add("edist=0");
						argList.add("ktrimexclusive=t");
					}
				}else{
					throw new RuntimeException("Null clip linker.");
				}
			}else if(libType==CLRS){
				refs.add(clrsLinker);
			}else{
				throw new RuntimeException("Unknown library type.");
			}

			StringBuilder refstring=new StringBuilder();
			for(String ref : refs){
				if(ref!=null){
					refstring.append(refstring.length()==0 ? "ref=" : ",");
					refstring.append(ref);
				}
			}

			if(refstring!=null && refstring.length()>0){
				argList.add(refstring.toString());
			}
		}
		
		String[] dukargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbduk.sh", dukargs, false, overwrite);
		}
		
		{//run BBDuk
			BBDukF duk=new BBDukF(dukargs);
			try {
				duk.process();
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("ktrim finish", true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Log a message in the log file
	 * @param message Message to log
	 * @param append True to append, false to overwrite
	 */
	private void log(String message, boolean append){
		if(logName!=null){
			ReadWrite.writeString(message+", "+timeString()+"\n", logName, append);
		}
	}
	
	
	/**
	 * Delete all non-null filenames.
	 * @param prefix Append this prefix to filenames before attempting to delete them
	 * @param names Filenames to delete
	 */
	private void delete(String prefix, String...names){
		log("delete temp files start", true);
		if(names!=null){
			final String pre=(prefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+prefix);
			for(String s : names){
				if(s!=null){
					s=pre+s;
					if(verbose){System.err.println("Trying to delete "+s);}
					File f=new File(s);
					if(f.exists()){
						f.delete();
						writeReproduceFile(reproduceName, "rm", new String[] {s}, true, overwrite);
					}
				}
			}
		}
		log("delete temp files finish", true);
	}
	
	/**
	 * @return String of symbols indicating which processes were applied to the input reads
	 */
	private String abbreviation(){
		StringBuilder sb=new StringBuilder();
		
		if(mainArtifactFile!=null || (rnaFlag ? artifactFileRna!=null : artifactFileDna!=null)){sb.append("a");}
		
		if(maxNs>=0){sb.append("n");}
//		if(qtrim!=null && !qtrim.equalsIgnoreCase("f") && !qtrim.equalsIgnoreCase("false")){sb.append("q");}
		if(maq>0){sb.append("q");}
		
		if(rnaFlag){sb.append("r");}
		else{sb.append("d");}
		
		if(libType==CLIP){sb.append("c");}
		else if(libType==LFPE){sb.append("l");}
		else if(libType==CLRS){sb.append("s");}
		
		if(phixFlag){sb.append("p");}
		
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * TODO:  Some machines are set to UTC rather than PST
	 * @return Timestamp in RQC's format
	 */
	public static String timeString(){
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
//		sdf.setTimeZone(TimeZone.getTimeZone("PST"));
		sdf.setTimeZone(TimeZone.getDefault());
		return sdf.format(new Date());
	}
	
	/**
	 * Set permissions on these files to 777
	 * @param names List of filenames
	 */
	private static void setPermissions(String...names){
		if(names==null){return;}
		for(String name : names){
			if(name!=null && name.trim().length()>0 && new File(name).exists()){
				ReadWrite.setPermissions(name, true, true, true, false);
			}
		}
	}
	
	/**
	 * Write a string to the file containing steps needed to regenerate the output
	 * @param fname Filename to write, including path
	 * @param command Command to add to file
	 * @param args Arguments to the command
	 * @param append Append to existing file rather than overwriting
	 * @param overwrite Permission to overwrite
	 */
	private static void writeReproduceFile(String fname, String command, String[] args, boolean append, boolean overwrite){
		StringBuilder sb=new StringBuilder();
		if(!append){
			boolean b=Tools.canWrite(fname, overwrite);
			assert(b) : "Can't write to "+fname;
			sb.append("#!/bin/bash\n");
		}
		sb.append(command);
		if(args!=null){
			for(String s : args){
				sb.append(' ').append(s);
			}
		}
		sb.append('\n');
		ReadWrite.writeString(sb, fname, append);
	}
	
	/**
	 * @param s String representation of library type
	 * @return Numeric code for library type
	 */
	private static int toLibType(String s){
		if(s==null){return FRAG;}
		s=s.trim().toLowerCase();
		if(s.equals("lfpe")){return LFPE;}
		if(s.equals("clip")){return CLIP;}
		if(s.equals("clrs")){return CLRS;}
		if(s.equals("frag") || s.equals("fragment")){return FRAG;}
		throw new RuntimeException("Unknown library type "+s);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Symbols to insert in output filename to denote operations performed */
	private final String symbols;
	
	/** Type of library; controls processing methods and references to use */
	private int libType=FRAG;
	/** True for rna artifacts, false for dna artifacts */
	private boolean rnaFlag=false;
	/** True if phix should be filtered out */
	private boolean phixFlag=false;
	/** Unused */
	private String jointSeq=null;
	/** Toss reads shorter than this */
	private int minLen=25;
	/** Toss reads shorter than this fraction of initital length, after trimming */
	private float minLenFraction=0.333f;
	/** Trim bases at this quality or below */
	private byte trimq=6;
	/** Throw away reads below this average quality before trimming.  Default: 6 */
	private byte maq=5;
	/** Quality-trimming mode */
	private String qtrim="f";//"rl";
	/** Kmer-trimming mode */
	private String ktrim="r";
	/** Shortest kmer to use for trimming */
	private int mink=8;
	/** Throw away reads containing more than this many Ns.  Default: 0 (toss reads with any Ns) */
	private int maxNs=0;
	
	/** Trim fragment adapters from right side of reads */
	private boolean fragAdapterFlag=false;
	
	/** Remove reads mapping to human with high identity */
	private boolean humanFlag=false;
	/** Use indexed version of human reference, rather than regenerating from fasta */
	private boolean humanRefIndexedFlag=true;
	
	private boolean verbose=false;
	private boolean overwrite=true;
	private boolean compress=true;
	
	private boolean writeTempToTmpdir=false;
	
	/** Arguments to pass to BBDuk */
	private ArrayList<String> primaryArgList=new ArrayList<String>();
	/** References to pass to BBDuk for artifact removal */
	private ArrayList<String> refs=new ArrayList<String>();
	
	/*--------------------------------------------------------------*/
	/*----------------        Read Data Files       ----------------*/
	/*--------------------------------------------------------------*/

	private final String tempSalt;
	
	private final String trimPrefix;
	private final String humanPrefix;
	private final String filterPrefix;
	
	/** Directory in which to write all files */
	private String outDir="";
	
	/** Directory in which to write all temp files */
	private String tmpDir=Shared.TMPDIR;
	
	/** Primary input reads file (required) */
	private String in1=null;
	/** Secondary input reads file */
	private String in2=null;
	/** Primary output reads file (required) */
	private String out1=null;
	/** Secondary output reads file */
	private String out2=null;
	/** Primary input qual file */
	private String qfin1=null;
	/** Secondary input qual file */
	private String qfin2=null;
	/** Primary output qual file */
	private String qfout1=null;
	/** Secondary output qual file */
	private String qfout2=null;
	
	/*--------------------------------------------------------------*/
	/*----------------           Log Files          ----------------*/
	/*--------------------------------------------------------------*/
	
	private String logName="status.log";
	private String reproduceName="reproduce.sh";
	private String fileListName="file-list.txt";
	
	private String rqcStatsName="filterStats.txt";
	private String kmerStatsName="kmerStats.txt";
	private String scaffoldStatsName="scaffoldStats.txt";
	
	/** ktrim phase rqc stats file */
	private String rqcStatsName_kt;
	/** ktrim phase stats file */
	private String kmerStatsName_kt;
	/** ktrim phase scaffold stats file */
	private String scaffoldStatsName_kt;
	
	/*--------------------------------------------------------------*/
	/*----------------        Reference Files       ----------------*/
	/*--------------------------------------------------------------*/
	
	private String mainArtifactFile = "/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/Illumina.artifacts.2013.12.no_DNA_RNA_spikeins.fa";
	private String artifactFileRna = "/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/RNA_spikeins.artifacts.2012.10.NoPolyA.fa";
	private String artifactFileDna = "/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/DNA_spikeins.artifacts.2012.10.fa";
	private String phixRef = "/global/dna/shared/rqc/ref_databases/qaqc/databases/phix174_ill.ref.fa";
	private String lfpeLinker = "/global/dna/shared/rqc/ref_databases/qaqc/databases/lfpe.linker.fa";
	private String clrsLinker = "/global/dna/shared/rqc/ref_databases/qaqc/databases/crelox.fa";
	private String clipLinker = clipLinkerDefault; //A literal string; "CATG" is supposed to be the normal linker.
	
	private String allArtifactsLatest = "/global/projectb/sandbox/rqc/qcdb/illumina.artifacts/Illumina.artifacts.fa";
	private String fragAdapter = "/global/projectb/sandbox/gaag/bbtools/data/adapters.fa";
	private String humanPath = "/global/projectb/sandbox/gaag/bbtools/hg19/";
	private String humanRef = "/global/projectb/sandbox/gaag/bbtools/hg19/hg19.fa.gz";
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Library type codes */
	private static final int FRAG=0, LFPE=1, CLIP=2, CLRS=3;
	private static final String clipLinkerDefault = "CATG";
	
}
