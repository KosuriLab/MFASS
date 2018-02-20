package fileIO;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

/**
 * @author Brian Bushnell
 * @date Dec 19, 2012
 *
 */
public final class FileFormat {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static FileFormat testInput(String fname, String overrideExtension, boolean allowSubprocess){
		return testInput(fname, FASTQ, overrideExtension, allowSubprocess, true);
	}
	
	public static FileFormat testInput(String fname, int defaultFormat, String overrideExtension, boolean allowSubprocess, boolean allowFileRead){
		if(fname==null){return null;}
		int overrideFormat=0;
		int overrideCompression=0;
		if(overrideExtension!=null && overrideExtension.length()>0){
			int[] a=testFormat(overrideExtension, false);
			if(a!=null){
				overrideFormat=a[0];
				if(a[1]!=RAW){overrideCompression=a[1];}
			}
		}
		return testInput(fname, defaultFormat, overrideFormat, overrideCompression, allowFileRead, allowSubprocess);
	}

	public static FileFormat testInput(String fname, int defaultFormat, int overrideFormat, int overrideCompression, boolean allowSubprocess, boolean allowFileRead){
		if(fname==null){return null;}
		return new FileFormat(fname, READ, defaultFormat, overrideFormat, overrideCompression, allowFileRead, false, allowSubprocess, false);
	}
	
	public static FileFormat testOutput(String fname, int defaultFormat, String overrideExtension, boolean allowSubprocess, boolean overwrite, boolean ordered){
		if(fname==null){return null;}
		int overrideFormat=0;
		int overrideCompression=0;
		if(overrideExtension!=null && overrideExtension.length()>0){
			int[] a=testFormat(overrideExtension, false);
			if(a!=null){
				overrideFormat=a[0];
				if(a[1]!=RAW){overrideCompression=a[1];}
			}
		}
		return testOutput(fname, defaultFormat, overrideFormat, overrideCompression, allowSubprocess, overwrite, ordered);
	}
	
	public static FileFormat testOutput(String fname, int defaultFormat, int overrideFormat, int overrideCompression, boolean allowSubprocess, boolean overwrite, boolean ordered){
		if(fname==null){return null;}
		return new FileFormat(fname, WRITE, defaultFormat, overrideFormat, overrideCompression, false, overwrite, allowSubprocess, ordered);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Constructor         ----------------*/
	/*--------------------------------------------------------------*/
	
	private FileFormat(String fname, int mode_, int defaultFormat, int overrideFormat, int overrideCompression, boolean allowFileRead
			, boolean overwrite_, boolean allowSubprocess_, boolean ordered_){
//			, boolean interleaved_, boolean colorspace_, long maxReads_){
		
		if(verbose){
			new Exception().printStackTrace(System.err);
			System.err.println("FileFormat(fname="+fname+", mode="+mode_+", dFormat="+defaultFormat+", oFormat="+overrideFormat+", oCompression="+overrideCompression+
					", allowRead="+allowFileRead+", ow="+overwrite_+", allowSub="+allowSubprocess_+", ordered="+ordered_+")");
		}
		
		assert(fname!=null);
		fname=fname.trim().replace('\\', '/');
		assert(fname.trim().length()>0) : fname;
		
		if(defaultFormat<1){defaultFormat=FQ;}
		allowFileRead&=(mode_==READ);
		int[] a=testFormat(fname, allowFileRead);
		
		if(verbose){System.err.println(Arrays.toString(a));}
		
		if(a[0]==UNKNOWN && overrideFormat<1){
			a[0]=defaultFormat;
			if(defaultFormat!=TEXT){
				System.err.println("Unspecified format for "+(mode_==READ ? "input" : "output")+" "+(fname==null ? "stream" : fname)+"; defaulting to "+FORMAT_ARRAY[a[0]]+".");
			}
		}
		if(verbose){System.err.println(Arrays.toString(a));}
		
		if(overrideFormat>0){a[0]=overrideFormat;}
		if(overrideCompression>0){a[1]=overrideCompression;}
		
		if(verbose){System.err.println(Arrays.toString(a));}
		
		name=fname;
		format=a[0];
		compression=a[1];
		type=a[2];
		mode=mode_;
		
		overwrite=overwrite_;
		allowSubprocess=allowSubprocess_;
		ordered=ordered_;
		
//		interleaved=interleaved_;
//		colorspace=colorspace_;
//		maxReads=write() ? -1 : maxReads_;

		assert(!unknownFormat()) : "Unknown file format for "+fname+"\n"+
			mode_+", "+defaultFormat+", "+overrideFormat+", "+overrideCompression+", "+allowFileRead+", "+overwrite_+", "+allowSubprocess_;
		assert(!unknownCompression()) : "Unknown compression for "+fname+"\n"+
			mode_+", "+defaultFormat+", "+overrideFormat+", "+overrideCompression+", "+allowFileRead+", "+overwrite_+", "+allowSubprocess_;
		assert(!unknownType()) : "Unknown stream type for "+fname+"\n"+
			mode_+", "+defaultFormat+", "+overrideFormat+", "+overrideCompression+", "+allowFileRead+", "+overwrite_+", "+allowSubprocess_;
		assert(!unknownMode()) : "Unknown I/O mode for "+fname+"\n"+
			mode_+", "+defaultFormat+", "+overrideFormat+", "+overrideCompression+", "+allowFileRead+", "+overwrite_+", "+allowSubprocess_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append(name).append(',');
		sb.append(format+"("+FORMAT_ARRAY[format]+")").append(',');
		sb.append(compression+"("+COMPRESSION_ARRAY[compression]+")").append(',');
		sb.append(type+"("+TYPE_ARRAY[type]+")").append(',');
		sb.append(mode+"("+MODE_ARRAY[mode]+")").append(',');
		sb.append("ow="+(overwrite ? "t" : "f")).append(',');
		sb.append("sub="+(allowSubprocess ? "t" : "f")).append(',');
		sb.append("ordered="+(ordered ? "t" : "f"));
		return sb.toString();
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static int[] testFormat(String fname, boolean allowFileRead){
		int[] r=new int[] {UNKNOWN, RAW, FILE};
		if(fname==null || fname.length()<1){
			r[2]=STDIO;
			return r;
		}
		String slc=fname.trim().toLowerCase();
		if(slc.indexOf('/')<0){slc=slc.substring(slc.lastIndexOf('/')+1);}
		if(slc.indexOf('.')<0){slc="."+slc;}
		String comp=ReadWrite.compressionType(slc);
		String ext=ReadWrite.rawExtension(slc);
		
		if(ext==null){}
		else if(ext.equals("fq") || ext.equals("fastq")){r[0]=FASTQ;}
		else if(ext.equals("fa") || ext.equals("fasta") || ext.equals("fas") || ext.equals("fna") || ext.equals("ffn") 
				|| ext.equals("frn") || ext.equals("seq")|| ext.equals("fsa")){r[0]=FASTA;}
		else if(/*ext.equals("txt") || */ext.equals("bread")){r[0]=BREAD;}
		else if(ext.equals("sam")){r[0]=SAM;}
		else if(ext.equals("csfasta")){r[0]=CSFASTA;}
		else if(ext.equals("qual")){r[0]=QUAL;}
		else if(ext.equals("bam")){r[0]=BAM;}
		else if(ext.equals("sites") || ext.equals("sitesonly")){r[0]=SITES;}
		else if(ext.equals("info") || ext.equals("attachment")){r[0]=ATTACHMENT;}
		else if(ext.equals("scarf")){r[0]=SCARF;}
		
		if(comp==null){}
		else if(comp.equals("gz")){r[1]=GZ;}
		else if(comp.equals("zip")){r[1]=ZIP;}
		else if(comp.equals("bz2")){r[1]=BZ2;}
		else if(comp.equals("xz")){r[1]=XZ;}
		
//		assert(false) : Arrays.toString(r);
		

		if(slc.length()>2 && slc.charAt(0)=='s' && slc.charAt(1)=='t'){
			if(slc.equals("stdin") || slc.startsWith("stdin.") || slc.equals("standardin")){r[2]=STDIO;}
			else if(slc.equals("stdout") || slc.startsWith("stdout.") || slc.equals("standardout")){r[2]=STDIO;}
		}else if("/dev/null".equalsIgnoreCase(slc)){
			r[2]=DEVNULL;
		}
		
		if(r[0]==UNKNOWN){
			File f=(allowFileRead && r[2]==FILE ? new File(fname) : null);
			if(f!=null && f.exists() && !f.isDirectory()){
				InputStream is=ReadWrite.getInputStream(fname, false, false);
				int b=-1;
				try {
					b=is.read();
					is.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(b=='>'){r[0]=FA;}
				else if(b=='@'){r[0]=FQ;} //TODO: Note - could be sam
				else{r[0]=BREAD;}
			}else{
				if(fname.equals("sequential")){r[0]=SEQUENTIAL;}
				else if(fname.equals("random")){r[0]=RANDOM;}
				else if(fname.equals("sitesonly")){r[0]=SITES;}
			}
		}
		
		
		if(r[2]==STDIO && allowFileRead){
			File f=new File(fname);
			if(f.exists() && !f.isDirectory()){r[2]=FILE;}
		}
//		else{
//			r[2]=FILE; //What is this for?
//		}
		
		return r;
	}
	
	public static boolean hasFastaExtension(String fname){
		int[] r=testFormat(fname, false);
		return r[0]==FA;
	}
	
	public static boolean hasFastqExtension(String fname){
		int[] r=testFormat(fname, false);
		return r[0]==FQ;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	public final String name(){return name;}
	public final int format(){return format;}
	public final int compression(){return compression;}
	public final int type(){return type;}
	public final int mode(){return mode;}

	public final boolean hasName(){return name!=null;}
	public final boolean canWrite(){
		assert(write());
		if(stdio() || devnull()){return true;}
		assert(hasName());
		File f=new File(name);
		if(!f.exists()){return true;}
		if(!f.canWrite()){return false;}
		return overwrite();
	}
	public final boolean canRead(){
		assert(read());
		if(stdio()){return true;}
		assert(hasName());
		File f=new File(name);
		return f.canRead();
	}
	
	public final boolean unknownField(){return unknownFormat() || unknownCompression() || unknownType() || unknownMode();}

	public final boolean unknownFormat(){return format<=UNKNOWN;}
	public final boolean fasta(){return format==FASTA;}
	public final boolean fastq(){return format==FASTQ;}
	public final boolean bread(){return format==BREAD;}
	public final boolean sam(){return format==SAM;}
	public final boolean samOrBam(){return format==SAM || format==BAM;}
	public final boolean csfasta(){return format==CSFASTA;}
	public final boolean qual(){return format==QUAL;}
	public final boolean sequential(){return format==SEQUENTIAL;}
	public final boolean random(){return format==RANDOM;}
	public final boolean sites(){return format==SITES;}
	public final boolean attachment(){return format==ATTACHMENT;}
	public final boolean bam(){return format==BAM;}
	public final boolean scarf(){return format==SCARF;}
	public final boolean text(){return format==TEXT;}

	public final boolean unknownCompression(){return compression<=UNKNOWN;}
	public final boolean raw(){return compression==RAW;}
	public final boolean gzip(){return compression==GZIP;}
	public final boolean zip(){return compression==ZIP;}
	public final boolean bz2(){return compression==BZ2;}
	public final boolean xz(){return compression==XZ;}
	public final boolean sevenz(){return compression==SEVENZ;}

	public final boolean unknownType(){return type<=UNKNOWN;}
	public final boolean file(){return type==FILE;}
	public final boolean stdio(){return type==STDIO;}
	public final boolean devnull(){return type==DEVNULL;}

	public final boolean unknownMode(){return mode<=UNKNOWN;}
	public final boolean read(){return mode==READ;}
	public final boolean write(){return mode==WRITE;}

	public final boolean overwrite(){return overwrite;}
	public final boolean allowSubprocess(){return allowSubprocess;}
	public final boolean ordered(){return ordered;}
	
//	public final boolean interleaved(){return interleaved;}
//	public final boolean colorspace(){return colorspace;}
//	public final long maxReads(){return maxReads;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String name;
	private final int format;
	private final int compression;
	private final int type;
	private final int mode;

	private final boolean overwrite;
	private final boolean allowSubprocess;
	private final boolean ordered;
	
//	private final boolean interleaved;
//	private final boolean colorspace;
//	private final long maxReads;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean verbose=false;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final int UNKNOWN=0;
	
	/* Format */
	
	public static final int FA=1, FASTA=1;
	public static final int FQ=2, FASTQ=2;
	public static final int BREAD=3;
	public static final int SAM=4;
	public static final int CSFASTA=5;
	public static final int QUAL=6;
	public static final int SEQUENTIAL=7;
	public static final int RANDOM=8;
	public static final int SITES=9;
	public static final int ATTACHMENT=10;
	public static final int BAM=11;
	public static final int SCARF=12;
	public static final int TEXT=13;
	
	private static final String[] FORMAT_ARRAY=new String[] {
		"unknown", "fasta", "fastq", "bread", "sam", "csfasta",
		"qual", "sequential", "random", "sites", "attachment",
		"bam", "scarf", "text"
	};
	
	/* Compression */
	
	public static final int RAW=1;
	public static final int GZ=2, GZIP=2;
	public static final int ZIP=3;
	public static final int BZ2=4;
	public static final int XZ=5;
	public static final int c4=6;
	public static final int SEVENZ=7;
	
	private static final String[] COMPRESSION_ARRAY=new String[] {
		"unknown", "raw", "gz", "zip", "bz2", "xz",
		"c4", "7z"
	};
	
	/* Type */
	
	public static final int FILE=1;
	public static final int STDIO=2, STDIN=2, STDOUT=2;
	public static final int DEVNULL=3;
//	public static final int NULL=4;
	
	private static final String[] TYPE_ARRAY=new String[] {
		"unknown", "file", "stdio", "devnull"
	};
	
	/* Mode */
	
	public static final int READ=1, WRITE=2;
	
	private static final String[] MODE_ARRAY=new String[] {
		"unknown", "read", "write"
	};
	
}
