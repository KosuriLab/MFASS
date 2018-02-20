package align2;

import java.lang.management.ManagementFactory;
import java.util.List;

import dna.Data;

public class Shared {
	
	public static int THREADS=SET_THREADS(-1);

	public static int READ_BUFFER_LENGTH=200;
	public static int READ_BUFFER_NUM_BUFFERS=Tools.max(4, (THREADS*3)/2);
	public static final long READ_BUFFER_MAX_DATA=500000;
	
	//TODO:  Actually... for some reason...  it seems as though GAPBUFFER must equal exactly 1/2 of GAPLEN.  Not good; 1/4 would be far better.
	
	public static final int GAPBUFFER=64; //TODO:  Seems to break less than 64, for some reason
	public static final int GAPBUFFER2=2*GAPBUFFER;
	public static final int GAPLEN=128; //TODO: May break when over 128
	public static final int MINGAP=GAPBUFFER2+GAPLEN;
	public static final int GAPCOST=Tools.max(1, GAPLEN/64);
	public static final byte GAPC='-';

	public static int BBMAP_VERSION=31;
	public static int BBMAP_VERSION_MINOR=27;
	public static String BBMAP_VERSION_STRING=BBMAP_VERSION+"."+BBMAP_VERSION_MINOR;

	public static boolean TRIM_READ_COMMENTS=false;

	public static String BBMAP_CLASS=null;
	public static String[] COMMAND_LINE=null;
	public static List<String> JVM_ARGS(){
		return ManagementFactory.getRuntimeMXBean().getInputArguments();
	}

	/** Directory in which to write temp files */
	public static String TMPDIR=(System.getenv("TMPDIR")==null ? null : (System.getenv("TMPDIR")+"/").replaceAll("//", "/"));
//	static{assert(false) : "TMPDIR="+TMPDIR;}
	
	/** Anomaly probably resolved as of v.20.1 
	 * This variable should be TRUE for normal users and FALSE for me. */
	public static boolean anomaly=!System.getProperty("user.dir").contains("/bushnell/") && !Data.WINDOWS;
	
	public static final char[] getTLCB(int len){
		char[] buffer=TLCB.get();
		if(buffer==null || buffer.length<len){
			buffer=new char[len];
			if(len<1000000){TLCB.set(buffer);}
		}
		return buffer;
	}
	private static final ThreadLocal<char[]> TLCB=new ThreadLocal<char[]>();

	public static int SET_THREADS(int x){
		if(x>0){
			THREADS=x;
		}else{
			THREADS=(Data.HOSTNAME()==null || !Data.HOSTNAME().startsWith("gpint") ? Data.LOGICAL_PROCESSORS : Tools.min(4, Data.LOGICAL_PROCESSORS));
		}
//		assert(false) : Data.HOSTNAME()+", "+THREADS;
		return THREADS;
	}
	
}
