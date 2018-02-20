package dna;
import java.io.Serializable;

import driver.Translator2;

import fileIO.ReadWrite;


public class CoverageArray1 extends CoverageArray implements Serializable {
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 6711045833114428632L;
	
	
	public static void main(String[] args){
		runSpeedTest(args);
		
//		translateGenomeBuild(args);
	}
	
	public static void runSpeedTest(String[] args){
		
		long time1=System.nanoTime();
		
		CoverageArray1 ca=(CoverageArray1)read(args[1]);
		ca.chromosome=Byte.parseByte(args[0]);
		long time2=System.nanoTime();
		
//		int dot=args[1].lastIndexOf(".");
//		String outfile=args[1].substring(0,dot)+".ca";
		
		args[1]=args[1].replace('\\', '/');
		int slash=args[1].lastIndexOf('/');
		String outfile;
		if(slash<1){
			outfile="coverage-chr"+Gene.chromCodes[ca.chromosome]+"-build"+Data.GENOME_BUILD+".ca";
		}else{
			outfile=args[1].substring(0,slash+1)+"coverage-chr"+Gene.chromCodes[ca.chromosome]+"-build"+Data.GENOME_BUILD+".ca";
		}
		
		System.out.println("minIndex="+ca.minIndex+", maxIndex="+ca.maxIndex+", length="+ca.array.length+
				"; time="+String.format("%.3f seconds", (time2-time1)/1000000000d));

		long time3=System.nanoTime();
		ReadWrite.write(ca, outfile, false);
		ca=null;
		System.gc();
		ca=(CoverageArray1)read(outfile);
		long time4=System.nanoTime();
		
		System.out.println("minIndex="+ca.minIndex+", maxIndex="+ca.maxIndex+", length="+ca.array.length+
				"; time="+String.format("%.3f seconds", (time4-time3)/1000000000d));
		
		
	}
	
	public static void translateGenomeBuild(String[] args){

		Timer t=new Timer();
		t.start();

		int inBuild=Integer.parseInt(args[0]);
		int outBuild=Integer.parseInt(args[1]);
		String root=args[2];
		
		translateGenomeBuild(inBuild, outBuild, root);
		
		t.stop();
		System.out.println("Time:\t"+t);
		
	}
	
	public static void translateGenomeBuild(int inBuild, int outBuild, String root){
		root=root.replace('\\', '/');
		if(!root.endsWith("/")){root+="/";}
		
		CoverageArray1[] out=new CoverageArray1[27];
		
		for(int chrom=1; chrom<out.length; chrom++){
			out[chrom]=new CoverageArray1(chrom);
		}
		
		final byte PLUS=Gene.PLUS;
		
		for(int chrom=1; chrom<=25; chrom++){
			String infile=root+"coverage-chr"+Gene.chromCodes[chrom]+"-build"+inBuild+".ca.zip";
			CoverageArray1 ca1=ReadWrite.read(CoverageArray1.class, infile);
			for(int loc1=ca1.minIndex; loc1<=ca1.maxIndex; loc1++){
				short cov=(short)ca1.get(loc1);
				int[] xform=Translator2.translate(inBuild, outBuild, chrom, PLUS, loc1);
				if(xform!=null){
					int chrom2=(int)xform[0];
					int loc2=xform[2];
					out[chrom2].set(loc2, cov);
				}
			}
			ca1=null;
			System.out.println("Read "+infile);
		}
		
		for(int chrom=1; chrom<=25; chrom++){
			String outfile=root+"coverage-chr"+Gene.chromCodes[chrom]+"-build"+outBuild+".ca.zip";
			out[chrom].resize(out[chrom].maxIndex+1);
			ReadWrite.write(out[chrom], outfile, false);
			out[chrom]=null;
			System.out.println("Wrote "+outfile);
		}
		
	}
	
	public CoverageArray1(){
		this((int)-1);
	}
	
	public CoverageArray1(int chrom){
		this(chrom, 1<<24);
	}
	
	public CoverageArray1(int chrom, int initialLen){
		super(chrom);
		array=new short[initialLen];
	}
	
	/**
	 * @param loc
	 * @param amt
	 */
	public void increment(int loc, int amt) {
		set(loc, get(loc)+amt);
	}
	
	/**
	 * @param loc
	 */
	public void increment(int loc) {
		set(loc, get(loc)+1);
	}

	public void incrementRange(int min, int max, int amt) {
		if(min<0){min=0;}
		if(max>=array.length){//Increase size
			int newlen=1+(7*max(array.length, max))/4;
			assert(newlen>max);
			resize(newlen);
			assert(array.length==newlen);
		}else if(max<0){max=-1;}
		for(int i=min; i<=max; i++){
			int val=array[i]+amt;
			if(val>Short.MAX_VALUE){
				val=Short.MAX_VALUE;
				 if(!OVERFLOWED){
					 System.err.println("Note: Coverage capped at "+Short.MAX_VALUE);
					 OVERFLOWED=true;
				 }
			}
			array[i]=(short)val;
		}
	}
	
	public void set(int loc, int val){
		
		if(loc>=array.length){//Increase size
			int newlen=1+(7*max(array.length, loc))/4;
			assert(newlen>loc);
			resize(newlen);
			assert(array.length==newlen);
		}else if(loc<0){
//			minIndex=min(0, minIndex);
//			maxIndex=max(0, maxIndex);
			return;
		}
		
		if(val>Short.MAX_VALUE && !OVERFLOWED){
			System.err.println("Note: Coverage capped at "+Short.MAX_VALUE);
			OVERFLOWED=true;
		}
		array[loc]=(val>Short.MAX_VALUE ? Short.MAX_VALUE : (short)val);
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}
	
	public int get(int loc){
		return loc>=array.length || loc<0 ? 0 : array[loc];
	}
	
	public void resize(int newlen){
		System.out.println("Resized chrom "+chromosome+" to "+newlen);
		short[] temp=new short[newlen];
		int lim=min(array.length, newlen);
		assert(lim>maxIndex) : lim+","+maxIndex;
		for(int i=0; i<lim; i++){
			temp[i]=array[i];
		}
		array=temp;
	}
	
	
	private static final long min(long x, long y){return x<y ? x : y;}
	private static final long max(long x, long y){return x>y ? x : y;}
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public short[] array;
	public int length(){return maxIndex-minIndex+1;}
	public int arrayLength(){return array.length;}
	
	private static boolean OVERFLOWED=false;
	/**
	 * 
	 */
//	private static final long serialVersionUID = -7493066925636540386L;
	
}
