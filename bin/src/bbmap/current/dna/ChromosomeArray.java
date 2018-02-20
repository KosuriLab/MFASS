package dna;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import align2.Tools;

import fileIO.ReadWrite;


public class ChromosomeArray implements Serializable {
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 3199182397853127842L;

	public static void main(String[] args){
		translateFile(args[1], Byte.parseByte(args[0]));
	}
	
	
	private static void translateFile(String fname, int chrom){
		
		long time1=System.nanoTime();
		
		ChromosomeArray cha=read(fname, chrom);
		cha.chromosome=chrom;
		long time2=System.nanoTime();
		
		int dot=fname.lastIndexOf(".fa");
		String outfile=fname.substring(0,dot).replace("hs_ref_", "")+".chrom";
		
		System.out.println("Writing to "+outfile);
		
		System.out.println("minIndex="+cha.minIndex+", maxIndex="+cha.maxIndex+", length="+cha.array.length+
				"; time="+String.format("%.3f seconds", (time2-time1)/1000000000d));

		long time3=System.nanoTime();
		ReadWrite.write(cha, outfile, false);
		cha=null;
		System.gc();
		cha=read(outfile);
		long time4=System.nanoTime();
		
		System.out.println("minIndex="+cha.minIndex+", maxIndex="+cha.maxIndex+", length="+cha.array.length+
				"; time="+String.format("%.3f seconds", (time4-time3)/1000000000d));
	}
	
	public static ChromosomeArray read(String fname, int chrom){
		ChromosomeArray cha=read(fname);
		assert(cha.chromosome<1);
		cha.chromosome=chrom;
		return cha;
	}
	
	public static ChromosomeArray read(String fname){
		
		if(fname.endsWith(".chrom") || fname.endsWith(".chrom.gz")){
			ChromosomeArray ca=ReadWrite.read(ChromosomeArray.class, fname);
			return ca;
		}else{
			assert(fname.endsWith(".chromC") || fname.endsWith(".chromC.gz"));
			
			ChromosomeArrayCompressed cac=ReadWrite.read(ChromosomeArrayCompressed.class, fname);
			return cac.toChromosomeArray();
		}
	}
	
	public ChromosomeArray(){
		this((byte)-1, Gene.PLUS);
	}
	
	public ChromosomeArray toColorspace(){
		assert(!colorspace);
		ChromosomeArray ca=new ChromosomeArray(chromosome, strand, 0, maxIndex, true);
		
		for(int i=0; i<maxIndex; i++){
			byte a=AminoAcid.baseToColor(get(i), get(i+1));
			ca.set(i, a);
		}
		
//		System.err.println(maxIndex+", "+ca.maxIndex+", "+ca.array.length);
		
//		assert(ca.maxIndex==maxIndex-1);
//		assert(ca.array.length==ca.maxIndex+1);
		
		return ca;
	}
	
	/** Actually does reverse complement */
	public ChromosomeArray complement(){
		assert(!colorspace) : "Needs different method";
		byte otherStrand=(strand==Gene.MINUS ? Gene.PLUS : Gene.MINUS);
		ChromosomeArray ca=new ChromosomeArray(chromosome, otherStrand, 0, maxIndex);
		for(int i=0; i<=maxIndex; i++){
			int pos=maxIndex-i;
			byte b=AminoAcid.baseToComplementExtended[array[i]];
			ca.array[pos]=b;
		}
		return ca;
	}
	
	public ChromosomeArray(int chrom, byte strnd){
		this(chrom, strnd, false);
	}
	
	public ChromosomeArray(int chrom, byte strnd, boolean cs){
		chromosome=chrom;
		strand=strnd;
		array=new byte[1<<25];
		colorspace=cs;
	}
	
	public ChromosomeArray(int chrom, byte strnd, int min, int max){
		this(chrom, strnd, min, max, false);
	}
	
	public ChromosomeArray(int chrom, byte strnd, int min, int max, boolean cs){
		chromosome=chrom;
		strand=strnd;
		array=new byte[max+1];
		minIndex=min;
		maxIndex=max;
		colorspace=cs;
	}
	
	public ChromosomeArray(int chrom, byte strnd, String s){
		this(chrom, strnd, s, false);
	}
	
	public ChromosomeArray(int chrom, byte strnd, String s, boolean cs){
		chromosome=chrom;
		strand=strnd;
		array=s.getBytes();
		minIndex=0;
		maxIndex=s.length()-1;
		colorspace=cs;
	}
	
	
	public void set(int loc, int val){
		
		if(loc>=array.length){//Increase size
			int newlen=(int)(1+(3L*max(array.length, loc))/2);
			assert(newlen>loc) : newlen+", "+loc+", "+array.length;
			resize(newlen);
			assert(array.length==newlen);
//			System.err.println("Resized array to "+newlen);
		}
		char c=Character.toUpperCase((char)val);
		if(AminoAcid.baseToNumberExtended[c]<0){c='N';}
		array[loc]=(val>Byte.MAX_VALUE ? Byte.MAX_VALUE : (byte)val);
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}
	
	
	public void set(int loc, CharSequence s){
		int loc2=loc+s.length();
		if(loc2>array.length){//Increase size
			int newlen=(int)(1+(3L*max(array.length, loc2))/2);
			assert(newlen>loc2) : newlen+", "+loc2+", "+array.length;
			resize(newlen);
			assert(array.length==newlen);
//			System.err.println("Resized array to "+newlen);
		}
		
		for(int i=0; i<s.length(); i++, loc++){
			char c=Character.toUpperCase(s.charAt(i));
			if(AminoAcid.baseToNumberExtended[c]<0){c='N';}
			assert(Character.isLetter(c));
			assert(c<=Byte.MAX_VALUE);
			array[loc]=(byte)c;
		}
		loc--;
		assert(loc==loc2-1) : "loc="+loc+", loc2="+loc2+", s.len="+s.length();
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}
	
	/** Returns the letter (IUPAC) representation of the base, as a byte */
	public byte get(int loc){
		return loc<minIndex || loc>=maxIndex ? (byte)'N' : array[loc];
	}
	
	public String getString(int a, int b){
		StringBuilder sb=new StringBuilder(b-a+1);
		for(int i=a; i<=b; i++){
			sb.append((char)get(i));
		}
		return sb.toString();
	}
	
	/** Returns FASTA format bytes.  Same as getString, but faster. */
	public byte[] getBytes(int a, int b){
		byte[] out=Arrays.copyOfRange(array, a, b+1);
//		assert(out[0]>0 && out[out.length-1]>0) : a+", "+b+", "+minIndex+", "+maxIndex+", "+array.length;
		if(a<minIndex || b>maxIndex){
			for(int i=0; i<out.length; i++){
				if(out[i]==0){out[i]='N';}
			}
		}
		return out;
	}
	
	public byte getNumberACGTN(int loc){
		return AminoAcid.baseToNumberACGTN[array[loc]];
	}
	
	public byte getNumber(int loc){
		return AminoAcid.baseToNumber[array[loc]];
	}
	
	public byte getNumber(int loc, boolean colorspace){
		final byte b=array[loc];
		if(colorspace){
			return b>3 ? -1 : b;
		}else{
			return AminoAcid.baseToNumber[b];
		}
	}
	
	public boolean isFullyDefined(int a, int b){
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[array[i]];
			if(x<0){return false;}
		}
		return true;
	}
	
	public boolean isFullyUndefined(int a, int b){
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[array[i]];
			if(x>=0){return false;}
		}
		return true;
	}
	
	public int countDefinedBases(){
		return countDefinedBases(minIndex, maxIndex);
	}
	
	public int countDefinedBases(int a, int b){
		int sum=0;
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[array[i]];
			if(x>=0){sum++;}
		}
		return sum;
	}
	
	public int getNumber(int a, int b){
		return toNumber(a, b, array);
	}
	
	public int getNumber(int a, int b, boolean colorspace){
		return colorspace ? toNumberColorspace(a, b, array) : toNumber(a, b, array);
	}
	
	public static int toNumberColorspace(int a, int b, byte[] bases){
		assert(b>=a);
		assert(b-a<17); //<17 for unsigned, <16 for signed
		int out=0;
		for(int i=a; i<=b; i++){
			int x=bases[i];
			if(x<0 || x>3){return -1;}
			out=((out<<2)|x);
		}
		return out;
	}
	
	public static int toNumber(int a, int b, byte[] bases){
		assert(b>=a);
		assert(b-a<17); //<17 for unsigned, <16 for signed
		int out=0;
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[bases[i]];
			if(x<0){return -1;}
			out=((out<<2)|x);
		}
		return out;
	}
	
	public static int toNumber(int a, int b, String bases){
		int out=0;
		for(int i=a; i<=b; i++){
			int x=AminoAcid.baseToNumber[bases.charAt(i)];
			if(x<0){return -1;}
			out=((out<<2)|x);
		}
		return out;
	}
	
	public void resize(int newlen){
		byte[] temp=new byte[newlen];
		int lim=min(array.length, newlen);
		assert(lim>=maxIndex) : lim+","+maxIndex;
		for(int i=0; i<lim; i++){
			temp[i]=array[i];
		}
		array=temp;
	}
	
	public String toBaseString(){
		String s=new String(array);
		return s;
	}
	
	public char[] nearestDefinedBase(){
		char[] r=new char[array.length];
		final char max=Character.MAX_VALUE;
		
		char dist=max;
		for(int i=0; i<r.length; i++){
			byte b=array[i];
			if(b=='A' || b=='C' || b=='G' || b=='T'){
				dist=0;
			}else{
				dist=(dist==max ? max : (char)(dist+1));
			}
			r[i]=dist;
		}
		
		dist=r[r.length-1];
		for(int i=r.length-1; i>=0; i--){
			byte b=array[i];
			if(b=='A' || b=='C' || b=='G' || b=='T'){
				dist=0;
			}else{
				dist=(dist==max ? max : (char)(dist+1));
			}
			r[i]=Tools.min(dist, r[i]);
		}
		return r;
	}
	
	public ArrayList<Range> toContigRanges(final int nBlockSize){
		assert(nBlockSize>0);
		ArrayList<Range> list=new ArrayList<Range>();
		
		int start=-1;
		int stop=-1;
		int ns=nBlockSize+1;
		
		boolean contig=false;
		
		for(int i=minIndex; i<=maxIndex; i++){
			byte b=array[i];
			if(b=='N' || b=='X'){
				ns++;
				if(contig && (b=='X' || ns>=nBlockSize)){
					Range r=new Range(start, stop);
					list.add(r);
					contig=false;
				}
			}else{
				ns=0;
				if(!contig){start=i;}
				contig=true;
				stop=i;
			}
		}
		if(contig){
			Range r=new Range(start, stop);
			list.add(r);
		}
		return list;
	}
	
	
	public boolean equalsIgnoreCase(ChromosomeArray other){
		if(minIndex!=other.minIndex){System.err.println("a");return false;}
		if(maxIndex!=other.maxIndex){System.err.println("b");return false;}
		if(chromosome!=other.chromosome){System.err.println("c");return false;}
		if(array.length!=other.array.length){System.err.println("d");return false;}
		for(int i=minIndex; i<=maxIndex; i++){
			if(Character.toLowerCase(array[i])!=Character.toLowerCase(other.array[i])){
				System.err.println("e");
				return false;
			}
		}
		return true;
	}
	
	private static final long min(long x, long y){return x<y ? x : y;}
	private static final long max(long x, long y){return x>y ? x : y;}
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public final byte strand;
	public int chromosome;
	public byte[] array;
	public int maxIndex=-1;
	public int minIndex=Integer.MAX_VALUE;
	
	public final boolean colorspace;
	
	
}
