package dna;
import java.io.Serializable;

import align2.Tools;

import fileIO.ReadWrite;


public class ChromosomeArrayCompressed implements Serializable {
	
	
	public static void main(String[] args){
		
//		ChromosomeArray ca;
//		
//		
//		ca=Data.getChromosome(21);
//		ca=new ChromosomeArray((byte)1, "acagtgca");
//		
//		
//		ChromosomeArrayCompressed cac=new ChromosomeArrayCompressed(ca);
//		ChromosomeArray ca2=cac.toChromosomeArray();
//
//		assert(ca.minIndex==ca2.minIndex);
//		assert(ca.maxIndex==ca2.maxIndex);
//		assert(ca.chromosome==ca2.chromosome);
//		assert(ca.array.length==ca2.array.length);
//
//		System.out.println("Old: "+ca.toBaseString());
//		System.out.println("New: "+ca2.toBaseString());
//		
//		for(int i=0; i<=ca.maxIndex; i++){
//			if(Character.toLowerCase(ca.array[i])!=Character.toLowerCase(ca2.array[i])){
//				System.out.println("Error at "+i);
//				System.exit(1);
//			}
//		}

		if(args.length>2){
			Data.setGenome(Integer.parseInt(args[2]));
		}

		byte minChrom=1;
		byte maxChrom=26;
		
		String root=args[0].replace('\\', '/');
		if(!root.endsWith("/")){root+="/";}
		
		String fname=args[1];

		if(fname.contains("#")){
			for(int chrom=minChrom; chrom<=maxChrom; chrom++){
				
				try {
					translateFile(root, chrom, fname.replace("#", Gene.chromCodes[chrom]));
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}else{
//			translateFile(root, fname); //chrom is unknown
			assert(false);
		}
		
		if(args.length>2){
			for(int chrom=minChrom; chrom<=maxChrom; chrom++){
				System.out.println("Loading chr"+Gene.chromCodes[chrom]+" colorspace.");
				Data.getChromosome(chrom);
				Data.unload(chrom, true);
			}
		}
		
	}
	
	
	private static void translateFile(String root, int chrom, String fname){
		System.out.print(chrom+":\t");
		
		ChromosomeArray ca=ChromosomeArray.read(root+fname, (byte) chrom);
		System.out.print("Loaded\t");
//		System.out.println("\n"+ca.minIndex+", "+ca.maxIndex+"\n"+ca.getString(19999992,20000081)+"\n");
//		System.out.print(root+fname);
		
		assert(chrom==ca.chromosome);
		
		ChromosomeArrayCompressed cac=new ChromosomeArrayCompressed(ca);
		System.out.print("Translated\t");
		assert(chrom==cac.chromosome);
		ReadWrite.write(cac, root+"chr"+Gene.chromCodes[chrom]+".chromC", false);
		System.out.print("Wrote\t");
		
		cac=null;
		cac=ReadWrite.read(ChromosomeArrayCompressed.class, root+"chr"+Gene.chromCodes[chrom]+".chromC");
		System.out.print("Reloaded\t");
		assert(chrom==cac.chromosome);
		
		ChromosomeArray ca2=cac.toChromosomeArray();

//		System.out.println("\n"+ca2.minIndex+", "+ca2.maxIndex+"\n"+ca2.getString(19999992,20000081)+"\n");
		
		boolean success=ca2.equalsIgnoreCase(ca);
		
		System.out.println(success ? "Success" : "Fail");
		
		assert(success) : "\n"+ca2.getString(0, 100)+"\n"+ca.getString(0, 100);
		
//		Data.unload(chrom, true); //Why is this here?
	}
	
	
	public ChromosomeArrayCompressed(ChromosomeArray ca){
		this(ca.chromosome, ca.minIndex, ca.maxIndex, ca.array);
	}
	
	
	public ChromosomeArrayCompressed(int chrom, int min, int max, byte[] letters){
		
		
		array=new byte[Tools.max(max, letters.length-1)/2+1];
		minIndex=min;
		maxIndex=max;
		chromosome=chrom;
		
		for(int i=min; i<=max; i++){
			byte letter=letters[i];
			write(i, letter);
		}
		
	}
	
	
	public void translate(byte[] dest){
		
//		int min=minIndex/2;
//		int max=
		
		byte[] map=AminoAcid.numberToBaseExtended2;
		
		int max=dest.length/2;
		
		
//		int min=minIndex/2;
//		int min2=min*2;
		
		
//		for(int i=0, j=0, k=1; i<max; i++, j+=2, k+=2){
//			int b=array[i]&0xFF;
//			
//			byte n0=(byte)(b>>4);
//			dest[j]=map[n0];
//			
//			byte n1=(byte)(b&0xF);
//			dest[k]=map[n1];
//		}
		
		
		for(int i=0, j=0; i<max; i++, j+=2){
			int b=array[i]&0xFF;
			
			byte n0=(byte)(b>>4);
			dest[j]=map[n0];
			
			byte n1=(byte)(b&0xF);
			dest[j+1]=map[n1];
		}
		
		if((dest.length&1) == 1){
			dest[dest.length-1]=readLetter(dest.length-1);
		}
		
	}
	
	
	public ChromosomeArray toChromosomeArray(){
		
		//TODO Store strand data
		ChromosomeArray ca=new ChromosomeArray(chromosome, Gene.PLUS, minIndex, maxIndex);
		
//		byte[] caarray=ca.array;
//		
//		for(int i=ca.minIndex; i<=ca.maxIndex; i++){
//			caarray[i]=readLetter(i);
//		}
		
		translate(ca.array);
		
		return ca;
	}
	
	public byte readNumber(int pos){
		int remap=pos/2;
		byte bit1=(byte) (pos&1);
		int old=array[remap]&0xFF;
		
		if(bit1==0){
			old=(old>>>4);
		}else{
			old=(old&0x0F);
		}
		return (byte)old;
	}
	
	public byte readLetter(int pos){
		int remap=pos/2;
		byte bit1=(byte) (pos&1);
		int old=array[remap]&0xFF;
		
		if(bit1==0){
			old=(old>>>4);
		}else{
			old=(old&0x0F);
		}
		
//		System.out.println("Reading "+Integer.toHexString(array[remap])+" at position "+pos+" -> "+(char)AminoAcid.numberToBaseExtended[old]);
		
		return AminoAcid.numberToBaseExtended[old];
	}
	
	private void write(int pos, byte letter){
		int remap=pos/2;
		byte bit1=(byte) (pos&1);
		int old=array[remap];
		int number=AminoAcid.baseToNumberExtended[letter];
		assert(number<16);
		
		if(bit1==0){
			
//			System.out.println("01: "+padZeroes(Integer.toBinaryString(old), 8));
			old=(old&0x0F);
//			System.out.println("02: "+padZeroes(Integer.toBinaryString(old), 8));
			old=(old|(number<<4));
//			System.out.println("03: "+padZeroes(Integer.toBinaryString(old), 8));
		}else{
//			System.out.println("11: "+padZeroes(Integer.toBinaryString(old), 8));
			old=(old&0xF0);
//			System.out.println("12: "+padZeroes(Integer.toBinaryString(old), 8));
			old=(old|number);
//			System.out.println("13: "+padZeroes(Integer.toBinaryString(old), 8));
		}
		
//		System.out.println("Writing "+(char)letter+" at position "+pos
//				+":  "+padZeroes(Integer.toBinaryString(array[remap]), 8)+" -> "+padZeroes(Integer.toBinaryString(old), 8));
		
		array[remap]=(byte)old;
	}
	
	
	private static String padZeroes(String s, int len){
		while(s.length()<len){
			s="0"+s;
		}
		
		if(s.length()>len){
			int chop=s.length()-len;
			s=s.substring(chop);
		}
		return s;
	}
	
	public ChromosomeArrayCompressed(){
		this((byte)-1);
	}
	
	public ChromosomeArrayCompressed(int chrom){
		chromosome=chrom;
	}
	
	
	public void set(int loc, int val){
		
		if(loc>=array.length){//Increase size
			int newlen=1+(3*max(array.length, loc))/2;
			assert(newlen>loc);
			resize(newlen);
			assert(array.length==newlen);
		}
		
		array[loc]=(val>Byte.MAX_VALUE ? Byte.MAX_VALUE : (byte)val);
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}
	
	public byte get(int loc){
		return array[loc];
	}
	
	public String get(int a, int b){
		StringBuilder sb=new StringBuilder(b-a+1);
		for(int i=a; i<=b; i++){
			sb.append((char)get(i));
		}
		return sb.toString();
	}
	
	public void resize(int newlen){
		byte[] temp=new byte[newlen];
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
	
	public boolean extended=false;
	
	public final int chromosome;
	public byte[] array;
	public int maxIndex=-1;
	public int minIndex=Integer.MAX_VALUE;
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -8836873912811373713L;
	
}
