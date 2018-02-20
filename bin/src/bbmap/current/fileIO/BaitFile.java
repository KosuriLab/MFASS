package fileIO;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import align2.Tools;

import dna.Data;
import dna.Gene;
import driver.Translator2;

public class BaitFile extends TextFile {
	
	
	public static void main(String[] args){

		String root=ReadWrite.parseRoot(args[0]);
		BaitFile bf=new BaitFile(args[0], true);
		

		int buildIn=37;
		int buildOut=37;
		
		buildIn=Integer.parseInt(args[1]);
		buildOut=Integer.parseInt(args[2]);
		
		boolean makeInt3d=true;
		boolean makeBed=!makeInt3d;

		if(makeInt3d){
			String baitName=args[3];

			int[][][] lines=bf.toLines(buildIn, buildOut, true);

			bf.close();
			ReadWrite.write(lines, root+"baits_"+baitName+"_build"+buildOut+".int3d", false);
		}else{
			translateToBed(bf, buildIn, buildOut);
		}
	}
	
	public static void translateToBed(BaitFile bf, int buildIn, int buildOut){
//		int[][][] lines=bf.toLines(buildIn, buildOut, true);
		int[][][] lines=Data.getBaits();
		bf.close();
		
		int x=0;
		for(int chrom=0; chrom<lines.length; chrom++){
			int[][] set=lines[chrom];
			if(set!=null){
				for(int i=0; i<set[0].length; i++){
					x++;
					System.out.println("chr"+Gene.chromCodes[chrom]+"\t"+set[0][i]+"\t"+set[1][i]+"\tzone"+x+"\t1000\t+");
				}
			}
		}
	}
	

	public BaitFile(String fname, boolean tryAllExtensions) {
		super(fname, false, tryAllExtensions);
	}
	
	/** Format: [chrom][start array, stop array, numBaits array] */
	public int[][][] toLines(int buildIn, int buildOut, boolean condense){
		
		String s=null;
		HashMap<BaitLine, BaitLine> table=new HashMap<BaitLine, BaitLine>(8000);
		
		for(s=nextLine(); s!=null; s=nextLine()){
			BaitLine b=new BaitLine(s);
			
			
			if(buildIn!=buildOut && b.chrom<26){
				b=b.translate(buildIn, buildOut);
			}

			if(b!=null && b.chrom<26){

				BaitLine old=table.get(b);
				if(old==null){
					table.put(b, b);
				}else{
					old.add(b);
				}
			}
		}
		
		ArrayList<BaitLine> master=new ArrayList<BaitLine>(table.size());
		master.addAll(table.values());
		table=null;
		
		if(condense){condense(master);}
		
		ArrayList<BaitLine>[] lists=new ArrayList[26];
		for(BaitLine b : master){
			if(lists[b.chrom]==null){lists[b.chrom]=new ArrayList<BaitLine>();}
			lists[b.chrom].add(b);
		}

		int[][][] out=new int[26][3][];
		for(int i=0; i<lists.length; i++){
			if(lists[i]!=null){
				Collections.sort(lists[i]);
				out[i][0]=new int[lists[i].size()];
				out[i][1]=new int[lists[i].size()];
				out[i][2]=new int[lists[i].size()];
				
				for(int j=0; j<lists[i].size(); j++){
					BaitLine b=lists[i].get(j);
					out[i][0][j]=b.start;
					out[i][1][j]=b.stop;
					out[i][2][j]=b.names.size();
				}
				
			}else{
				out[i][0]=out[i][1]=out[i][2]=new int[0];
			}
		}
		
		return out;
	}
	
	public static void condense(ArrayList<BaitLine> master){
		Collections.sort(master);
		
		int merged=0;
		BaitLine prev=null;
		for(int i=0; i<master.size(); i++){
			BaitLine b=master.get(i);
			if(prev==null){prev=b;}
			else{
				if(prev.touches(b)){
					assert(b.touches(prev));
					assert(prev.chrom==b.chrom && touch(prev.start, prev.stop, b.start, b.stop));
					prev.merge(b);
					master.set(i, null);
					merged++;
				}else{
					assert(!b.touches(prev));
					assert(prev.chrom!=b.chrom || !touch(prev.start, prev.stop, b.start, b.stop));
					prev=b;
				}
			}
		}
		System.err.println("Merged "+merged);
		Tools.condense(master);
		
//		merged=0;
//		Collections.sort(master);
//		prev=null;
//		for(int i=0; i<master.size(); i++){
//			BaitLine b=master.get(i);
//			if(prev==null){prev=b;}
//			else{
//				assert(prev!=b);
//				if(prev.touches(b)){
//					assert(b.touches(prev));
//					assert(prev.chrom==b.chrom && touch(prev.start, prev.stop, b.start, b.stop));
//					prev.merge(b);
//					master.set(i, null);
//					merged++;
//				}else{
//					assert(!b.touches(prev));
//					assert(prev.chrom!=b.chrom || !touch(prev.start, prev.stop, b.start, b.stop));
//					prev=b;
//				}
//			}
//		}
//		System.err.println("Merged "+merged);
//		Tools.condense(master);
		
		Collections.sort(master);
		for(int i=1; i<master.size(); i++){
			assert(!master.get(i-1).touches(master.get(i))) : "\n"+i+", "+master.get(i-1)+", "+master.get(i)+"\n";
		}
	}
	
	public String nextLine(){
		String line=readLine();
		while(line!=null && (!line.startsWith("chr") || line.contains("random") || line.startsWith("#"))){
			line=readLine();
		}
		return line;
	}

	
	public static boolean touch(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=(b1+1) && b2>=(a1-1);
	}
	
	private static class BaitLine implements Comparable<BaitLine> {
		
		public BaitLine(String s){
			try{
			String[] split=s.split("\t", -1);
			chrom=Gene.toChromosome(split[0]);
//			assert(chrom<26) : s;
			start=Integer.parseInt(split[1])-1;
			stop=Integer.parseInt(split[2])-1;
			if(split.length>3){names.add(split[3]);}
			}catch(Exception e){
				System.err.println(s);
				throw new RuntimeException(e);
			}
		}
		
		public BaitLine(int chr, int sta, int sto, ArrayList<String> nam){
			chrom=chr;
			start=sta;
			stop=sto;
			names.addAll(nam);
		}
		
		private BaitLine translate(int buildIn, int buildOut){

			int[] startTrans=Translator2.translate(buildIn, buildOut, chrom, Gene.PLUS, start);
			int[] stopTrans=Translator2.translate(buildIn, buildOut, chrom, Gene.PLUS, stop);
			
			if(startTrans==null || stopTrans==null){return null;}
			
			if(startTrans[0]!=stopTrans[0]){return null;}//different chromosomes
			
			int chrom2=startTrans[0];
			int start2=startTrans[2];
			int stop2=stopTrans[2];
			if(start2>stop2){
				int temp=start2;
				start2=stop2;
				stop2=temp;
			}
			
			int len=stop-start;
			int len2=stop2-start2;
			int dif=(len>len2 ? len-len2 : len2-len);
			
//			assert(len>0 && len<3000) : this; //Baits should be 120 long, IIRC. (**No longer true!**)
			if(dif>50){return null;}
			
			return new BaitLine(chrom2, start2, stop2, names);
		}
		
		public int[] toInt(){
			return new int[] {chrom, start, stop, names.size()};
		}
		
		public String toString(){
			return chrom+"\t"+start+"\t"+stop;
		}
		
		final int chrom;
		int start;
		int stop;
		final ArrayList<String> names=new ArrayList<String>(2);
		
		public void add(BaitLine other){
			assert(this.equals(other));
			assert(this!=other);
			names.addAll(other.names);
		}

		@Override
		public int compareTo(BaitLine other) {
			int r;
			
			r=chrom-other.chrom;
			if(r!=0){return r;}
			
			r=start-other.start;
			if(r!=0){return r;}
			
			r=stop-other.stop;
			if(r!=0){return r;}
			
			return 0;
		}
		
		public void merge(BaitLine b){
			assert(touches(b));
			assert(chrom==b.chrom);
			assert(start<=b.start) : this+", "+b;
			start=Data.min(start, b.start);
			stop=Data.max(stop, b.stop);
		}
		
		public boolean touches(BaitLine b){
			return (chrom==b.chrom && b.stop>=start-1 && b.start<=stop+1);
		}
		
		@Override
		public boolean equals(Object other){
			return equals((BaitLine)other);
		}
		
		public boolean equals(BaitLine other){
			return compareTo(other)==0;
		}
		
		@Override
		public int hashCode(){
			return start^chrom;
		}
		
	}
	

}
