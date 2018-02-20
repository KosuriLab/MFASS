package align2;

import java.util.ArrayList;

import stream.ConcurrentReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.SiteScore;

import dna.Gene;
import fileIO.ReadWrite;

public class MakeInsertHistogram {
	
	public static void main(String[] args){
		
		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
		
		long maxReads=0;
		RTextInputStream rtis=new RTextInputStream(fname1, fname2, maxReads);
		ConcurrentReadInputStream cris=new ConcurrentReadInputStream(rtis, maxReads);
		
		int[] counts=process(cris);
		printMappedHistogram(counts);
//		System.out.println("*** main() finished ***");
	}
	
	public static void printMappedHistogram(int[] counts){
		System.out.println("#Insert Length Histogram");
		System.out.println("#Reads:    \t"+totalReads);
		System.out.println("#Used:     \t"+used+String.format("\t%.3f", (used*100d/totalReads)));
		
		long wtSum=0;
		for(int i=0; i<counts.length; i++){
			wtSum+=(counts[i]*(long)i);
		}
		System.out.println("#Average:  \t"+String.format("%.3f", (wtSum*1d/used)));
		
		
		System.out.println("Length\tCount\tPercent");
		
		long sum=Tools.sum(counts);
		
		for(int i=0; i<counts.length; i++){
			int c=counts[i];
			float f=c*100f/sum;
			System.out.println(i+"\t"+c+"\t"+String.format("%.3f", f));
		}
	}
	
	public static int[] process(ConcurrentReadInputStream cris){
		
		new Thread(cris).start();

		int[] counts=new int[1000];
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> readlist=ln.list;
		while(!readlist.isEmpty()){
			
			processList(readlist, counts);
			
			cris.returnList(ln, readlist.isEmpty());
			//System.err.println("Waiting on a list...");
			ln=cris.nextList();
			readlist=ln.list;
		}
		
		//System.err.println("Returning a list... (final)");
		assert(readlist.isEmpty());
		cris.returnList(ln, readlist.isEmpty());
		ReadWrite.closeStream(cris);
		
		return counts;
	}

	private static void processList(ArrayList<Read> list, int[] counts) {
		for(Read r : list){
			processRead(r, counts);
//			if(r.mate!=null){
//				processRead(r.mate, mapped, paired);
//			}
		}
	}

	private static void processRead(Read r, int[]counts) {
		totalReads++;
		if(!r.paired()){return;}
		Read r2=r.mate;
		if(r.match==null || r2.match==null || r.invalid() || r2.invalid()){return;}

		if(r.containsIndels()){return;}
		if(r2.containsIndels()){return;}
//		
		if(r.countMismatches()>5){return;}
		if(r.avgQuality()<12){return;}
		if(r2.countMismatches()>5){return;}
		if(r2.avgQuality()<12){return;}
		
		if(r.chrom<1 && r.numSites()>0){
			assert(false) : r.toText(false);
			SiteScore ss=r.topSite(); //Should not be necessary
			r.start=ss.start;
			r.stop=ss.stop;
			r.chrom=ss.chrom;
			r.setStrand(ss.strand);
		}
		
		if(r2.chrom<1 && r2.numSites()>0){
			assert(false) : r2.toText(false);
			SiteScore ss=r2.topSite(); //Should not be necessary
			r2.start=ss.start;
			r2.stop=ss.stop;
			r2.chrom=ss.chrom;
			r2.setStrand(ss.strand);
		}
		
		if(r.chrom!=r2.chrom || Tools.absdif(r.start, r2.start)>2000){
			return;
		}
		
		if(r.chrom<1 || r2.chrom<1){return;}
		
//		int insert;
//		if(r.start<=r2.start){
//			insert=r2.start-r.stop;
//		}else{
//			insert=r.start-r2.stop;
//		}
		
//		int insert=Tools.max(r.stop, r2.stop)-Tools.min(r.start, r2.start);
		int insert;
		if(r.strand()==Gene.PLUS){
			insert=r2.stop-r.start;
		}else{
			insert=r.stop-r2.start;
		}
		
		if(insert<0){insert=0;}
		if(insert>=counts.length){insert=counts.length-1;}
		counts[insert]++;
		used++;
	}

	public static long totalReads=0;
	public static long used=0;
	
}
