package align2;

import java.util.ArrayList;

import stream.ConcurrentReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.SiteScore;

import dna.Gene;
import fileIO.ReadWrite;

public class MakeErrorQualityHistogram {
	
	public static void main(String[] args){
		
		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
		
		long maxReads=0;
		RTextInputStream rtis=new RTextInputStream(fname1, fname2, maxReads);
		ConcurrentReadInputStream cris=new ConcurrentReadInputStream(rtis, maxReads);
		
		int[][] errors=process(cris);
		printHistogram(errors);
//		System.out.println("*** main() finished ***");
	}
	
	public static void printHistogram(int[][] errors){
		System.out.println("#Error Quality Histogram");
		System.out.println("Quality\tErrors\tMatches\tPercent Errors");
		for(int i=0; i<errors[0].length; i++){
			int e=errors[0][i];
			int m=errors[1][i];
			float percent=e*100f/(e+m);
			System.out.println(i+"\t"+e+"\t"+m+"\t"+String.format("%.3f", percent));
		}
	}
	
	public static int[][] process(ConcurrentReadInputStream cris){
		
		new Thread(cris).start();
		
		int[][] counts=new int[2][50];
		
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

	private static void processList(ArrayList<Read> list, int[][] counts) {
		for(Read r : list){
			processRead(r, counts);
			if(r.mate!=null){
				processRead(r.mate, counts);
			}
		}
	}

	private static void processRead(Read r, int[][] counts) {
		if(!r.mapped() || r.invalid()){return;}
		if(r.mate!=null){
			if(!r.mate.mapped() || r.mate.invalid() || r.mate.match==null || r.mate.containsIndels()){return;}
			int len=Tools.max(r.stop, r.mate.stop)-Tools.min(r.start, r.mate.start);
			if(len<r.bases.length+r.mate.bases.length){return;}
		}
		if(r.match==null){return;}
		

		if(!r.paired()){return;}
		
		if(r.containsIndels()){return;}
		
		if(r.countMismatches()>(2+r.bases.length/10)){return;}
		if(r.avgQuality()<8){return;}
		
		if(r.chrom<1 && r.numSites()>0){
			assert(false) : r.toText(false);
			SiteScore ss=r.topSite(); //Should not be necessary
			r.start=ss.start;
			r.stop=ss.stop;
			r.chrom=ss.chrom;
			r.setStrand(ss.strand);
		}
		
		if(r.strand()==Gene.MINUS){Tools.reverseInPlace(r.match);}
		
		for(int i=0; i<r.match.length; i++){
			byte m=r.match[i];
			byte q=r.quality[i];
			byte b=r.bases[i];
			if(q<0){q=0;}
			if(m=='m'){
				counts[1][q]++;
			}else if(m=='N'){
				assert(q==0 || b!='N') : "\nq="+q+"\n"+r.toText(false); //Could be a noref, though
				//TODO: I need a special symbol for no-ref 
				if(b=='N'){counts[0][q]++;}
			}else{
				assert(m=='S') : r.toText(false);
				counts[0][q]++;
			}
		}
		if(r.strand()==Gene.MINUS){Tools.reverseInPlace(r.match);}
		
	}
	
}
