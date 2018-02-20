package align2;

import java.util.ArrayList;

import stream.ConcurrentReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.SiteScore;

import fileIO.ReadWrite;

public class VerifyPairedReads {
	
	
	public static void main(String[] args){

		String fname1=args[0];
		String fname2=args[1];
		
		processPaired(fname1, fname2);
	}
	
	
	
	public static void processPaired(String fname1, String fname2){
		RTextInputStream rtis=new RTextInputStream(fname1, fname2, -1);
		load(rtis);
	}
	
	
	
	public static void load(RTextInputStream rtis){
		
		ConcurrentReadInputStream cris=(USE_CRIS ? new ConcurrentReadInputStream(rtis, -1) : null);
		
		if(cris!=null){
			new Thread(cris).start();
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(reads!=null && reads.size()>0){
				for(Read r : reads){
					if(r.paired()){
						Read r2=r.mate;

						if(r.chrom<1 && r.numSites()>0){
							SiteScore ss=r.topSite(); //Should not be necessary
							r.start=ss.start;
							r.stop=ss.stop;
							r.chrom=ss.chrom;
							r.setStrand(ss.strand);
						}
						if(r2.chrom<1 && r2.numSites()>0){
							SiteScore ss=r2.topSite(); //Should not be necessary
							r2.start=ss.start;
							r2.stop=ss.stop;
							r2.chrom=ss.chrom;
							r2.setStrand(ss.strand);
						}

						assert(r.paired());
						assert(r2.paired());
						assert(r.numericID==r2.numericID);
						assert(r.chrom==r2.chrom) : "\n\n"+r.toText(false)+"\n\n"+r2.toText(false)+"\n\n";
						assert(Tools.absdif(r.start, r2.start)<100000) : "\n\n"+r.toText(false)+"\n\n"+r2.toText(false)+"\n\n";
					}
				}
				
				cris.returnList(ln, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln, ln.list.isEmpty());
			ReadWrite.closeStream(cris);
		}else{
			ArrayList<Read> reads=rtis.nextList();
			while(reads!=null && reads.size()>0){
				for(Read r : reads){
					if(r.paired()){
						Read r2=r.mate;

						if(r.chrom<1 && r.numSites()>0){
							SiteScore ss=r.topSite(); //Should not be necessary
							r.start=ss.start;
							r.stop=ss.stop;
							r.chrom=ss.chrom;
							r.setStrand(ss.strand);
						}
						if(r2.chrom<1 && r2.numSites()>0){
							SiteScore ss=r2.topSite(); //Should not be necessary
							r2.start=ss.start;
							r2.stop=ss.stop;
							r2.chrom=ss.chrom;
							r2.setStrand(ss.strand);
						}

						assert(r.paired());
						assert(r2.paired());
						assert(r.numericID==r2.numericID);
						assert(r.chrom==r2.chrom);
						assert(Tools.absdif(r.start, r2.start)<100000);
					}
				}
				reads=rtis.nextList();
			}
			rtis.close();
		}
	}
	
	
	public static boolean USE_CRIS=true; //Similar speed either way.  "true" may be better with many threads.
	
	
}
