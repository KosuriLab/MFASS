package align2;

import java.util.ArrayList;

import stream.Read;

import fileIO.TextStreamWriter;


/**
 * @author Brian Bushnell
 * @date Mar 18, 2013
 *
 */
public class ReadStats {
	
	public ReadStats(){this(true);}
		
	public ReadStats(boolean addToList){
		if(addToList){
			synchronized(ReadStats.class){
				objectList.add(this);
			}
		}

		if(COLLECT_QUALITY_STATS){
			qualLength=new long[2][MAXLEN];
			qualSum=new long[2][MAXLEN];
		}else{
			qualLength=null;
			qualSum=null;
		}

		if(COLLECT_MATCH_STATS){
			matchSum=new long[2][MAXLEN];
			delSum=new long[2][MAXLEN];
			insSum=new long[2][MAXLEN];
			subSum=new long[2][MAXLEN];
			nSum=new long[2][MAXLEN];
			clipSum=new long[2][MAXLEN];
			otherSum=new long[2][MAXLEN];
		}else{
			matchSum=null;
			delSum=null;
			insSum=null;
			subSum=null;
			nSum=null;
			clipSum=null;
			otherSum=null;
		}

		if(COLLECT_INSERT_STATS){
			insertHist=new LongList(MAXLEN);
		}else{
			insertHist=null;
		}
	}
	
	public static ReadStats mergeAll(){
		ReadStats x=new ReadStats(false);
		for(ReadStats rs : objectList){
			if(COLLECT_QUALITY_STATS){
				for(int i=0; i<MAXLEN; i++){
					x.qualLength[0][i]+=rs.qualLength[0][i];
					x.qualLength[1][i]+=rs.qualLength[1][i];
					x.qualSum[0][i]+=rs.qualSum[0][i];
					x.qualSum[1][i]+=rs.qualSum[1][i];
				}
			}
			if(COLLECT_MATCH_STATS){
				for(int i=0; i<MAXLEN; i++){
					x.matchSum[0][i]+=rs.matchSum[0][i];
					x.matchSum[1][i]+=rs.matchSum[1][i];
					x.delSum[0][i]+=rs.delSum[0][i];
					x.delSum[1][i]+=rs.delSum[1][i];
					x.insSum[0][i]+=rs.insSum[0][i];
					x.insSum[1][i]+=rs.insSum[1][i];
					x.subSum[0][i]+=rs.subSum[0][i];
					x.subSum[1][i]+=rs.subSum[1][i];
					x.nSum[0][i]+=rs.nSum[0][i];
					x.nSum[1][i]+=rs.nSum[1][i];
					x.clipSum[0][i]+=rs.clipSum[0][i];
					x.clipSum[1][i]+=rs.clipSum[1][i];
					x.otherSum[0][i]+=rs.otherSum[0][i];
					x.otherSum[1][i]+=rs.otherSum[1][i];
				}
			}
			if(COLLECT_INSERT_STATS){
				for(int i=0; i<rs.insertHist.size; i++){
					long y=rs.insertHist.get(i);
					if(y>0){x.insertHist.increment(i, y);}
				}
			}
			
		}
		return x;
	}
	
	public void addToQualityHistogram(final Read r){
		if(r==null){return;}
		addToQualityHistogram(r, r.obj, 0);
		if(r.mate!=null){addToQualityHistogram(r.mate, r.obj, 1);}
	}
	
	private void addToQualityHistogram(final Read r, Object obj, final int pairnum){
		if(r==null || r.quality==null || r.quality.length<1){return;}
		final byte[] qual;
		if(obj!=null && obj.getClass()==TrimRead.class){
			qual=(pairnum==0 ? ((TrimRead)obj).qual1 : ((TrimRead)obj).qual2);
		}else{
			qual=r.quality;
		}
		final int limit=Tools.min(qual.length, MAXLEN);
		final long[] ql=qualLength[pairnum], qs=qualSum[pairnum];
		ql[limit-1]++;
		for(int i=0; i<limit; i++){qs[i]+=qual[i];}
	}
	
	public void addToMatchHistogram(final Read r){
		if(r==null){return;}
		addToMatchHistogram(r, 0);
		if(r.mate!=null){addToMatchHistogram(r.mate, 1);}
	}
	
	private void addToMatchHistogram(final Read r, final int pairnum){
		if(r==null || r.bases==null || r.bases.length<1 || !r.mapped()){return;}
		final byte[] bases=r.bases, match=r.match;
		final int limit=Tools.min(bases.length, MAXLEN);
		final long[] ms=matchSum[pairnum], ds=delSum[pairnum], is=insSum[pairnum],
				ss=subSum[pairnum], ns=nSum[pairnum], cs=clipSum[pairnum], os=otherSum[pairnum];
		
		if(match==null){
			for(int i=0; i<limit; i++){
				byte b=bases[i];
				if(b=='N'){ns[i]++;}
				else{os[i]++;}
			}
		}else{
			final boolean plus=(r.strand()==0);
			int rpos=0;
			byte lastm='A';
			for(int mpos=0; mpos<match.length && rpos<limit; mpos++){
				byte b=bases[plus ? rpos : bases.length-rpos-1];
				byte m=match[mpos];
				if(b=='N'){
					if(m=='D'){
						if(lastm!=m){ds[rpos]++;}
						rpos--;
					}else{ns[rpos]++;}
				}else{
					if(m=='m'){
						ms[rpos]++;
					}else if(m=='S'){
						ss[rpos]++;
					}else if(m=='I'){
						is[rpos]++;
					}else if(m=='N'){
//						assert(false) : "\n"+r+"\n"+new String(Data.getChromosome(r.chrom).getBytes(r.start, r.stop))+"\nrpos="+rpos+", mpos="+mpos;
						os[rpos]++;
					}else if(m=='C'){
//						assert(false) : r;
						cs[rpos]++;
					}else if(m=='D'){
						if(lastm!=m){ds[rpos]++;}
						rpos--;
					}else{
						os[rpos]++;
						assert(false) : "For read "+r.numericID+", unknown symbol in match string: ASCII "+m+" = "+(char)m;
					}
				}
				rpos++;
				lastm=m;
			}
		}
	}
	
	public void addToInsertHistogram(final Read r, boolean ignoreMappingStrand){
		if(verbose){
			System.err.print(r.numericID);
			if(r==null || r.mate==null || !r.mapped() || !r.mate.mapped() || !r.paired()){
				System.err.println("\n");
			}else{
				System.err.println("\t"+r.strand()+"\t"+r.insertSizeMapped(ignoreMappingStrand)+"\t"+r.mate.insertSizeMapped(ignoreMappingStrand));
			}
		}
		if(r==null || r.mate==null || !r.mapped() || !r.mate.mapped() || !r.paired()){return;}
		int x=Tools.min(MAXINSERTLEN, r.insertSizeMapped(ignoreMappingStrand));
		if(x>0){insertHist.increment(x, 1);}
//		assert(x!=1) : "\n"+r+"\n\n"+r.mate+"\n";
//		System.out.println("Incrementing "+x);
	}
	
	public void writeQualityToFile(String fname, boolean writePaired){
		TextStreamWriter tsw=new TextStreamWriter(fname, OVERWRITE, false, false);
		tsw.start();
		tsw.print("#BaseNum\tRead1"+(writePaired ? "\tRead2" : "")+"\n");
		
		final long[] qs1=qualSum[0], qs2=qualSum[1], ql1=qualLength[0], ql2=qualLength[1];
		
		for(int i=MAXLEN-2; i>=0; i--){
			ql1[i]+=ql1[i+1];
			ql2[i]+=ql2[i+1];
		}
		
		if(writePaired){
			for(int i=0; i<MAXLEN && (ql1[i]>0 || ql2[i]>0); i++){
				int a=i+1;
				double b=qs1[i]/(double)Tools.max(1, ql1[i]);
				double c=qs2[i]/(double)Tools.max(1, ql2[i]);
				tsw.print(String.format("%d\t%.3f\t%.3f\n", a, b, c));
			}
		}else{
			for(int i=0; i<MAXLEN && ql1[i]>0; i++){
				int a=i+1;
				double b=qs1[i]/(double)Tools.max(1, ql1[i]);
				tsw.print(String.format("%d\t%.3f\n", a, b));
			}
		}
		tsw.poison();
		tsw.waitForFinish();
	}
	
	public void writeMatchToFile(String fname, boolean writePaired){
		if(!writePaired){
			writeMatchToFileUnpaired(fname);
			return;
		}
		TextStreamWriter tsw=new TextStreamWriter(fname, OVERWRITE, false, false);
		tsw.start();
		tsw.print("#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1\tMatch2\tSub2\tDel2\tIns2\tN2\tOther2\n");
		
		final long[] ms1=matchSum[0], ds1=delSum[0], is1=insSum[0],
				ss1=subSum[0], ns1=nSum[0], cs1=clipSum[0], os1=otherSum[0];
		final long[] ms2=matchSum[1], ds2=delSum[1], is2=insSum[1],
				ss2=subSum[1], ns2=nSum[1], cs2=clipSum[1], os2=otherSum[1];
		
		for(int i=0; i<MAXLEN; i++){
			int a=i+1;
			long sum1=ms1[i]+is1[i]+ss1[i]+ns1[i]+cs1[i]+os1[i]; //no deletions
			long sum2=ms2[i]+is2[i]+ss2[i]+ns2[i]+cs2[i]+os2[i]; //no deletions
			if(sum1==0 && sum2==0){break;}
			double inv1=1.0/(double)Tools.max(1, sum1);
			double inv2=1.0/(double)Tools.max(1, sum2);

			tsw.print(String.format("%d", a));
			tsw.print(String.format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f", 
					ms1[i]*inv1, ss1[i]*inv1, ds1[i]*inv1, is1[i]*inv1, ns1[i]*inv1, (os1[i]+cs1[i])*inv1));
			tsw.print(String.format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f", 
					ms2[i]*inv2, ss2[i]*inv2, ds2[i]*inv2, is2[i]*inv2, ns2[i]*inv2, (os2[i]+cs2[i])*inv2)
//					+", "+ms2[i]+", "+is2[i]+", "+ss2[i]+", "+ns2[i]+", "+cs2[i]+", "+os2[i]
					);
			tsw.print("\n");
		}
		tsw.poison();
		tsw.waitForFinish();
	}
	
	public void writeMatchToFileUnpaired(String fname){
		TextStreamWriter tsw=new TextStreamWriter(fname, OVERWRITE, false, false);
		tsw.start();
		tsw.print("#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1\n");
		
		final long[] ms1=matchSum[0], ds1=delSum[0], is1=insSum[0],
				ss1=subSum[0], ns1=nSum[0], cs1=clipSum[0], os1=otherSum[0];
		
		for(int i=0; i<MAXLEN; i++){
			int a=i+1;
			long sum1=ms1[i]+is1[i]+ss1[i]+ns1[i]+cs1[i]+os1[i]; //no deletions
			if(sum1==0){break;}
			double inv1=1.0/(double)Tools.max(1, sum1);

			tsw.print(String.format("%d", a));
			tsw.print(String.format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f", 
					ms1[i]*inv1, ss1[i]*inv1, ds1[i]*inv1, is1[i]*inv1, ns1[i]*inv1, (os1[i]+cs1[i])*inv1)
//					+", "+ms1[i]+", "+is1[i]+", "+ss1[i]+", "+ns1[i]+", "+cs1[i]+", "+os1[i]
					);
			tsw.print("\n");
		}
		tsw.poison();
		tsw.waitForFinish();
	}
	
	public void writeInsertToFile(String fname){
		TextStreamWriter tsw=new TextStreamWriter(fname, OVERWRITE, false, false);
		tsw.start();
		tsw.print("#InsertSize\tCount\n");
		
		for(int i=0; i<insertHist.size; i++){
			long x=insertHist.get(i);
			if(x>0 || !skipZeroInsertCount){
				tsw.print(i+"\t"+x+"\t"+"\n");
			}
		}
		tsw.poison();
		tsw.waitForFinish();
	}
	
	public final long[][] qualLength;
	public final long[][] qualSum;
	
	public final long[][] matchSum;
	public final long[][] delSum;
	public final long[][] insSum;
	public final long[][] subSum;
	public final long[][] nSum;
	public final long[][] clipSum;
	public final long[][] otherSum;
	
	public final LongList insertHist;

	public static final int MAXLEN=2000;
	public static final int MAXINSERTLEN=24000;
	
	public static ArrayList<ReadStats> objectList=new ArrayList<ReadStats>();
	public static boolean COLLECT_QUALITY_STATS=false;
	public static boolean COLLECT_MATCH_STATS=false;
	public static boolean COLLECT_INSERT_STATS=false;
	public static String QUAL_HIST_FILE=null;
	public static String MATCH_HIST_FILE=null;
	public static String INSERT_HIST_FILE=null;
	public static boolean OVERWRITE=false;
	public static final boolean verbose=false;
	
	public static boolean skipZeroInsertCount=true;
	
}
