package stream;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import align2.MSA;
import align2.Shared;
import align2.Tools;

import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;



public final class SiteScore implements Comparable<SiteScore>, Cloneable{
	
	public SiteScore(int chrom_, byte strand_, int start_, int stop_, int hits_, int quickScore_){
		start=start_;
		stop=stop_;
		hits=hits_;
		quickScore=quickScore_;
		score=quickScore_;
		chrom=chrom_;
		strand=strand_;
//		assert(chrom_>=0) : this.toText()+"\nchrom_="+chrom_+", strand_="+strand_+", start_="+start_+", stop_="+stop_+", hits_="+hits_+", quickScore_="+quickScore_;
		assert(start_<=stop_) : this.toText()+"\nchrom_="+chrom_+", strand_="+strand_+", start_="+start_+", stop_="+stop_+", hits_="+hits_+", quickScore_="+quickScore_;
	}
	
	public SiteScore(int chrom_, byte strand_, int start_, int stop_, int hits_, int quickScore_, boolean rescued_, boolean perfect_){
		start=start_;
		stop=stop_;
		hits=hits_;
		quickScore=quickScore_;
		score=quickScore_;
		chrom=chrom_;
		strand=strand_;
		rescued=rescued_;
		perfect=perfect_;
		semiperfect=perfect;
		assert(start_<=stop_) : this.toText();
	}
	
	@Override
	public int compareTo(SiteScore other) {
		int x=other.score-score;
		if(x!=0){return x;}
		
		x=other.slowScore-slowScore;
		if(x!=0){return x;}
		
		x=other.pairedScore-pairedScore;
		if(x!=0){return x;}
		
		x=other.quickScore-quickScore;
		if(x!=0){return x;}
		
		x=chrom-other.chrom;
		if(x!=0){return x;}
		
		x=start-other.start;
		return x;
	}
	
	public boolean equals(Object other){
		return compareTo((SiteScore)other)==0;
	}
	
	public String toString(){
		return toText().toString();
	}
	
//	9+2+1+9+9+1+1+4+4+4+4+gaps
	public CharSequence toText(){
		StringBuilder sb=new StringBuilder(53+(gaps==null ? 0 : gaps.length*10));
		sb.append(chrom);
		sb.append(',');
		sb.append(strand);
		sb.append(',');
		sb.append(start);
		sb.append(',');
		sb.append(stop);
		sb.append(',');
		sb.append((rescued ? 1 : 0));
		sb.append(',');
		sb.append((semiperfect ? 1 : 0));
		sb.append((perfect ? 1 : 0));
		sb.append(',');
		sb.append(hits);
		sb.append(',');
		sb.append(quickScore);
		sb.append(',');
		sb.append(slowScore);
		sb.append(',');
		sb.append(pairedScore);
		sb.append(',');
		sb.append(score);
		
		if(gaps!=null){
			sb.append(',');
			for(int i=0; i<gaps.length; i++){
				if(i>0){sb.append('~');}
				sb.append(gaps[i]);
			}
		}
		
		if(match!=null){
			sb.append(',');
			final char[] buffer=Shared.getTLCB(match.length);
			for(int i=0; i<match.length; i++){buffer[i]=(char)match[i];}
			sb.append(buffer);
		}
		
		return sb;
//		chrom+","+strand+","+start+","+stop+","+(rescued ? 1 : 0)+","+
//		(perfect ? 1 : 0)+","+quickScore+","+slowScore+","+pairedScore+","+score;
	}
	
//	9+2+1+9+9+1+1+4+4+4+4+gaps
	public ByteBuilder toBytes(ByteBuilder sb){
		if(sb==null){sb=new ByteBuilder(53+(gaps==null ? 0 : gaps.length*10));}
		sb.append(chrom);
		sb.append(',');
		sb.append((int)strand);
		sb.append(',');
		sb.append(start);
		sb.append(',');
		sb.append(stop);
		sb.append(',');
		sb.append((rescued ? 1 : 0));
		sb.append(',');
		sb.append((semiperfect ? 1 : 0));
		sb.append((perfect ? 1 : 0));
		sb.append(',');
		sb.append(hits);
		sb.append(',');
		sb.append(quickScore);
		sb.append(',');
		sb.append(slowScore);
		sb.append(',');
		sb.append(pairedScore);
		sb.append(',');
		sb.append(score);
		
		if(gaps!=null){
			sb.append(',');
			for(int i=0; i<gaps.length; i++){
				if(i>0){sb.append('~');}
				sb.append(gaps[i]);
			}
		}
		
		if(match!=null){
			sb.append(',');
			sb.append(match);
		}
		
		return sb;
//		chrom+","+strand+","+start+","+stop+","+(rescued ? 1 : 0)+","+
//		(perfect ? 1 : 0)+","+quickScore+","+slowScore+","+pairedScore+","+score;
	}
	
	public boolean isSemiPerfect(byte[] bases){
		if(bases.length!=stop-start+1){return false;}
		byte[] ref=Data.getChromosome(chrom).array;

		//This block handles cases where the read runs outside the reference
		//Of course, padding the reference with 'N' would be better, but...
		int readStart=0;
		int readStop=bases.length;
		final int refStop=start+bases.length;
		int maxNoref=bases.length/2;

		if(start<0){
			readStart=0-start;
		}
		if(refStop>ref.length){
			int dif=(refStop-ref.length);
			readStop-=dif;
		}

		for(int i=readStart; i<readStop; i++){
			byte c=bases[i];
			byte r=ref[start+i];
			
//			assert(Character.isUpperCase(c) && Character.isUpperCase(r));
			if(c=='N'){return false;}
			if(c!=r){
				maxNoref--;
				if(maxNoref<0 || r!='N'){return false;}
			}
		}
		return true;
	}
	
	public boolean isPerfect(byte[] bases){
		if(bases.length!=stop-start+1 || start<0){return false;}
		byte[] ref=Data.getChromosome(chrom).array;
		if(stop>=ref.length){return false;}
		
		for(int i=0; i<bases.length; i++){
			byte c=bases[i];
			byte r=ref[start+i];
			assert(Character.isUpperCase(c) && Character.isUpperCase(r));

			if((c!=r /* && (Character.toUpperCase(c)!=Character.toUpperCase(r))*/) || c=='N'){
				return false;
			}
		}
		return true;
	}
	
	
	public boolean setPerfectFlag(int maxScore, byte[] bases){
		if(maxScore==slowScore){
			assert(isPerfect(bases));
			return perfect=semiperfect=true;
		}
		return setPerfect(bases, false);
	}
	
	/** Sets "perfect" and "semiperfect" flags */
	public boolean setPerfect(byte[] bases){return setPerfect(bases, false);}
	
	/** Sets "perfect" and "semiperfect" flags, optionally assuming "perfect" flag is correct. */
	public boolean setPerfect(byte[] bases, boolean assumePerfectCorrect){
		if(bases.length!=stop-start+1){
			assert(!perfect || !assumePerfectCorrect) : perfect+", "+toString()+", "+
					new String(Data.getChromosome(chrom).array, Tools.max(0, start), (Tools.min(Data.chromLengths[chrom], stop)-start));
			perfect=false;
			semiperfect=false;
			assert(Read.CHECKSITE(this, bases, 0)); //123
			return perfect;
		}
		byte[] ref=Data.getChromosome(chrom).array;
		
		perfect=semiperfect=true;
		int refloc=start, readloc=0, N=0, max=Tools.min(stop, ref.length-1), nlimit=bases.length/2;
		if(start<0){
			N-=start;
			readloc-=start;
			refloc-=start;
			assert(!perfect || !assumePerfectCorrect);
			perfect=false;
		}
		if(stop>=ref.length){
			N+=(stop-ref.length+1);
			assert(!perfect || !assumePerfectCorrect);
			perfect=false;
		}
		if(N>nlimit){
			perfect=semiperfect=false;
			assert(Read.CHECKSITE(this, bases, 0)); //123
			return perfect;
		}
		
		final byte bn=(byte)'N';
		for(; refloc<=max; refloc++, readloc++){
			final byte c=bases[readloc];
			final byte r=ref[refloc];
			assert(Character.isUpperCase(r) && Character.isUpperCase(c)) :
				"\nAn input read appears to contain a non-upper-case base.  Please rerun with the 'touppercase' flag.\n"+
				r+", "+c+"\n";
			if(c!=r || c==bn){
				perfect=false;
				if(c==bn){semiperfect=false;}
				if(r!=bn || (N=N+1)>nlimit){
					semiperfect=false;
					assert(Read.CHECKSITE(this, bases, 0)); //123
					return semiperfect;
				}
			}
		}
		
		semiperfect=(semiperfect && (N<=nlimit));
		perfect=(perfect && semiperfect && (N==0));
		assert(Read.CHECKSITE(this, bases, 0)); //123
		return perfect;
	}
	
	public final boolean overlaps(SiteScore ss){
		return chrom==ss.chrom && strand==ss.strand && overlap(start, stop, ss.start, ss.stop);
	}
	public final boolean overlaps(SiteScore ss, boolean ignoreStrand){
		return chrom==ss.chrom && (ignoreStrand || strand==ss.strand) && overlap(start, stop, ss.start, ss.stop);
	}
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1; 
	}
	
	public static String header() {
		return "chrom,strand,start,stop,rescued,semiperfect+perfect,hits,quickScore,slowScore,pairedScore,score,match";
	}
	
	public static SiteScore fromText(String s){
//		System.err.println("Trying to make a SS from "+s);
		String line[]=s.split(",");
		
		SiteScore ss;

		assert(line.length>=11 && line.length<=13) : "\n"+line.length+"\n"+s+"\n"+Arrays.toString(line)+"\n";
		int chrom=Byte.parseByte(line[0].charAt(0)=='*' ? line[0].substring(1) : line[0]);
		byte strand=Byte.parseByte(line[1]);
		int start=Integer.parseInt(line[2]);
		int stop=Integer.parseInt(line[3]);
		boolean rescued=Integer.parseInt(line[4])==1;
//		[1, 1, 9397398, 9398220, 0, 00, 20, 8701, 9084, 0, 9084, 9397398~9397471~9398145~9398220]
		int p=Integer.parseInt(line[5], 2);
//		assert(false) : line[5]+"->"+p;
		boolean perfect=(p&1)==1;
		boolean semiperfect=(p&2)==2;
		int hits=Integer.parseInt(line[6]);
		int quickScore=Integer.parseInt(line[7]);
		int swscore=Integer.parseInt(line[8]);
		int pairedScore=Integer.parseInt(line[9]);
		int score=Integer.parseInt(line[10]);
		ss=new SiteScore(chrom, strand, start, stop, hits, quickScore, rescued, perfect);
		ss.score=score;
		ss.slowScore=swscore;
		ss.pairedScore=pairedScore;
		ss.semiperfect=semiperfect;
		
		if(line.length>11){
			String[] gstring=line[11].split("~");
			ss.gaps=new int[gstring.length];
			for(int i=0; i<gstring.length; i++){
				ss.gaps[i]=Integer.parseInt(gstring[i]);
			}
		}
		
		if(line.length>12){
			ss.match=line[12].getBytes();
		}
		
		return ss;
	}
	
	public boolean positionalMatch(SiteScore b, boolean testGaps){
//		return chrom==b.chrom && strand==b.strand && start==b.start && stop==b.stop;
		if(chrom!=b.chrom || strand!=b.strand || start!=b.start || stop!=b.stop){
			return false;
		}
		if(!testGaps || (gaps==null && b.gaps==null)){return true;}
		if((gaps==null) != (b.gaps==null)){return false;}
		if(gaps.length!=b.gaps.length){return false;}
		for(int i=0; i<gaps.length; i++){
			if(gaps[i]!=b.gaps[i]){return false;}
		}
		return true;
	}
	
	public byte[] getScaffoldName(boolean requireSingleScaffold){
		byte[] name=null;
		if(!requireSingleScaffold || Data.isSingleScaffold(chrom, start, stop)){
			int idx=Data.scaffoldIndex(chrom, (start+stop)/2);
			name=Data.scaffoldNames[chrom][idx];
			//				int scaflen=Data.scaffoldLengths[chrom][idx];
			//				a1=Data.scaffoldRelativeLoc(chrom, start, idx);
			//				b1=a1-start1+stop1;
		}
		return name;
	}
	
	public static class PositionComparator implements Comparator<SiteScore>{
		
		private PositionComparator(){}
		
		@Override
		public int compare(SiteScore a, SiteScore b) {
			if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
			if(a.start!=b.start){return a.start-b.start;}
			if(a.stop!=b.stop){return a.stop-b.stop;}
			if(a.strand!=b.strand){return a.strand-b.strand;}
			if(a.score!=b.score){return b.score-a.score;}
			if(a.slowScore!=b.slowScore){return b.slowScore-a.slowScore;}
			if(a.quickScore!=b.quickScore){return b.quickScore-a.quickScore;}
			if(a.perfect!=b.perfect){return a.perfect ? -1 : 1;}
			if(a.rescued!=b.rescued){return a.rescued ? 1 : -1;}
			return 0;
		}
		
		public void sort(List<SiteScore> list){
			if(list==null || list.size()<2){return;}
			Collections.sort(list, this);
		}
		
		public void sort(SiteScore[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}
	
	public SiteScore copy(){
		SiteScore ss2=this.clone();
		if(gaps!=null){ss2.gaps=ss2.gaps.clone();}
		return ss2;
	}
	
	public SiteScore clone(){
		try {
			return (SiteScore)super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException();
	}
	
	public boolean isInBounds(){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return (start>=0 && stop<=cha.maxIndex);
	}
	
	protected boolean matchContainsXY(){
		if(match==null || match.length<1){return false;}
		final byte a=match[0], b=match[match.length-1];
		return (a=='X' ||a=='Y' || b=='X' || b=='Y');
	}
	
	public boolean isCorrect(int chrom_, byte strand_, int start_, int stop_, int thresh){
		if(chrom_!=chrom || strand_!=strand){return false;}
		if(thresh<=0){return start_==start && stop_==stop;}
		return Tools.absdif(start_, start)<=thresh || Tools.absdif(stop_, stop)<=thresh;
	}
	
	/** TODO: Test
	 * Attempt to extend match/N symbols where there are X and Y symbols
	 * */
	public boolean fixXY(byte[] bases, boolean nullifyOnFailure, MSA msa){
		if(verbose){System.err.println("ss.fixXY()");}
		
		if(!matchContainsXY()){return true;}
		
		boolean disable=false;
		if(disable){
			if(nullifyOnFailure){
				match=null;
			}
//			else if(clipOnFailure){
//				for(int i=0; i<match.length; i++){
//					if(match[i]=='X' || match[i]=='Y'){match[i]='C';}
//				}
//			}
			return false;
		}
		
		if(match==null || match.length<1){return false;}
		final ChromosomeArray ca=Data.getChromosome(chrom);
		final int tip=3;
		boolean success=true;
		
		{//Process left side
			int mloc=0;
			while(mloc<match.length && (match[mloc]=='X' || match[mloc]=='Y')){mloc++;}
			if(mloc>=match.length || mloc>=bases.length){success=false;}
			else if(mloc>0){
				mloc--;
				int rloc=start+mloc, cloc=mloc;
				while(mloc>=0){
					byte m=match[mloc];
					byte c=bases[cloc];
					byte r=ca.get(rloc);
					assert(m=='X' || m=='Y') : (char)m+", "+mloc+", "+(char)c+", "+(char)r+"\n"+new String(bases)+"\n"+this.toString();
					if(r=='N' || c=='N'){match[mloc]='N';}
					else if(c==r){match[mloc]='m';}
					else if(mloc<=tip){match[mloc]='S';}
					else{
						success=false;
						break;
					}
					mloc--;
					rloc--;
					cloc--;
				}
			}
		}
		
		if(success){//Process right side
			int mloc=match.length-1;
			while(mloc>=0 && (match[mloc]=='X' || match[mloc]=='Y')){mloc--;}
			int dif=match.length-1-mloc;
			if(mloc<0){success=false;}
			else if(dif>0){
				mloc++;
				int rloc=stop-dif+1, cloc=bases.length-dif;
				if(cloc<0){success=false;}
				else{
					final int tip2=match.length-tip;
					while(mloc<match.length){
						byte m=match[mloc];
						byte c=bases[cloc];
						byte r=ca.get(rloc);
						assert(m=='X' || m=='Y') : (char)m+", "+mloc+", "+(char)c+", "+(char)r+"\n"+new String(bases)+"\n"+this.toString();
						if(r=='N' || c=='N'){match[mloc]='N';}
						else if(c==r){match[mloc]='m';}
						else if(mloc>=tip2){match[mloc]='S';}
						else{
							success=false;
							break;
						}
						mloc++;
						rloc++;
						cloc++;
					}
				}
			}
		}
		
		success=success && !matchContainsXY();
		if(!success && nullifyOnFailure){match=null;}
		
//		assert(false) : "TODO: Alter score to reflect changes"; //TODO
		if(match!=null){slowScore=msa.score(match);}
		
		return success;
	}

//	public boolean plus(){return strand()==Gene.PLUS;}
//	public boolean minus(){return strand()==Gene.MINUS;}
//	
//	public final byte strand(){return (byte)(flags&strandMask);}
//	public boolean rescued(){return (flags&rescuedMask)!=0;}
//	public boolean perfect(){return (flags&perfectMask)!=0;}
//	public boolean semiperfect(){return (flags&semiperfectMask)!=0;}
//	
//	public final int setStrand(int x){
//		assert(x==0 || x==1);
//		if(x==0){flags=(flags&~strandMask);}
//		else{flags=(flags|strandMask);}
//		assert(strand()==x);
//		return x;
//	}
//	public boolean setRescued(boolean b){
//		if(b){flags=(flags|rescuedMask);}
//		else{flags=(flags&~rescuedMask);}
//		assert(rescued()==b);
//		return b;
//	}
//	public boolean setPerfect(boolean b){
//		if(b){flags=(flags|semiperfectMask);}
//		else{flags=(flags&~semiperfectMask);}
//		assert(perfect()==b);
//		return b;
//	}
//	public boolean setSemiperfect(boolean b){
//		if(b){flags=(flags|semiperfectMask);}
//		else{flags=(flags&~semiperfectMask);}
//		assert(semiperfect()==b);
//		return b;
//	}

	public boolean plus(){return strand==Gene.PLUS;}
	public boolean minus(){return strand==Gene.MINUS;}
	public boolean perfect(){return perfect;}
	public boolean semiperfect(){return semiperfect;}
	public boolean rescued(){return rescued;}
	public byte strand(){return strand;}
	
	public final byte strand;
	public boolean rescued=false;
	public boolean perfect=false;
	public boolean semiperfect=false;
	
	public int start;
	public int stop;
	public int quickScore;
	public int score;
	public int slowScore;
	public int pairedScore;
	public int hits;
	public final int chrom;
	
	public long flags; //TODO Use this instead of fields
	
	public int[] gaps; //Limits of large gaps
	public byte[] match;
	

	public static final PositionComparator PCOMP=new PositionComparator();
	public static final long strandMask=(1L<<0);
	public static final long rescuedMask=(1L<<1);
	public static final long perfectMask=(1L<<2);
	public static final long semiperfectMask=(1L<<3);
	public static boolean verbose=false;
	
}
