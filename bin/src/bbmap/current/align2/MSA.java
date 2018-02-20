package align2;

import java.util.Arrays;

import stream.Read;
import stream.SiteScore;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;

/**
 * @author Brian Bushnell
 * @date Jun 20, 2013
 *
 */
public abstract class MSA {
	
	public static final float minIdToMinRatio(double minid, String classname){
		if("MultiStateAligner9ts".equalsIgnoreCase(classname)){
			return MultiStateAligner9ts.minIdToMinRatio(minid);
		}else if("MultiStateAligner10ts".equalsIgnoreCase(classname)){
			return MultiStateAligner10ts.minIdToMinRatio(minid);
		}else if("MultiStateAligner11ts".equalsIgnoreCase(classname)){
			return MultiStateAligner11ts.minIdToMinRatio(minid);
		}else if("MultiStateAligner9PacBio".equalsIgnoreCase(classname)){
			return MultiStateAligner9PacBio.minIdToMinRatio(minid);
		}else{
			assert(false) : "Unhandled MSA type: "+classname;
			return MultiStateAligner11ts.minIdToMinRatio(minid);
		}
	}
	
	public static final MSA makeMSA(int maxRows_, int maxColumns_, boolean colorspace_, String classname){
		if("MultiStateAligner9ts".equalsIgnoreCase(classname)){
			return new MultiStateAligner9ts(maxRows_, maxColumns_, colorspace_);
		}else if("MultiStateAligner10ts".equalsIgnoreCase(classname)){
			return new MultiStateAligner10ts(maxRows_, maxColumns_, colorspace_);
		}else if("MultiStateAligner11ts".equalsIgnoreCase(classname)){
			return new MultiStateAligner11ts(maxRows_, maxColumns_, colorspace_);
		}else if("MultiStateAligner9PacBio".equalsIgnoreCase(classname)){
			return new MultiStateAligner9PacBio(maxRows_, maxColumns_, colorspace_);
		}else{
			assert(false) : "Unhandled MSA type: "+classname;
			return new MultiStateAligner11ts(maxRows_, maxColumns_, colorspace_);
		}
	}
	
	public MSA(int maxRows_, int maxColumns_, boolean colorspace_){
		maxRows=maxRows_;
		maxColumns=maxColumns_;
		colorspace=colorspace_;
		
		
//		packed=new int[3][maxRows+1][maxColumns+1];
//		grefbuffer=new byte[maxColumns+2];
//
//		vertLimit=new int[maxRows+1];
//		horizLimit=new int[maxColumns+1];
//		Arrays.fill(vertLimit, BADoff);
//		Arrays.fill(horizLimit, BADoff);
//		
////		for(int i=0; i<maxColumns+1; i++){
////			scores[0][i]=0-i;
////		}
//		
//		for(int matrix=0; matrix<packed.length; matrix++){
//			for(int i=1; i<=maxRows; i++){
//				for(int j=0; j<packed[matrix][i].length; j++){
//					packed[matrix][i][j]|=BADoff;
//				}
////				packed[matrix][i][0]|=MODE_INS;
//			}
////			for(int i=0; i<maxRows+1; i++){
////				scores[matrix][i][0]=(i*POINTSoff_NOREF);
////			}
//			for(int i=0; i<=maxRows; i++){
//				
//				int prevScore=(i<2 ? 0 : packed[matrix][i-1][0]);
//				int score=(i<2 ? (i*POINTSoff_INS) : 
//					(i<LIMIT_FOR_COST_3 ? prevScore+POINTSoff_INS2 : 
//						(i<LIMIT_FOR_COST_4 ? prevScore+POINTSoff_INS3 : prevScore+POINTSoff_INS4)));
//				
//				packed[matrix][i][0]=score;
//			}
////			for(int i=1; i<maxColumns+1; i++){
////				prevState[matrix][0][i]=MODE_DEL;
////			}
////			for(int i=0; i<=maxColumns; i++){
////				packed[matrix][0][i]|=MODE_MS;
////			}
//		}
	}
	
	
	/** return new int[] {rows, maxC, maxS, max}; 
	 * Will not fill areas that cannot match minScore */
	public abstract int[] fillLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore, int[] gaps);
	
	
	/** return new int[] {rows, maxC, maxS, max}; 
	 * Will not fill areas that cannot match minScore */
	public abstract int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int[] gaps);
	
	@Deprecated
	/** return new int[] {rows, maxC, maxS, max}; */
	public abstract int[] fillQ(byte[] read, byte[] ref, byte[] baseScores, int refStartLoc, int refEndLoc);

	
	/** @return {score, bestRefStart, bestRefStop} */
	/** Generates the match string */
	public abstract byte[] traceback(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state, boolean gapped);
	
	
	/** Generates the match string */
	public abstract byte[] traceback2(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state);
	
	/** @return {score, bestRefStart, bestRefStop} */
	public abstract int[] score(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc, 
			final int maxRow, final int maxCol, final int maxState, boolean gapped);
	
	/** @return {score, bestRefStart, bestRefStop}, or {score, bestRefStart, bestRefStop, padLeft, padRight} if more padding is needed */
	public abstract int[] score2(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc, 
			final int maxRow, final int maxCol, final int maxState);
	
	
	/** Will not fill areas that cannot match minScore.
	 * @return {score, bestRefStart, bestRefStop}  */
	public final int[] fillAndScoreLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore, int[] gaps){
		int a=Tools.max(0, refStartLoc);
		int b=Tools.min(ref.length-1, refEndLoc);
		assert(b>=a);
		
		int[] score;
		
		if(gaps==null){
			if(verbose){
				System.err.println("no gaps");
			}
			if(b-a>=maxColumns){
				System.err.println("Warning: Max alignment columns exceeded; restricting range. "+(b-a+1)+" > "+maxColumns);
				assert(false) : refStartLoc+", "+refEndLoc;
				b=Tools.min(ref.length-1, a+maxColumns-1);
			}
			int[] max=fillLimited(read, ref, a, b, minScore, gaps);
			score=(max==null ? null : score(read, ref, a, b, max[0], max[1], max[2], false));
		}else{
			if(verbose){System.err.println("\ngaps: "+Arrays.toString(gaps)+"\n"+new String(read)+"\ncoords: "+refStartLoc+", "+refEndLoc);}
			int[] max=fillLimited(read, ref, a, b, minScore, gaps);
			if(verbose){System.err.println("max: "+Arrays.toString(max));}
//			score=(max==null ? null : score(read, grefbuffer, 0, greflimit, max[0], max[1], max[2], true));
			score=(max==null ? null : score(read, ref, a, b, max[0], max[1], max[2], true));
		}
		return score;
	}
	
	public final int[] fillAndScoreLimited(byte[] read, SiteScore ss, int thresh, int minScore){
		return fillAndScoreLimited(read, ss.chrom, ss.start, ss.stop, thresh, minScore, ss.gaps);
	}
	
//	public final int[] translateScoreFromGappedCoordinate(int[] score)
	
	public final int[] fillAndScoreLimited(byte[] read, int chrom, int start, int stop, int thresh, int minScore, int[] gaps){
		return fillAndScoreLimited(read, Data.getChromosome(chrom).array, start-thresh, stop+thresh, minScore, gaps);
	}
	
	@Deprecated
	public final int[] fillAndScoreQ(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, byte[] baseScores){
		int a=Tools.max(0, refStartLoc);
		int b=Tools.min(ref.length-1, refEndLoc);
		assert(b>=a);
		if(b-a>=maxColumns){
			System.err.println("Warning: Max alignment columns exceeded; restricting range. "+(b-a+1)+" > "+maxColumns);
			b=Tools.min(ref.length-1, a+maxColumns-1);
		}
		int[] max=fillQ(read, ref, baseScores, a, b);
//		int[] score=score(read, ref,  a, b, max[0], max[1], max[2]);
//		return score;
		return null;
	}
	
	@Deprecated
	public final int[] fillAndScoreQ(byte[] read, SiteScore ss, int thresh, byte[] baseScores){
		return fillAndScoreQ(read, ss.chrom, ss.start, ss.stop, thresh, baseScores);
	}
	
	@Deprecated
	public final int[] fillAndScoreQ(byte[] read, int chrom, int start, int stop, int thresh, byte[] baseScores){
		return fillAndScoreQ(read, Data.getChromosome(chrom).array, start-thresh, stop+thresh, baseScores);
	}
	
//	public final int scoreNoIndels(byte[] read, SiteScore ss){
//		
//		ChromosomeArray cha=Data.getChromosome(ss.chrom);
//		final int refStart=ss.start;
//		
//		int score=0;
//		int mode=MODE_START;
//		int timeInMode=0;
//		if(refStart<0 || refStart+read.length>cha.maxIndex+1){return -99999;} //TODO: Partial match
//		
//		for(int i=0; i<read.length; i++){
//			byte c=read[i];
//			byte r=cha.get(refStart+i);
//			
//			if(c==r){
//				if(mode==MODE_MS){
//					timeInMode++;
//					score+=POINTSoff_MATCH2;
//				}else{
//					timeInMode=0;
//					score+=POINTSoff_MATCH;
//				}
//				mode=MODE_MS;
//			}else if(c<0 || c=='N'){
//				score+=POINTSoff_NOCALL;
//			}else if(r<0 || r=='N'){
//				score+=POINTSoff_NOREF;
//			}else{
//				if(mode==MODE_SUB){timeInMode++;}
//				else{timeInMode=0;}
//				
//				if(timeInMode==0){score+=POINTSoff_SUB;}
//				else if(timeInMode<LIMIT_FOR_COST_3){score+=POINTSoff_SUB2;}
//				else{score+=POINTSoff_SUB3;}
//			}
//		}
//		
//		return score;
//	}
	

	public final int scoreNoIndels(byte[] read, SiteScore ss){
		ChromosomeArray cha=Data.getChromosome(ss.chrom);
		return scoreNoIndels(read, cha.array, ss.start);
	}

	public final int scoreNoIndels(byte[] read, final int chrom, final int refStart){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return scoreNoIndels(read, cha.array, refStart);
	}
	
	public final int scoreNoIndels(byte[] read, SiteScore ss, byte[] baseScores){
		ChromosomeArray cha=Data.getChromosome(ss.chrom);
		return scoreNoIndels(read, cha.array, baseScores, ss.start);
	}

	public final int scoreNoIndels(byte[] read, final int chrom, final int refStart, byte[] baseScores){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return scoreNoIndels(read, cha.array, baseScores, refStart);
	}

//	public final int scoreNoIndels(byte[] read, final int chrom, final int refStart){
	
	/** Calculates score based on an array from Index. */
	public abstract int calcAffineScore(int[] locArray, byte[] baseScores, byte[] bases);

	/** Calculates score based on an array from Index using a kfilter.  Slightly slower. */
	public abstract int calcAffineScore(int[] locArray, byte[] baseScores, byte[] bases, int minContig);

	public abstract int scoreNoIndels(byte[] read, byte[] ref, final int refStart);
	public int scoreNoIndels(byte[] read, byte[] ref, final int refStart, final SiteScore ss){
		throw new RuntimeException("Unimplemented method in class "+this.getClass());
	}

	public abstract byte[] genMatchNoIndels(byte[] read, byte[] ref, final int refStart);

	public abstract int scoreNoIndels(byte[] read, byte[] ref, byte[] baseScores, final int refStart);
	public int scoreNoIndels(byte[] read, byte[] ref, byte[] baseScores, final int refStart, SiteScore ss){
		throw new RuntimeException("Unimplemented method in class "+this.getClass());
	}
	
	public abstract int scoreNoIndelsAndMakeMatchString(byte[] read, byte[] ref, byte[] baseScores, final int refStart, byte[][] matchReturn);
	
	public abstract int scoreNoIndelsAndMakeMatchString(byte[] read, byte[] ref, final int refStart, byte[][] matchReturn);

//	public final boolean toLocalAlignment(Read r, SiteScore ss, int minToClip){return toLocalAlignment(r, ss, minToClip, 1);}
	
	/** Assumes match string is in long format */
	public final boolean toLocalAlignment(Read r, SiteScore ss, int minToClip, float matchPointsMult){
		final byte[] match=r.match, bases=(r.strand()==Gene.PLUS ? r.bases : AminoAcid.reverseComplementBases(r.bases));
		if(match==null || match.length<1){return false;}
		
		assert(match==ss.match);
		assert(match==r.match);
		
		if(r.containsXY2()){
			if(verbose){System.err.println("\nInitial0:");}
			if(verbose){System.err.println("0: match="+new String(match));}
			if(verbose){System.err.println("0: r.start="+r.start+", r.stop="+r.stop+"; len="+bases.length+"; reflen="+(r.stop-r.start+1));}
			ss.fixXY(bases, false, this);
			assert(match==ss.match);
			assert(match==r.match);
		}
		
		int maxScore=-1;
		
		int startLocC=-1;
		int stopLocC=-1;
		int lastZeroC=0;
		
		int startLocM=-1;
		int stopLocM=-1;
		int lastZeroM=0;
		
		int startLocR=-1;
		int stopLocR=-1;
		int lastZeroR=0;
		
		byte mode=match[0], prevMode='0';
		int current=0, prevStreak=0;
		int cpos=0;
		int rpos=r.start;
		int score=0;

		if(verbose){System.err.println("\nInitial:");}
		if(verbose){System.err.println("A: match="+new String(match));}
		if(verbose){System.err.println("A: r.start="+r.start+", r.stop="+r.stop+"; rpos="+rpos+"; len="+bases.length+"; reflen="+(r.stop-r.start+1));}
		if(verbose){System.err.println(new String(bases));}
		if(verbose){System.err.println(Data.getChromosome(r.chrom).getString(r.start, r.stop));}
		
		if(verbose){
			int calcscore=score(match);
			System.err.println("A: score="+r.mapScore+", ss.slowScore="+ss.slowScore+", calcscore="+calcscore);
//			assert(ss.slowScore<=calcscore); //May be lower due to ambig3.  I found a case where this line fails, possibly due to long deletions?
		}
		
		for(int mpos=0; mpos<match.length; mpos++){
			byte c=match[mpos];
			
			if(mode==c){
				current++;
			}else{
				if(mode=='m'){
					if(score<=0){
						score=0;
						lastZeroC=cpos;
						lastZeroM=mpos-current;
						lastZeroR=rpos;
					}
					int add=calcMatchScore(current);
					score+=(matchPointsMult*add);
//					if(prevMode=='N' || prevMode=='R'){score=score+POINTS_MATCH2()-POINTS_MATCH();} //Don't penalize first match after N
					cpos+=current;
					rpos+=current;
					if(score>maxScore){
						maxScore=score;
						startLocC=lastZeroC;
						startLocM=lastZeroM;
						startLocR=lastZeroR;
						stopLocC=cpos-1;
						stopLocM=mpos-1;
						stopLocR=rpos-1;
					}
				}else if(mode=='S'){
					score+=calcSubScore(current);
					if(prevMode=='N' || prevMode=='R'){score=score+POINTS_SUB2()-POINTS_SUB();} //Don't penalize first sub after N
					else if(prevMode=='m' && prevStreak<2){score=score+POINTS_SUBR()-POINTS_SUB();}
					cpos+=current;
					rpos+=current;
				}else if(mode=='D'){
					score+=calcDelScore(current, true);
					rpos+=current;
				}else if(mode=='I'){
					score+=calcInsScore(current);
					cpos+=current;
				}else if(mode=='C'){
					cpos+=current;
				}else if(mode=='X' || mode=='Y'){
					score+=calcInsScore(current);//TODO: Consider changing XY to subs
					cpos+=current;
					rpos+=current;
				}else if(mode=='N'){
					score+=calcNocallScore(current);
					cpos+=current;
					rpos+=current;
				}else if(mode=='R'){
					score+=calcNorefScore(current);
					cpos+=current;
					rpos+=current;
				}else{
					assert(false) : "Unhandled symbol "+mode+"\n"+(char)mode+"\n"+new String(match)+"\n"+new String(bases);
				}
				if(verbose){System.err.println("mode "+(char)mode+"->"+(char)c+"; rpos="+rpos);}
				prevMode=mode;
				prevStreak=current;
				mode=c;
				current=1;
			}
		}
		if(current>0){
			assert(mode==match[match.length-1]);
			if(mode=='m'){
				if(score<=0){
					score=0;
					lastZeroC=cpos;
					lastZeroM=match.length-current;
					lastZeroR=rpos;
				}
				int add=calcMatchScore(current);
				score+=(matchPointsMult*add);
//				if(prevMode=='N' || prevMode=='R'){score=score+POINTS_MATCH2()-POINTS_MATCH();} //Don't penalize first match after N
				cpos+=current;
				rpos+=current;
				if(score>maxScore){
					maxScore=score;
					startLocC=lastZeroC;
					startLocM=lastZeroM;
					startLocR=lastZeroR;
					stopLocC=cpos-1;
					stopLocM=match.length-1;
					stopLocR=rpos-1;
				}
			}else if(mode=='S'){
				score+=calcSubScore(current);
				if(prevMode=='N' || prevMode=='R'){score=score+POINTS_SUB2()-POINTS_SUB();} //Don't penalize first sub after N
				else if(prevMode=='m' && prevStreak<2){score=score+POINTS_SUBR()-POINTS_SUB();}
				cpos+=current;
				rpos+=current;
			}else if(mode=='D'){
				score+=calcDelScore(current, true);
				rpos+=current;
			}else if(mode=='I'){
				score+=calcInsScore(current);
				cpos+=current;
			}else if(mode=='C'){
				cpos+=current;
			}else if(mode=='X' || mode=='Y'){
				score+=calcInsScore(current);
				cpos+=current;
			}else if(mode=='N'){
				score+=calcNocallScore(current);
				cpos+=current;
				rpos+=current;
			}else if(mode=='R'){
				score+=calcNorefScore(current);
				cpos+=current;
				rpos+=current;
			}else if(mode!=0){
				assert(false) : "Unhandled symbol "+mode+"\n"+(char)mode+"\n"+new String(match)+"\n"+new String(bases);
			}
			if(verbose){System.err.println("mode "+(char)mode+"->end; rpos="+rpos);}
		}
		
		if(startLocC<0 || stopLocC<0){
			assert(false) : "Failed.";
			return false;
		}
		

		if(verbose){System.err.println("A: r.start="+r.start+", r.stop="+r.stop+"; rpos="+rpos+"; len="+bases.length+"; reflen="+(r.stop-r.start+1));}
		
		assert(rpos==r.stop+1) : rpos+"!="+r.start+"\n"+r;
		
		if(verbose){System.err.println("B: rpos="+rpos+", startLocR="+startLocR+", stopLocR="+stopLocR);}
		
		int headTrimR=startLocC;
		int headTrimM=startLocM;
		int tailTrimR=bases.length-stopLocC-1;
		int tailTrimM=match.length-stopLocM-1;
		
		if(verbose){System.err.println("C: headTrimR="+headTrimR+", headTrimM="+headTrimM+", tailTrimR="+tailTrimR+", tailTrimM="+tailTrimM);}
		
		if(headTrimR<=minToClip && headTrimM<=minToClip){
			headTrimR=headTrimM=0;
		}
		if(tailTrimR<=minToClip && tailTrimM<=minToClip){
			tailTrimR=tailTrimM=0;
		}
		if(headTrimR==0 && headTrimM==0 && tailTrimR==0 && tailTrimM==0){
			return false;
		}
		//Do trimming
		final int headDelta=headTrimR-headTrimM;
		final int tailDelta=tailTrimR-tailTrimM;
		final byte[] match2;
		
		if(verbose){System.err.println("D: headTrimR="+headTrimR+", headTrimM="+headTrimM+", tailTrimR="+tailTrimR+", tailTrimM="+tailTrimM);}
		if(verbose){System.err.println("D: headDelta="+headDelta+", tailDelta="+tailDelta);}
		
		if(headDelta==0 && tailDelta==0){
			//Length-neutral trimming
			match2=match;
			for(int i=0; i<headTrimM; i++){match[i]='C';}
			for(int i=match.length-tailTrimM; i<match.length; i++){match[i]='C';}
		}else{
			final int newlen=match.length-headTrimM-tailTrimM+headTrimR+tailTrimR;
			match2=new byte[newlen];
			for(int i=0; i<headTrimR; i++){match2[i]='C';}
			for(int i=match2.length-tailTrimR; i<match2.length; i++){match2[i]='C';}
			for(int i=headTrimM, i2=headTrimR, lim=match2.length-tailTrimR; i2<lim; i++, i2++){
				match2[i2]=match[i];
			}
		}
		
		assert(ss==null || ((ss.start==r.start) && (ss.stop==r.stop) && (ss.strand==r.strand()) && (ss.chrom==r.chrom) && (ss.match==r.match))) : 
			"\nr="+r+"\nr2="+r.mate+"\nss=\n"+ss+"\n"+(ss==null ? "ss is null" : ((ss.start==r.start)+", "+(ss.stop==r.stop)+", "+
			(ss.strand==r.strand())+", "+(ss.chrom==r.chrom)+", "+(ss.match==r.match)));
		
		if(headTrimR!=0){r.start=startLocR-headTrimR;}
		if(tailTrimR!=0){r.stop=stopLocR+tailTrimR;}
		r.match=match2;
		
		if(matchPointsMult!=1f){
			maxScore=score(match);
		}
		if(ss!=null){maxScore=Tools.max(maxScore, ss.slowScore);}
		r.mapScore=maxScore;

		if(verbose){System.err.println("E: r.start="+r.start+", r.stop="+r.stop);}
		
		if(ss!=null){
			assert(maxScore>=ss.slowScore) : maxScore+", "+ss.slowScore+"\n"+r.toFastq();
			ss.start=r.start;
			ss.stop=r.stop;
			ss.match=r.match;
			int pairedScore=ss.pairedScore>0 ? Tools.max(ss.pairedScore+(maxScore-ss.slowScore), 0) : 0;
		}
		
		return true;
	}
	

	/** Assumes match string is in long format. */
	public final int score(byte[] match){
		if(match==null || match.length<1){return 0;}
		
		byte mode=match[0], prevMode='0';
		int current=0, prevStreak=0;
		int score=0;
		
		for(int mpos=0; mpos<match.length; mpos++){
			byte c=match[mpos];
			
			if(mode==c){
				current++;
			}else{
				if(mode=='m'){
					score+=calcMatchScore(current);
//					if(prevMode=='N' || prevMode=='R'){score=score+POINTS_MATCH2()-POINTS_MATCH();} //Don't penalize first match after N
				}else if(mode=='S'){
					score+=calcSubScore(current);
					if(prevMode=='N' || prevMode=='R'){score=score+POINTS_SUB2()-POINTS_SUB();} //Don't penalize first sub after N
					else if(prevMode=='m' && prevStreak<2){score=score+POINTS_SUBR()-POINTS_SUB();}
				}else if(mode=='D'){
					score+=calcDelScore(current, true);
				}else if(mode=='I'){
					score+=calcInsScore(current);
				}else if(mode=='C'){
					//do nothing
				}else if(mode=='X' || mode=='Y'){
					score+=calcInsScore(current);
				}else if(mode=='N'){
					score+=calcNocallScore(current);
				}else if(mode=='R'){
					score+=calcNorefScore(current);
				}else{
					assert(false) : "Unhandled symbol "+mode+"\n"+(char)mode+"\n"+new String(match);
				}
				if(verbose){System.err.println("mode "+(char)mode+"->"+(char)c+"\tcurrent="+current+"\tscore="+score);}
				prevMode=mode;
				prevStreak=current;
				mode=c;
				current=1;
			}
		}
		if(current>0){
			assert(mode==match[match.length-1]);
			if(mode=='m'){
				score+=calcMatchScore(current);
//				if(prevMode=='N' || prevMode=='R'){score=score+POINTS_MATCH2()-POINTS_MATCH();} //Don't penalize first match after N
			}else if(mode=='S'){
				score+=calcSubScore(current);
				if(prevMode=='N' || prevMode=='R'){score=score+POINTS_SUB2()-POINTS_SUB();} //Don't penalize first sub after N
				else if(prevMode=='m' && prevStreak<2){score=score+POINTS_SUBR()-POINTS_SUB();}
			}else if(mode=='D'){
				score+=calcDelScore(current, true);
			}else if(mode=='I'){
				score+=calcInsScore(current);
			}else if(mode=='C'){
				//do nothing
			}else if(mode=='X' || mode=='Y'){
				score+=calcInsScore(current);
			}else if(mode=='N'){
				score+=calcNocallScore(current);
			}else if(mode=='R'){
				score+=calcNorefScore(current);
			}else if(mode!=0){
				assert(false) : "Unhandled symbol "+mode+"\n"+(char)mode+"\n"+new String(match);
			}
			if(verbose){System.err.println("mode "+(char)mode+"->end; score="+score);}
		}
		
		return score;
	}
	
//	//TODO
//	public final byte[] softClipBoundsShortmatch(byte[] match, byte[] bases, int minToClip){
//		if(match==null || match.length<1){return null;}
//		int[] score=new int[bases.length];
//
//		byte mode='0', c='0';
//		int current=0;
//		int rpos=0;
//		long currentScore;
//		for(int i=0; i<match.length; i++){
//			c=match[i];
//			if(Character.isDigit(c)){
//				current=(current*10)+(c-'0');
//			}else{
//				if(mode==c){
//					current=Tools.max(current+1, 2);
//				}else{
//					current=Tools.max(current, 1);
//
//					if(mode=='m'){
//						msdicn[0]+=current;
//					}else if(mode=='S'){
//						msdicn[1]+=current;
//					}else if(mode=='D'){
//						msdicn[2]+=current;
//					}else if(mode=='I'){
//						msdicn[3]+=current;
//					}else if(mode=='C' || mode=='X' || mode=='Y'){
//						msdicn[4]+=current;
//					}else if(mode=='N' || mode=='R'){
//						msdicn[5]+=current;
//					}
//					mode=c;
//					current=0;
//				}
//			}
//		}
//		if(current>0 || !Character.isDigit(c)){
//			current=Tools.max(current, 1);
//			if(mode=='m'){
//				msdicn[0]+=current;
//			}else if(mode=='S'){
//				msdicn[1]+=current;
//			}else if(mode=='D'){
//				msdicn[2]+=current;
//			}else if(mode=='I'){
//				msdicn[3]+=current;
//			}else if(mode=='C' || mode=='X' || mode=='Y'){
//				msdicn[4]+=current;
//			}else if(mode=='N' || mode=='R'){
//				msdicn[5]+=current;
//			}
//		}
//		return msdicn;
//	}
	
	public abstract int maxQuality(int numBases);
	
	public abstract int maxQuality(byte[] baseScores);
	
	public abstract int maxImperfectScore(int numBases);
	
	public abstract int maxImperfectScore(byte[] baseScores);
	
	public final static String toString(int[] a){
		
		int width=7;
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num : a){
			String s=" "+num;
			int spaces=width-s.length();
			assert(spaces>=0) : width+", "+s.length()+", "+s+", "+num+", "+spaces;
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	public final static String toTimePacked(int[] a, int TIMEMASK){
		int width=7;
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num_ : a){
			int num=num_&TIMEMASK;
			String s=" "+num;
			int spaces=width-s.length();
			assert(spaces>=0) : width+", "+s.length()+", "+s+", "+num+", "+spaces;
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	public final static String toScorePacked(int[] a, int SCOREOFFSET){
		int width=7;

		String minString=" -";
		String maxString="  ";
		while(minString.length()<width){minString+='9';}
		while(maxString.length()<width){maxString+='9';}
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num_ : a){
			int num=num_>>SCOREOFFSET;
			String s=" "+num;
			if(s.length()>width){s=num>0 ? maxString : minString;}
			int spaces=width-s.length();
			assert(spaces>=0) : width+", "+s.length()+", "+s+", "+num+", "+spaces;
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	public final static String toString(byte[] a){
		
		int width=6;
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num : a){
			String s=" "+num;
			int spaces=width-s.length();
			assert(spaces>=0);
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	public final static String toString(byte[] ref, int startLoc, int stopLoc){
		StringBuilder sb=new StringBuilder(stopLoc-startLoc+1);
		for(int i=startLoc; i<=stopLoc; i++){sb.append((char)ref[i]);}
		return sb.toString();
	}
	
	public final int calcMatchScore(int len){
		assert(len>0) : len;
		return POINTS_MATCH()+(len-1)*POINTS_MATCH2();
	}
	
	public final int calcSubScore(int len){
		assert(len>0) : len;
		final int lim3=LIMIT_FOR_COST_3();
		int score=POINTS_SUB();
		if(len>lim3){
			score+=(len-lim3)*POINTS_SUB3();
			len=lim3;
		}
		if(len>1){
			score+=(len-1)*POINTS_SUB2();
		}
		return score;
	}
	
	public final int calcNorefScore(int len){return len*POINTS_NOREF();}
	
	public final int calcNocallScore(int len){return len*POINTS_NOCALL();}
	
	public abstract int calcDelScore(int len, boolean approximateGaps);
	
//	private static int calcDelScoreOffset(int len){
//		if(len<=0){return 0;}
//		int score=POINTSoff_DEL;
//		
//		if(len>LIMIT_FOR_COST_5){
//			score+=((len-LIMIT_FOR_COST_5+MASK5)/TIMESLIP)*POINTSoff_DEL5;
//			len=LIMIT_FOR_COST_5;
//		}
//		if(len>LIMIT_FOR_COST_4){
//			score+=(len-LIMIT_FOR_COST_4)*POINTSoff_DEL4;
//			len=LIMIT_FOR_COST_4;
//		}
//		if(len>LIMIT_FOR_COST_3){
//			score+=(len-LIMIT_FOR_COST_3)*POINTSoff_DEL3;
//			len=LIMIT_FOR_COST_3;
//		}
//		if(len>1){
//			score+=(len-1)*POINTSoff_DEL2;
//		}
//		return score;
//	}
	
	public abstract int calcInsScore(int len);
	
//	private static int calcInsScoreOffset(int len){
//		if(len<=0){return 0;}
//		int score=POINTSoff_INS;
//		if(len>LIMIT_FOR_COST_4){
//			score+=(len-LIMIT_FOR_COST_4)*POINTSoff_INS4;
//			len=LIMIT_FOR_COST_4;
//		}
//		if(len>LIMIT_FOR_COST_3){
//			score+=(len-LIMIT_FOR_COST_3)*POINTSoff_INS3;
//			len=LIMIT_FOR_COST_3;
//		}
//		if(len>1){
//			score+=(len-1)*POINTSoff_INS2;
//		}
//		return score;
//	}
	
	static final int GAPBUFFER=Shared.GAPBUFFER;
	static final int GAPBUFFER2=Shared.GAPBUFFER2;
	static final int GAPLEN=Shared.GAPLEN;
	static final int MINGAP=Shared.MINGAP;
	static final int GAPCOST=Shared.GAPCOST;
	static final byte GAPC=Shared.GAPC;
	
	/** Seemingly to clear out prior data from the gref.  Not sure what else it's used for. */
	static final int GREFLIMIT2_CUSHION=128; //Tools.max(GAPBUFFER2, GAPLEN);
	
	
	/**DO NOT MODIFY*/
	public abstract byte[] getGrefbuffer();

//	public final int[] vertLimit;
//	public final int[] horizLimit;

	public abstract CharSequence showVertLimit();
	public abstract CharSequence showHorizLimit();

////	public static final int MODEBITS=2;
//	public static final int TIMEBITS=11;
//	public static final int SCOREBITS=32-TIMEBITS;
//	public static final int MAX_TIME=((1<<TIMEBITS)-1);
//	public static final int MAX_SCORE=((1<<(SCOREBITS-1))-1)-2000;
//	public static final int MIN_SCORE=0-MAX_SCORE; //Keeps it 1 point above "BAD".
//
////	public static final int MODEOFFSET=0; //Always zero.
////	public static final int TIMEOFFSET=0;
	public abstract int SCOREOFFSET();
//
////	public static final int MODEMASK=~((-1)<<MODEBITS);
////	public static final int TIMEMASK=(~((-1)<<TIMEBITS))<<TIMEOFFSET;
//	public static final int TIMEMASK=~((-1)<<TIMEBITS);
//	public static final int SCOREMASK=(~((-1)<<SCOREBITS))<<SCOREOFFSET;
	
	static final byte MODE_MS=0;
	static final byte MODE_DEL=1;
	static final byte MODE_INS=2;
	static final byte MODE_SUB=3;
	
	public abstract int POINTS_NOREF();
	public abstract int POINTS_NOCALL();
	public abstract int POINTS_MATCH();
	public abstract int POINTS_MATCH2();
	public abstract int POINTS_COMPATIBLE();
	public abstract int POINTS_SUB();
	public abstract int POINTS_SUBR();
	public abstract int POINTS_SUB2();
	public abstract int POINTS_SUB3();
	public abstract int POINTS_MATCHSUB();
	public abstract int POINTS_INS();
	public abstract int POINTS_INS2();
	public abstract int POINTS_INS3();
	public abstract int POINTS_INS4();
	public abstract int POINTS_DEL();
	public abstract int POINTS_DEL2();
	public abstract int POINTS_DEL3();
	public abstract int POINTS_DEL4();
	public abstract int POINTS_DEL5();
	public abstract int POINTS_DEL_REF_N();
	public abstract int POINTS_GAP();

	public abstract int TIMESLIP();
	public abstract int MASK5();
	
	abstract int BARRIER_I1();
	abstract int BARRIER_D1();

	public abstract int LIMIT_FOR_COST_3();
	public abstract int LIMIT_FOR_COST_4();
	public abstract int LIMIT_FOR_COST_5();
	
	public abstract int BAD();
	
	
//	public static final int POINTSoff_NOREF=(POINTS_NOREF<<SCOREOFFSET);
//	public static final int POINTSoff_NOCALL=(POINTS_NOCALL<<SCOREOFFSET);
//	public static final int POINTSoff_MATCH=(POINTS_MATCH<<SCOREOFFSET);
//	public static final int POINTSoff_MATCH2=(POINTS_MATCH2<<SCOREOFFSET);
//	public static final int POINTSoff_COMPATIBLE=(POINTS_COMPATIBLE<<SCOREOFFSET);
//	public static final int POINTSoff_SUB=(POINTS_SUB<<SCOREOFFSET);
//	public static final int POINTSoff_SUBR=(POINTS_SUBR<<SCOREOFFSET);
//	public static final int POINTSoff_SUB2=(POINTS_SUB2<<SCOREOFFSET);
//	public static final int POINTSoff_SUB3=(POINTS_SUB3<<SCOREOFFSET);
//	public static final int POINTSoff_MATCHSUB=(POINTS_MATCHSUB<<SCOREOFFSET);
//	public static final int POINTSoff_INS=(POINTS_INS<<SCOREOFFSET);
//	public static final int POINTSoff_INS2=(POINTS_INS2<<SCOREOFFSET);
//	public static final int POINTSoff_INS3=(POINTS_INS3<<SCOREOFFSET);
//	public static final int POINTSoff_INS4=(POINTS_INS4<<SCOREOFFSET);
//	public static final int POINTSoff_DEL=(POINTS_DEL<<SCOREOFFSET);
//	public static final int POINTSoff_DEL2=(POINTS_DEL2<<SCOREOFFSET);
//	public static final int POINTSoff_DEL3=(POINTS_DEL3<<SCOREOFFSET);
//	public static final int POINTSoff_DEL4=(POINTS_DEL4<<SCOREOFFSET);
//	public static final int POINTSoff_DEL5=(POINTS_DEL5<<SCOREOFFSET);
//	public static final int POINTSoff_GAP=(POINTS_GAP<<SCOREOFFSET);
//	public static final int POINTSoff_DEL_REF_N=(POINTS_DEL_REF_N<<SCOREOFFSET);
//	public static final int BADoff=(BAD<<SCOREOFFSET);
//	public static final int MAXoff_SCORE=MAX_SCORE<<SCOREOFFSET;
//	public static final int MINoff_SCORE=MIN_SCORE<<SCOREOFFSET;
	
	
	public final int maxRows;
	public final int maxColumns;
	
	public final boolean colorspace;

	public long iterationsLimited=0;
	public long iterationsUnlimited=0;

	public boolean verbose=false;
	public boolean verbose2=false;

	public static int bandwidth=0;
	public static float bandwidthRatio=0;

}
