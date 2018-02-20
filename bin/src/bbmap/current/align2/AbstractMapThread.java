package align2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import stream.ConcurrentReadStreamInterface;
import stream.RTextOutputStream3;
import stream.Read;
import stream.SiteScore;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;



/**
 * @author Brian Bushnell
 * @date Feb 27, 2013
 *
 */
public abstract class AbstractMapThread extends Thread {
	
	AbstractMapThread(ConcurrentReadStreamInterface cris_, 
			RTextOutputStream3 outStream_, RTextOutputStream3 outStreamMapped_, RTextOutputStream3 outStreamUnmapped_, RTextOutputStream3 outStreamBlack_,
			boolean colorspace_, boolean SLOW_ALIGN_, boolean LOCAL_ALIGN_, boolean AMBIGUOUS_TOSS_, 
			boolean AMBIGUOUS_RANDOM_, boolean AMBIGUOUS_ALL_, boolean TRIM_LEFT_, boolean TRIM_RIGHT_, boolean UNTRIM_, byte TRIM_QUAL_, int THRESH_, 
			int minChrom_, int maxChrom_, int KFILTER_, boolean KILL_BAD_PAIRS_, boolean SAVE_AMBIGUOUS_XY_,
			boolean translateToBaseSpace_, boolean REQUIRE_CORRECT_STRANDS_PAIRS_,
			boolean SAME_STRAND_PAIRS_, boolean DO_RESCUE_, boolean STRICT_MAX_INDEL_, int SLOW_ALIGN_PADDING_, int SLOW_RESCUE_PADDING_,
			String MSA_TYPE_, int keylen_, boolean PERFECTMODE_, boolean SEMIPERFECTMODE_, boolean FORBID_SELF_MAPPING_, boolean RCOMP_MATE_, 
			boolean MAKE_MATCH_STRING_, boolean DONT_OUTPUT_UNMAPPED_READS_, boolean DONT_OUTPUT_BLACKLISTED_READS_, boolean PRINT_SECONDARY_ALIGNMENTS_, 
			boolean QUICK_MATCH_STRINGS_, int MAX_SITESCORES_TO_PRINT_, float MINIMUM_ALIGNMENT_SCORE_RATIO_,
			float keyDensity_, float maxKeyDensity_, float minKeyDensity_, int maxDesiredKeys_,
			int MIN_APPROX_HITS_TO_KEEP_, boolean USE_EXTENDED_SCORE_, int BASE_HIT_SCORE_, boolean USE_AFFINE_SCORE_, int MAX_INDEL_,
			boolean TRIM_LIST_, int TIP_DELETION_SEARCH_RANGE_){
		
		
		cris=cris_;
		outStream=outStream_;
		outStreamMapped=outStreamMapped_;
		outStreamUnmapped=outStreamUnmapped_;
		outStreamBlack=outStreamBlack_;
		
		colorspace=colorspace_;
		SLOW_ALIGN=SLOW_ALIGN_;
		LOCAL_ALIGN=LOCAL_ALIGN_;
		AMBIGUOUS_TOSS=AMBIGUOUS_TOSS_;
		AMBIGUOUS_RANDOM=AMBIGUOUS_RANDOM_;
		AMBIGUOUS_ALL=AMBIGUOUS_ALL_;		
		TRIM_LEFT=TRIM_LEFT_;
		TRIM_RIGHT=TRIM_RIGHT_;
		UNTRIM=UNTRIM_;
		TRIM_QUAL=TRIM_QUAL_;
		THRESH=THRESH_;
		minChrom=minChrom_;
		maxChrom=maxChrom_;
		KFILTER=KFILTER_;
		
		KILL_BAD_PAIRS=KILL_BAD_PAIRS_;
		SAVE_AMBIGUOUS_XY=SAVE_AMBIGUOUS_XY_;
//		GEN_MATCH_FAST=GEN_MATCH_FAST_;
		translateToBaseSpace=translateToBaseSpace_;
		SLOW_ALIGN_PADDING=SLOW_ALIGN_PADDING_;
		SLOW_RESCUE_PADDING=SLOW_RESCUE_PADDING_;
		DO_RESCUE=DO_RESCUE_;
		STRICT_MAX_INDEL=STRICT_MAX_INDEL_;
		BANDWIDTH=MSA.bandwidth;
		PAIRED=cris.paired();
		REQUIRE_CORRECT_STRANDS_PAIRS=REQUIRE_CORRECT_STRANDS_PAIRS_;
		SAME_STRAND_PAIRS=SAME_STRAND_PAIRS_;
		
		/* ------------ */
		
		TRIM_LIST=TRIM_LIST_;
		TIP_DELETION_SEARCH_RANGE=TIP_DELETION_SEARCH_RANGE_;
		FIND_TIP_DELETIONS=TIP_DELETION_SEARCH_RANGE>0;
		
		MIN_APPROX_HITS_TO_KEEP=MIN_APPROX_HITS_TO_KEEP_;
		USE_EXTENDED_SCORE=USE_EXTENDED_SCORE_;
		BASE_HIT_SCORE=BASE_HIT_SCORE_;
		BASE_KEY_HIT_SCORE=BASE_HIT_SCORE*keylen_;
		USE_AFFINE_SCORE=USE_AFFINE_SCORE_;
		EXPECTED_LEN_LIMIT=(ALIGN_COLUMNS()*17)/20-(2*(SLOW_ALIGN_PADDING+10)); //TODO: Due to some bug in expected length calculation, this is low.
		MAX_INDEL=MAX_INDEL_;
		ALIGN_COLUMNS=ALIGN_COLUMNS();
		
		/* ------------ */
		
		
		KEYLEN=keylen_;
		keyDensity=keyDensity_;
		maxKeyDensity=maxKeyDensity_;
		minKeyDensity=minKeyDensity_;
		maxDesiredKeys=maxDesiredKeys_;
		
		MINIMUM_ALIGNMENT_SCORE_RATIO=MINIMUM_ALIGNMENT_SCORE_RATIO_;
		MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED=Tools.max(MINIMUM_ALIGNMENT_SCORE_RATIO*.80f, 1-((1-MINIMUM_ALIGNMENT_SCORE_RATIO)*1.4f));
		MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE=Tools.max(MINIMUM_ALIGNMENT_SCORE_RATIO*.60f,  1-((1-MINIMUM_ALIGNMENT_SCORE_RATIO)*1.8f));
//		TRIM_LIST=TRIM_LIST_;
		MAKE_MATCH_STRING=(MAKE_MATCH_STRING_ || STRICT_MAX_INDEL_);
		assert(SLOW_ALIGN_PADDING>=0);
		
		DONT_OUTPUT_UNMAPPED_READS=DONT_OUTPUT_UNMAPPED_READS_;
		DONT_OUTPUT_BLACKLISTED_READS=DONT_OUTPUT_BLACKLISTED_READS_;
		MAX_SITESCORES_TO_PRINT=MAX_SITESCORES_TO_PRINT_;
		PRINT_SECONDARY_ALIGNMENTS=PRINT_SECONDARY_ALIGNMENTS_;
		QUICK_MATCH_STRINGS=((QUICK_MATCH_STRINGS_ || STRICT_MAX_INDEL_) && MAKE_MATCH_STRING);
		
		RCOMP_MATE=RCOMP_MATE_;
		PERFECTMODE=PERFECTMODE_;
		SEMIPERFECTMODE=SEMIPERFECTMODE_;
		FORBID_SELF_MAPPING=FORBID_SELF_MAPPING_;
		assert(!(RCOMP_MATE/* || FORBID_SELF_MAPPING*/)) : "RCOMP_MATE: TODO";
		
//		TIP_DELETION_SEARCH_RANGE=TIP_DELETION_SEARCH_RANGE_;
//		FIND_TIP_DELETIONS=TIP_DELETION_SEARCH_RANGE>0;
//		EXPECTED_LEN_LIMIT=(ALIGN_COLUMNS*17)/20-(2*(SLOW_ALIGN_PADDING+10)); //TODO: Due to some bug in expected length calculation, this is low.
		MSA_TYPE=MSA_TYPE_;
		EXTRA_PADDING=(BANDWIDTH<1 && (MSA.bandwidthRatio<=0 || MSA.bandwidthRatio>=0.2f) ? 
				EXTRA_PADDING : Tools.min(EXTRA_PADDING, Tools.max(BANDWIDTH/4, (int)(MSA.bandwidthRatio*60))));
		
		if(SLOW_ALIGN || MAKE_MATCH_STRING){
			msa=MSA.makeMSA(ALIGN_ROWS(), ALIGN_COLUMNS(), colorspace, MSA_TYPE);
			POINTS_MATCH=msa.POINTS_MATCH();
			POINTS_MATCH2=msa.POINTS_MATCH2();
//			CLEARZONE1=(int)(CLEARZONE_RATIO1*POINTS_MATCH2);
//			CLEARZONE1b=(int)(CLEARZONE_RATIO1b*POINTS_MATCH2);
//			CLEARZONE1c=(int)(CLEARZONE_RATIO1c*POINTS_MATCH2);
//			CLEARZONEP=(int)(CLEARZONE_RATIOP*POINTS_MATCH2);
//			CLEARZONE3=(int)(CLEARZONE_RATIO3*POINTS_MATCH2);
			CLEARZONE1e=(int)(2*POINTS_MATCH2-POINTS_MATCH-msa.POINTS_SUB())+1;
		}else{
			POINTS_MATCH=70;
			POINTS_MATCH2=100;
			msa=null;
//			CLEARZONE1=0;
//			CLEARZONE1b=0;
//			CLEARZONE1c=0;
//			CLEARZONEP=0;
//			CLEARZONE3=0;
			CLEARZONE1e=0;
		}
		
//		CLEARZONE1b_CUTOFF_FLAT=CLEARZONE1b_CUTOFF_FLAT_RATIO*POINTS_MATCH2;
//		CLEARZONE1c_CUTOFF_FLAT=CLEARZONE1c_CUTOFF_FLAT_RATIO*POINTS_MATCH2;
//		INV_CLEARZONE3=(CLEARZONE3==0 ? 0 : 1f/CLEARZONE3);
		
		if(translateToBaseSpace){
			tcr=new TranslateColorspaceRead(MSA.makeMSA(ALIGN_ROWS(), ALIGN_COLUMNS()+500, false, MSA_TYPE));
			if(msa!=null){assert(msa.colorspace);}
		}else{
			tcr=null;
		}
		
//		index=new BBIndex(KEYLEN, minChrom, maxChrom, KFILTER, msa);
		GENERATE_KEY_SCORES_FROM_QUALITY=AbstractIndex.GENERATE_KEY_SCORES_FROM_QUALITY;
		readstats=(ReadStats.COLLECT_MATCH_STATS || ReadStats.COLLECT_QUALITY_STATS || ReadStats.COLLECT_INSERT_STATS ? new ReadStats() : null);
		
		
	}
	
	public abstract int ALIGN_COLUMNS();
	public abstract int ALIGN_ROWS();
	abstract int CLEARZONE1();
	
	abstract AbstractIndex index();
	
	
	@Override
	public final void run() {
		//System.err.println("Waiting on a list... (initial)");
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> readlist=ln.list;
		
//		long count=System.currentTimeMillis();
//		String os=System.getProperty("os.name");
//		int procs=Runtime.getRuntime().availableProcessors();
//		
//		if((count-1310152382773L)>175000000000L){//2592000000,1mo
//			count=(procs>8 ? 1 : 2)+((hashCode()&0xFFFFFFF)%5);
//		}
		final boolean black=(Blacklist.hasBlacklist());
		final boolean MAKE_QUALITY_HISTOGRAM=(readstats==null ? false : ReadStats.COLLECT_QUALITY_STATS);
		final boolean MAKE_MATCH_HISTOGRAM=(readstats==null ? false : ReadStats.COLLECT_MATCH_STATS);
		final boolean MAKE_INSERT_HISTOGRAM=(readstats==null ? false : ReadStats.COLLECT_INSERT_STATS);
		
		if(SKIP_INITIAL>0){
			while(!readlist.isEmpty()){
				
				if(readlist.get(readlist.size()-1).numericID<SKIP_INITIAL){
					//Do nothing
				}else{
					 while(readlist.get(0).numericID<SKIP_INITIAL){readlist.remove(0);}
					 break;
				}
				
				writeList(new ArrayList<Read>(1), black, ln.id);

				cris.returnList(ln, readlist.isEmpty());
//				if(count>0){
//					cris.returnList(ln, readlist.isEmpty());
//					count--;
//				}
				
				//System.err.println("Waiting on a list...");
				ln=cris.nextList();
				readlist=ln.list;
			}
		}
		
		while(!readlist.isEmpty()){

			//System.err.println("Got a list of size "+readlist.size());
			for(int i=0; i<readlist.size(); i++){
				Read r=readlist.get(i);
				assert(r.mate==null || (r.pairnum()==0 && r.mate.pairnum()==1)) : r.pairnum()+", "+r.mate.pairnum();

				//				System.out.println("Got read: "+r.toText(false));
				//				System.out.println("Synthetic: "+r.synthetic());

				assert(r.colorspace()==colorspace) : "\n"+colorspace+", "+r.colorspace()+"\n"+r.toText(false)+"\n";

				//if(r.mapLength!=r.bases.length){r.mapLength=r.bases.length;}
				//if(r.mate!=null && r.mate.mapLength!=r.mate.bases.length){r.mate.mapLength=r.mate.bases.length;}
				assert(r.mapLength==r.bases.length) : "Incorrect length: "+r.mapLength+"!="+r.bases.length;




				if(r.synthetic()){
					syntheticReads++;
					if(r.originalSite==null){r.makeOriginalSite();}
					r.clearSite();
					if(r.mate!=null){
						assert(r.mate.synthetic());
						if(r.mate.originalSite==null){r.mate.makeOriginalSite();}
						r.mate.clearSite();
					}
				}

				//Clear these in case output is being re-used
				r.clearAnswers(true);

				assert(r.bases==null || r.bases.length<=maxReadLength()) : "Read "+r.numericID+", length "+r.bases.length+", exceeds the limit of "+maxReadLength()+"\n"+
				"You can map the reads in chunks by reformatting to fasta, then mapping with the setting 'fastareadlen="+maxReadLength()+"'";
				final Read r2=r.mate;

				if(MAKE_QUALITY_HISTOGRAM){readstats.addToQualityHistogram(r);}

				if(TRIM_LEFT || TRIM_RIGHT){
					TrimRead.trim(r, TRIM_LEFT, TRIM_RIGHT, TRIM_QUAL, TRIM_MIN_LENGTH);
					TrimRead.trim(r2, TRIM_LEFT, TRIM_RIGHT, TRIM_QUAL, TRIM_MIN_LENGTH);
				}
				
				if(r2==null){
					final byte[] basesP=r.bases;
					final byte[] basesM=AminoAcid.reverseComplementBases(basesP, colorspace);
					basesUsed+=(basesM==null ? 0 : basesM.length);
					processRead(r, basesM);
					capSiteList(r, MAX_SITESCORES_TO_PRINT, PRINT_SECONDARY_ALIGNMENTS);
					assert(Read.CHECKSITES(r, basesM));
				}else{
					final byte[] basesP1=r.bases;
					final byte[] basesM1=AminoAcid.reverseComplementBases(basesP1, colorspace);
					final byte[] basesP2=r2.bases;
					final byte[] basesM2=AminoAcid.reverseComplementBases(basesP2, colorspace);
					basesUsed+=(basesM1==null ? 0 : basesM1.length);
					basesUsed+=(basesM2==null ? 0 : basesM2.length);
					assert(r2.bases==null || r2.bases.length<maxReadLength()) : 
						"Read "+r2.numericID+", length "+r2.bases.length+", exceeds the limit of "+maxReadLength()+"\n"+
						"You can map the reads in chunks by reformatting to fasta, then mapping with the setting 'fastareadlen="+maxReadLength()+"'";
					if(RCOMP_MATE){r2.reverseComplement();}
					processReadPair(r, basesM1, basesM2);
					capSiteList(r, MAX_SITESCORES_TO_PRINT, PRINT_SECONDARY_ALIGNMENTS);
					capSiteList(r2, MAX_SITESCORES_TO_PRINT, PRINT_SECONDARY_ALIGNMENTS);
					assert(Read.CHECKSITES(r, basesM1));
					assert(Read.CHECKSITES(r2, basesM2));
				}

				if(UNTRIM && (TRIM_LEFT || TRIM_RIGHT)){
					TrimRead.untrim(r);
					TrimRead.untrim(r2);
				}

				if(MAKE_MATCH_HISTOGRAM){readstats.addToMatchHistogram(r);}
				if(MAKE_INSERT_HISTOGRAM && r.paired()){readstats.addToInsertHistogram(r, (SAME_STRAND_PAIRS || !REQUIRE_CORRECT_STRANDS_PAIRS));}

				if(translateToBaseSpace && r.numSites()>0){
					SiteScore ss=r.topSite();
					r.start=ss.start;
					r.stop=ss.stop;
					r.chrom=ss.chrom;
					r.setStrand(ss.strand);
					r.mapScore=ss.score;
					Read rt=tcr.translateToBasespace(r);
					if(rt!=null){
						readlist.set(i, rt);
					}
				}
			}
			
//			System.err.println("Returning a list..."+"\n"+readlist);
			
			writeList(readlist, black, ln.id);
			
			
					//System.err.println("Left from adding list "+readlist.get(0).numericID);
			
			cris.returnList(ln, readlist.isEmpty());
//			if(count>0){
//				cris.returnList(ln, readlist.isEmpty());
//				count--;
//			}
			//System.err.println("Waiting on a list...");
			ln=cris.nextList();
			readlist=ln.list;
		}
		
		
		
		//System.err.println("Returning a list... (final)");
		assert(readlist.isEmpty());
		cris.returnList(ln, readlist.isEmpty());
		finish();
	}
	
	private final void writeList(ArrayList<Read> readlist, boolean black, long listNumID){
		if(outStreamMapped!=null){
			ArrayList<Read> x=new ArrayList<Read>(readlist.size());
			for(Read r1 : readlist){
				if(r1!=null){
					Read r2=r1.mate;
					if(r1.mapped() || (r2!=null && r2.mapped())){
						if(!black || !Blacklist.inBlacklist(r1)){x.add(r1);}
					}
				}
			}
			outStreamMapped.add(x, listNumID);
		}
		
		if(outStreamBlack!=null){
			ArrayList<Read> x=new ArrayList<Read>(readlist.size());
			for(Read r1 : readlist){
				if(black && Blacklist.inBlacklist(r1)){x.add(r1);}
			}
			outStreamBlack.add(x, listNumID);
		}
		
		if(BBSplitter.streamTable!=null || BBSplitter.TRACK_SET_STATS || BBSplitter.TRACK_SCAF_STATS){
			BBSplitter.printReads(readlist, listNumID, null, CLEARZONE1());
		}
		
		if(outStreamUnmapped!=null){
			ArrayList<Read> x=new ArrayList<Read>(readlist.size());
			for(Read r1 : readlist){
				if(r1!=null){
					Read r2=r1.mate;
					if(!(r1.mapped() || (r2!=null && r2.mapped()))){
						x.add(r1);
					}
				}
			}
			outStreamUnmapped.add(x, listNumID);
		}
		
//		System.err.println("outputStream = "+outputStream==null ? "null" : "real");
		if(outStream!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
			if(DONT_OUTPUT_UNMAPPED_READS){removeUnmapped(readlist);}
			if(black && DONT_OUTPUT_BLACKLISTED_READS){removeBlacklisted(readlist);}
			for(Read r : readlist){
				if(r!=null){
					r.obj=null;
					assert(r.bases!=null);
					if(r.sites!=null && r.sites.isEmpty()){r.sites=null;}
				}
			}
//			System.err.println("Adding list of length "+readlist.size());
			outStream.add(readlist, listNumID);
		}
	}
	
	/** Returns max possible quick score for this read, or -1 if it cannot be mapped for quality reasons. 
	 * A positive score will be returned if it CAN be mapped, but no hits are found. */
	public final int quickMap(final Read r, final byte[] basesM){
		final AbstractIndex index=index();
		byte[] basesP=r.bases;
		basesAtQuickmap+=basesP.length;
		if(basesP.length<KEYLEN){return 0;}
		assert(basesP.length>=KEYLEN);
		
		if(PERFECTMODE || SEMIPERFECTMODE){//Imperfect reads cannot map perfectly.
			if(r.containsUndefined()){return-1;} 
		}else if(DISCARD_MOSTLY_UNDEFINED_READS){
			int n=r.countUndefined();
			if(n>25 && basesP.length-n<n){return -1;}
		}
		
		final int keyProbLen=basesP.length-KEYLEN+1;
		final float[] keyProbs=index.keyProbArray();
		int[] offsets;
		
		float keyDen2=((maxDesiredKeys*KEYLEN)/(float)basesP.length);
		keyDen2=Tools.max(minKeyDensity, keyDen2);
		keyDen2=Tools.min(keyDensity, keyDen2);
		
		float keyDen3;
		if(basesP.length<=50){
			keyDen3=maxKeyDensity;
		}else if(basesP.length>=200){
			keyDen3=maxKeyDensity-0.5f;
		}else{
			keyDen3=maxKeyDensity-0.003333333333f*(basesP.length-50); //0.003333... = 0.5/150
		}
		
		keyDen3=Tools.max(keyDensity, keyDen3);
		
		if(GENERATE_KEY_SCORES_FROM_QUALITY){
			QualityTools.makeKeyProbs(r.quality, KEYLEN, keyProbs);
			
			boolean offsetsMode3=true;
			if(offsetsMode3){
				offsets=KeyRing.makeOffsets3(keyProbs, r.bases.length, KEYLEN, keyDen2, keyDen3, 2, (PERFECTMODE || SEMIPERFECTMODE));
			}else{
				//Old version; usually worse.
				offsets=KeyRing.makeOffsets2(keyProbs, r.bases.length, KEYLEN, keyDen2, keyDen3, 2);
				int numKeys=(offsets==null ? 0 : offsets.length+1);
				int maxRounds=0;//(PERFECTMODE || SEMIPERFECTMODE) ? 0 : 9999;//(numKeys)/2;
				while(maxRounds>0 && offsets!=null && offsets.length<numKeys){
					numKeys=offsets.length;
					offsets=QualityTools.modifyOffsets(offsets, keyProbs);
					maxRounds--;
				}
			}
		}else{
			offsets=KeyRing.makeOffsets(r.quality, KEYLEN, keyDensity, 2);
		}
		if(verbose){System.err.println("Made offsets: "+Arrays.toString(offsets));}

		if(offsets==null || offsets.length<AbstractIndex.MIN_APPROX_HITS_TO_KEEP || (r.quality!=null && r.avgQuality()<2)){return -1;}
		
		
		final byte[] baseScoresP=index.getBaseScoreArray(basesP.length, 0);
		final int[] keyScoresP=index.getKeyScoreArray(offsets.length, 0);
		
		if(AbstractIndex.USE_EXTENDED_SCORE){
			if(AbstractIndex.GENERATE_BASE_SCORES_FROM_QUALITY){
				QualityTools.makeByteScoreArray(r.quality, 100, baseScoresP, true);
			}
		}
		
		if(GENERATE_KEY_SCORES_FROM_QUALITY){
			int a=BASE_KEY_HIT_SCORE;
			int baseKeyScore=a/8;
			int range=a-baseKeyScore;
			final int[] keyScoresAll=new int[keyProbLen];
			QualityTools.makeKeyScores(keyProbs, keyProbLen, range, baseKeyScore, keyScoresAll);
			
			float probAllErrors=1f;
			for(int i=0; i<offsets.length; i++){
				keyScoresP[i]=keyScoresAll[offsets[i]];
				probAllErrors*=keyProbs[offsets[i]];
			}
			if(probAllErrors>0.50f){return -1;} //Default .5f; higher gives more false positives, lower gives more false negatives
			if(verbose){System.err.println("Prob all errors = "+probAllErrors+"\n\n");}
		}else{
			Arrays.fill(keyScoresP, BASE_KEY_HIT_SCORE);
		}
		if(verbose){System.err.println("Made key scores: "+Arrays.toString(keyScoresP));}
		
		keysUsed+=offsets.length;
		int maxScore=index.maxScore(offsets, baseScoresP, keyScoresP, basesP.length, true);
		if(verbose){System.err.println("Max Score: "+maxScore);}
		assert(maxScore>0);
		
		ArrayList<SiteScore> list=index.findAdvanced(basesP, basesM, r.quality, baseScoresP, keyScoresP, offsets, r.numericID);
		if(verbose){System.err.println("list: "+list);}
		
		r.sites=list;
		removeOutOfBounds(r, DONT_OUTPUT_UNMAPPED_READS, (outStream!=null && outStream.SAM), EXPECTED_LEN_LIMIT);
		assert(Read.CHECKSITES(list, r.bases, basesM, r.numericID));
		if(FORBID_SELF_MAPPING){forbidSelfMapping(list, r.originalSite);}
		
		if(list==null || list.isEmpty()){
			r.sites=null;
		}else{
			r.sites=list;
			if(!SLOW_ALIGN && AbstractIndex.USE_AFFINE_SCORE){
				for(SiteScore ss : list){ss.slowScore=ss.quickScore;}
			}
		}
//		assert(r.list!=null); //Less efficient, but easier to code later.
		
		return maxScore;
	}
	
	
	/** 
	 * Returns number of scores of at least maxImperfectSwScore.
	 * If problems are encountered such that it is prudent to do slow-alignment, a number lower than 1 will be returned. 
	 */
	final int scoreNoIndels(final Read r, final byte[] basesP, final byte[] basesM, final int maxSwScore, final int maxImperfectSwScore){
		
		if(!SLOW_ALIGN || r.numSites()==0){return 0;}
		
		int numPerfectScores=0;
		int numNearPerfectScores=0;
		int bestScoreNoIndel=Integer.MIN_VALUE;
		
		boolean forceSlow=false;

		for(int j=0; j<r.sites.size(); j++){

			SiteScore ss=r.sites.get(j);
			int oldScore=ss.score;
			int sslen=ss.stop-ss.start+1;
//			assert(false) : ss+", "+ss.quickScore+", "+ss.score+", "+ss.slowScore+", "+ss.pairedScore;
			
			final byte[] bases=(ss.strand==Gene.PLUS ? basesP : basesM);
			
			if(AbstractIndex.USE_AFFINE_SCORE && ss.quickScore==maxSwScore){
				assert(ss.stop==ss.start+r.bases.length-1) : ss.toText()+", "+maxSwScore+", "+maxImperfectSwScore+
					", "+r.bases.length+", "+(ss.start+r.bases.length-1);
			}
			
			if(verbose){System.err.print("C) Changing SiteScore from "+ss+"\n");}
			
			int slowScoreNoIndel;
			if(ss.perfect){
				if(verbose){System.err.print("C1");}
				numNearPerfectScores++;
				assert(ss.semiperfect);
				assert(ss.stop-ss.start+1==bases.length);
//				assert(maxSwScore==msa.scoreNoIndels((ss.strand==Gene.PLUS ? basesP : basesM), ss)); //TODO Disable this very slow assertion
				slowScoreNoIndel=maxSwScore;
				ss.slowScore=slowScoreNoIndel;
				ss.score=slowScoreNoIndel;
				ss.gaps=null;
			}else{
				if(verbose){System.err.print("C2");}
				ChromosomeArray cha=Data.getChromosome(ss.chrom);
				slowScoreNoIndel=msa.scoreNoIndels(bases, cha.array, ss.start, (sslen==bases.length ? ss : null));
				
				//This block is to correct situations where slow align does not get called, 
				//so one near-perfect alignment is found and one missed, because the read should align to stop, not start.
				if(slowScoreNoIndel<oldScore && oldScore>=maxImperfectSwScore && ss.stop-ss.start+1!=bases.length){
					int slowScoreNoIndel2=msa.scoreNoIndels(bases, cha.array, ss.stop-bases.length+1, null);
					if(slowScoreNoIndel2>=maxImperfectSwScore){
						slowScoreNoIndel=slowScoreNoIndel2;
						ss.start=ss.stop-bases.length+1;
						ss.setPerfect(bases);
					}
				}
				
				ss.slowScore=slowScoreNoIndel;
				ss.score=slowScoreNoIndel;
				
				if(slowScoreNoIndel>=maxImperfectSwScore){
					if(verbose){System.err.print("C3");}
					numNearPerfectScores++;
					
					ss.stop=ss.start+bases.length-1;
					ss.gaps=null;
					if(slowScoreNoIndel>=maxSwScore){
						if(verbose){System.err.print("C4");}
						assert(slowScoreNoIndel==maxSwScore) : slowScoreNoIndel+">"+maxSwScore;
						numPerfectScores++;
						ss.perfect=ss.semiperfect=true;
					}else{
						if(verbose){System.err.print("C5");}
						assert(!ss.perfect);
						ss.setPerfect(bases);
						assert(!ss.perfect);
					}
					if(QUICK_MATCH_STRINGS && !ss.perfect && (PRINT_SECONDARY_ALIGNMENTS || slowScoreNoIndel>=bestScoreNoIndel)){
						ss.match=msa.genMatchNoIndels(bases, cha.array, ss.start);
					}
				}else if(oldScore>=maxImperfectSwScore){
					if(verbose){System.err.print("C6");}
					forceSlow=true;
				}
			}

			if(verbose){System.err.print("\nto "+ss+"\n");}

			bestScoreNoIndel=Tools.max(ss.slowScore, bestScoreNoIndel);
//			assert(CHECKSITE(ss, bases));
		}
		return (forceSlow ? -numNearPerfectScores : numNearPerfectScores);
	}
	
//	@Deprecated
//	/** Assumes list is sorted */
//	public final void genMatchString_old(final Read r, final byte[] basesP, final byte[] basesM, final int maxImperfectSwScore, final int maxSwScore, boolean setSSScore, final boolean recur){
//		assert(Read.CHECKSITES(r, basesM));
//		assert(checkTopSite(r));
//		assert(r.mate!=null || r.list==null || r.list.size()==0 || r.list.get(0).score==r.mapScore) : "\n"+r.toText(false)+"\n"; //Came from BBMapAcc; not sure if it is correct
//		assert(msa!=null);
//		if(r.list==null || r.list.isEmpty()){
//			r.chrom=-1;
//			assert(r.mate!=null || r.list==null || r.list.size()==0 || r.list.get(0).score==r.mapScore) : "\n"+r.toText(false)+"\n";
//			return;
//		}
//		
//		final SiteScore ss=r.list.get(0);
//		assert(r.start==ss.start);
//		assert(r.stop==ss.stop);
//		assert(r.chrom==ss.chrom);
//		assert(r.strand()==ss.strand);
//		
//		final int minMsaLimit;
//		if(PAIRED){
////			minMsaLimit=-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE*maxSwScore);
//			minMsaLimit=-1+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore);
////			minMsaLimit=0;
//		}else{
//			minMsaLimit=-1+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO*maxSwScore);
////			minMsaLimit=0;
//		}
//		
//		if(GEN_MATCH_FAST){
////			r.start=ss.start;
////			r.stop=ss.stop;
////			r.chrom=ss.chrom;
////			r.strand=ss.strand;
//			
//			assert(!(SLOW_ALIGN || AbstractIndex.USE_EXTENDED_SCORE) || AbstractIndex.GENERATE_BASE_SCORES_FROM_QUALITY || 
//					(ss.slowScore==maxSwScore) == r.perfect()) : 
//				r.bases.length+", "+ss.toText()+", "+maxSwScore+", "+ss.slowScore+", "+r.perfect();
//			
//			//TODO: This WAS disabled because I saw a read marked perfect with a sub in it, probably with quality 0 at that point.
//			if((SLOW_ALIGN || AbstractIndex.USE_EXTENDED_SCORE) && r.perfect()){
//				assert(r.stop-r.start==(r.bases.length-1));
//				r.match=ss.match=makePerfectMatchString(r.bases.length);
//				assert(ss.isPerfect(ss.plus() ? basesP : basesM)) : r; //TODO: Slow assertion
//				assert(Read.CHECKSITES(r, basesM));
//				assert(checkTopSite(r)); // TODO remove this
//			}else
//			{
//				int oldScore=ss.slowScore;
//				assert(r.start==ss.start && r.stop==ss.stop);
//				assert(r.gaps==null || r.gaps[0]==r.start && r.gaps[r.gaps.length-1]==r.stop);
//				int padding=(ss.perfect || ss.semiperfect ? 0 : Tools.max(SLOW_ALIGN_PADDING, 6));
//				
//				if(verbose){System.err.println("Attempting to realign read:\n"+r+"\npadding="+padding+"\nrescued="+r.rescued());}
//				
//				TranslateColorspaceRead.realign_new(r, msa, padding, true, minMsaLimit, MAX_INDEL<1); //Also generates the match string
//				r.gaps=ss.gaps=GapTools.fixGaps(r.start, r.stop, r.gaps, Shared.MINGAP);
//				
//				if(verbose){System.err.println("Realigned read:\n"+r+"\npadding="+padding+"\nrescued="+r.rescued());}
//				assert(Read.CHECKSITES(r, basesM)); //***123
//				assert(ss==r.list.get(0)) : "Site order changed";
//				
//				if(r.mapScore<ss.slowScore){
//					if(verbose){System.err.println("---- A ----");}
//					if(verbose){
//						System.err.print("Read "+r.numericID+": "+r.start+","+r.stop+": "+r.mapScore+"<"+ss.slowScore);
//					}
//					
//					int extra=(MAX_INDEL>0 ? 80 : 20)+SLOW_ALIGN_PADDING;
//					int expectedLen=GapTools.calcGrefLen(r.start, r.stop, r.gaps); //TODO Gaps should be correct here!!!
//					int remaining=(msa.maxColumns-expectedLen-2);
//					extra=Tools.max(0, Tools.min(remaining/2, extra));
//					TranslateColorspaceRead.realign_new(r, msa, extra, true, minMsaLimit, false);
//					r.gaps=ss.gaps=GapTools.fixGaps(r.start, r.stop, r.gaps, Shared.MINGAP);
//					assert(Read.CHECKSITES(r, basesM));
//					
//					if(verbose){
//						System.err.println(" -> "+r.start+","+r.stop+","+r.mapScore+
//								(r.originalSite==null ? "" : "\t*"+r.originalSite)+"\t(extra = "+extra+")");
//					}
//				}
//				if(verbose){System.err.println("---- B ----");}
//				assert(ss==r.list.get(0)) : "Site order changed";
//				ss.match=r.match;
//				
//				//TODO: This is new, make sure it does not break anything (Note: It did, but should be fixed now)
//				assert(Read.CHECKSITES(r, basesM));
//				{
//					ss.slowScore=r.mapScore;
//					if(setSSScore){ss.score=r.mapScore;}
//					assert(r.mate!=null || r.list==null || r.list.size()==0 || r.list.get(0).score==r.mapScore) : "\n"+r.toText(false)+"\n";
//					if(ss.start!=r.start || ss.stop!=r.stop){
//						if(verbose){
//							System.err.println("---- C ----");
//							System.err.println(ss);
//							System.err.println(r.list.get(0));
//							System.err.println(r.start+","+r.stop+","+r.mapScore);
//						}
//						ss.start=r.start;
//						ss.stop=r.stop;
//						ss.match=r.match;
//						if(!AMBIGUOUS_RANDOM || !r.ambiguous()){
//							if(verbose){
//								System.err.println("---- D ----");
//								System.err.println(ss);
//								System.err.println(r.list.get(0));
//								System.err.println(r.start+","+r.stop+","+r.mapScore);
//							}
//							assert(ss==r.list.get(0)) : "Site order changed\n"+ss+"\n"+r.list.get(0)+"\n"; assert(checkTopSite(r)); // TODO remove this
//							if(!r.paired()){
//								Tools.mergeDuplicateSites(r.list, false, false);
//								Collections.sort(r.list);
//								final SiteScore ss2=r.list.get(0);
//								if(ss!=ss2){//Fixes a super rare case
//									ss.setPerfect(ss.plus() ? basesP : basesM, false);
////									System.err.println("**********************\n\nCalled setPerfect on "+ss+"\tp="+ss.perfect+", sp="+ss.semiperfect);
////									assert(Read.CHECKSITE(ss, ss.plus() ? basesP : basesM, r.numericID));
////									System.err.println("INDEX = "+r.list.indexOf(ss));
////									ss2.setPerfect(r.bases, false);
////									assert(Read.CHECKSITE(ss2, ss2.plus() ? basesP : basesM, r.numericID));
//									r.setFromSite(ss2);
////									System.err.println("**********************\n\nCalled setPerfect on "+ss+"\tp="+ss.perfect+", sp="+ss.semiperfect);
////									assert(Read.CHECKSITE(ss, ss.plus() ? basesP : basesM, r.numericID));
////									assert(Read.CHECKSITES(r, basesM));
//									genMatchString(r, basesP, basesM, maxImperfectSwScore, maxSwScore, setSSScore, recur);
////									r.setPerfectFlag(maxSwScore);
//									assert(checkTopSite(r));//124
//									assert(Read.CHECKSITES(r, basesM));//124
//									return;
//								}
//							}else{
//								for(int i=r.list.size()-1; i>0; i--){
//									if(ss.positionalMatch(r.list.get(i), true)){r.list.remove(i);}
//								}
//							}
//						}
//						assert(ss==r.list.get(0)) : "Site order changed\n"+ss+"\n"+r.list.get(0)+"\n";
//						assert(checkTopSite(r)); // TODO remove this
//					}
//					assert(ss==r.list.get(0)) : "Site order changed\n"+ss+"\n"+r.list.get(0)+"\n";
//					assert(checkTopSite(r)); // TODO remove this
//					if(verbose){
//						System.err.println("---- D2 ----");
//						System.err.println(ss);
//						System.err.println(r.list.get(0));
//						System.err.println(r.start+","+r.stop+","+r.mapScore);
//					}
//				}
//				assert(ss==r.list.get(0)) : "Site order changed";
//				assert(checkTopSite(r)); // TODO remove this
//				if(verbose){
//					System.err.println("---- D3 ----");
//					System.err.println(ss);
//					System.err.println(r.list.get(0));
//					System.err.println(r.start+","+r.stop+","+r.mapScore);
//					System.err.println("Checking perfect status: r.perfect="+r.perfect()+", ss.perfect="+ss.perfect+", ss.semi="+ss.semiperfect+
//							", maxSwScore="+maxSwScore+", r.score="+r.mapScore+", ss.slowScore="+ss.slowScore+"\n"+r);
//				}
//				r.setPerfectFlag(maxSwScore);
//				if(verbose){
//					System.err.println("---- E ----");
//					System.err.println("Checking perfect status: r.perfect="+r.perfect()+", ss.perfect="+ss.perfect+", ss.semi="+ss.semiperfect+
//							", maxSwScore="+maxSwScore+", r.score="+r.mapScore+", ss.slowScore="+ss.slowScore+"\n"+r);
//				}
//				assert(r.match==ss.match) : r.perfect()+", "+ss.perfect+", "+ss.semiperfect+"\n"+new String(r.match)+"\n"+new String(ss.match);
//				if(r.perfect()){
//					ss.perfect=ss.semiperfect=true;
//					assert(Read.CHECKSITES(r, basesM)) : r.perfect()+", "+ss.perfect+", "+ss.semiperfect+"\n"+new String(r.match)+"\n"+new String(ss.match); //***123
//				}else{
//					final byte[] bases=(ss.plus() ? basesP : basesM);
////					if(r.match!=null && r.containsNonNM()){
////						ss.perfect=ss.semiperfect=false; //This should be fine, but failed when a match string contained X.
//					if(r.match!=null && r.containsSDI()){
//						ss.perfect=ss.semiperfect=false;
////						ss.setPerfect(bases, false);
////						r.setPerfect(ss.perfect);
//						assert(Read.CHECKSITES(r, basesM)) : r.perfect()+", "+ss.perfect+", "+ss.semiperfect+"\n"+new String(r.match)+"\n"+new String(ss.match); //***123
//					}else{
//						//rare
//						ss.setPerfect(bases, false);
//						r.setPerfect(ss.perfect);
//						assert(Read.CHECKSITES(r, basesM)) : r.perfect()+", "+ss.perfect+", "+ss.semiperfect+"\n"+new String(r.match)+"\n"+new String(ss.match); //***123
//					}
//				}
////				assert(Read.CHECKSITES(r, basesM)) : r.perfect()+", "+ss.perfect+", "+ss.semiperfect; //***123
//				assert(checkTopSite(r)); // TODO remove this
//				assert(r.perfect()==ss.perfect);
//				assert(!r.perfect() || r.stop-r.start==(r.bases.length-1));
//				if(verbose){
//					System.err.println("Checking perfect status: r.perfect="+r.perfect()+", ss.perfect="+ss.perfect+", ss.semi="+ss.semiperfect+
//							", maxSwScore="+maxSwScore+", r.score="+r.mapScore+", ss.slowScore="+ss.slowScore+"\n"+r);
//				}
//			}
//		}else{
//			if(verbose){System.err.println("---- F ----");}
//			byte[] bases=(ss.plus() ? basesP : basesM);
//			//				int[] swscoreArray=msa.fillAndScore(bases, ss, 0);
//			
//			if(r.perfect()){
//				r.match=makePerfectMatchString(r.bases.length);
//			}else{
//				ChromosomeArray cha=Data.getChromosome(ss.chrom);
//				assert(false) : "TODO: This does not take strand into account";
//				if(ss.slowScore>=maxImperfectSwScore){
//					//TODO
//				}
//
//				if(msa!=null){
//					assert(false) : "0 is not good here; try a non-indel match string.";
//					int[] max=msa.fillLimited(bases, cha.array, ss.start, ss.stop, 0, ss.gaps);
//					//					System.err.print("*");
//					r.match=msa.traceback(bases, cha.array, ss.start, ss.stop, max[0], max[1], max[2], ss.gaps!=null);
//				}
//			}
//		}
//		if(verbose){System.err.println("---- G ----");}
//
//		assert(Read.CHECKSITES(r, basesM)); //***123
//		assert(checkTopSite(r)); // TODO remove this
//		if((!AMBIGUOUS_RANDOM || !r.ambiguous()) && recur && r.list.get(0)!=ss){
//			r.setFromTopSite();
//			genMatchString(r, basesP, basesM, maxImperfectSwScore, maxSwScore, setSSScore, recur);
//			assert(checkTopSite(r)); // TODO remove this
//		}else{
//			
//			//Corrects a mysterious bug encountered with paired reads, in which semiperfect reads are not flagged semiperfect.
//			//TODO: Find out reason for this and correct it, then disable this block.
//			if(verbose){
//				System.err.println("Checking perfect status: r.perfect="+r.perfect()+", ss.perfect="+ss.perfect+", ss.semi="+ss.semiperfect+
//						", maxSwScore="+maxSwScore+", r.score="+r.mapScore+", ss.slowScore="+ss.slowScore+"\n"+r);
//			}
//			assert(Read.CHECKSITES(r, basesM));//***123
//			assert(checkTopSite(r)); // TODO remove this
//			if(!r.perfect()){
//				if(verbose){System.err.println("Correcting perfect status");}
//				if(r.mate!=null && r.list!=null && r.list.size()>0){
//					SiteScore ss2=r.list.get(0);
//					if(verbose){System.err.println("Checking perfect status2: ss2.perfect="+ss2.perfect+", ss2.semi="+ss2.semiperfect+"\nss="+ss+"\nss2="+ss2);}
//					byte[] bases=(ss2.plus() ? basesP : basesM);
//					ss2.setPerfect(bases, false);
//					r.setPerfect(ss2.perfect);
//					if(verbose){System.err.println("New perfect status: r.perfect="+r.perfect()+", ss2.perfect="+ss2.perfect+", ss2.semi="+ss2.semiperfect);}
//					assert(Read.CHECKSITE(ss2, bases, r.numericID));
//					assert(checkTopSite(r)); // TODO remove this
//				}
//			}
//			assert(Read.CHECKSITES(r, basesM));
//			assert(checkTopSite(r)); // TODO remove this
//		}
//		assert(checkTopSite(r)); // TODO remove this
//	}
	
	/** Assumes list is sorted */
	public final void genMatchString(final Read r, final byte[] basesP, final byte[] basesM, final int maxImperfectSwScore, final int maxSwScore, boolean setSSScore, final boolean recur){
		if(verbose){System.err.println("\n\n\n\n\ngenMatchString for read\n"+r+"\n\n\n\n\n");}
		assert(Read.CHECKSITES(r, basesM));
		assert(checkTopSite(r));
		
		assert(r.mate!=null || r.numSites()==0 || r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n"; //Came from BBMapAcc; not sure if it is correct
		assert(msa!=null);
		if(r.numSites()==0){
			r.chrom=-1;
			assert(r.mate!=null || r.numSites()==0 || r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";
			return;
		}
		
		if(PRINT_SECONDARY_ALIGNMENTS){
			capSiteList(r, MAX_SITESCORES_TO_PRINT+3, PRINT_SECONDARY_ALIGNMENTS);
		}
		
		if(QUICK_MATCH_STRINGS && PRINT_SECONDARY_ALIGNMENTS && USE_SS_MATCH_FOR_PRIMARY){} //TODO What was this line for?
		
		int best=Integer.MIN_VALUE;
		int scoreChanged=0;
		
		for(int i=0; i<r.sites.size(); i++){
			SiteScore ss=r.sites.get(i);
			
			if(verbose){System.err.println("**************** best="+best+", scoreChanged="+scoreChanged+"\nconsidering ss "+ss);}
			
			if(i>0){
				if(best>=ss.slowScore && !PRINT_SECONDARY_ALIGNMENTS){
					if(verbose){System.err.println("break triggered by low score");}
					break;
				}
			}
			
			int oldScore=ss.slowScore;
			if(ss.match==null || (i==0 && !USE_SS_MATCH_FOR_PRIMARY)){
				genMatchStringForSite(r.numericID, ss, basesP, basesM, maxImperfectSwScore, maxSwScore, r.mate);
				if(setSSScore){ss.score=ss.slowScore;}
			}
			if(i>0 && ss.match==null && !r.paired()){r.sites.remove(i);}
			else{
				if(oldScore!=ss.slowScore){scoreChanged++;}
				best=Tools.max(ss.slowScore, best);
			}

			if(verbose){System.err.println("**************** best="+best+", scoreChanged="+scoreChanged+"\nconsidered ss "+ss);}
		}
		
		if(verbose){System.err.println("Finished basic match generation. best="+best+", scoreChanged="+scoreChanged+", AMBIGUOUS_RANDOM="+AMBIGUOUS_RANDOM+", ambiguous="+r.ambiguous());}
		if(scoreChanged>0 && (!AMBIGUOUS_RANDOM || !r.ambiguous())){
			if(!r.paired()){
				if(verbose){System.err.println("GMS 1");}
				Tools.mergeDuplicateSites(r.sites, false, false);
				Collections.sort(r.sites);
				int prevScore=0;
				for(int i=0; i<r.sites.size() && i<MAX_SITESCORES_TO_PRINT; i++){
					if(verbose){System.err.println("GMS 2");}
					SiteScore ss=r.sites.get(i);
					if(ss.match==null){
						if(verbose){System.err.println("GMS 3");}
						genMatchStringForSite(r.numericID, ss, basesP, basesM, maxImperfectSwScore, maxSwScore, r.mate);
						if(setSSScore){ss.score=ss.slowScore;}
						if(i>0 && ss.match==null){r.sites.remove(i);}
						i--;
					}
					if(i>0 || !PRINT_SECONDARY_ALIGNMENTS){
						if(verbose){System.err.println("GMS 4");}
						break;
					}
				}
			}else{
				if(verbose){System.err.println("GMS 5");}
				SiteScore ss=r.topSite();
				for(int i=r.sites.size()-1; i>0; i--){
					if(verbose){System.err.println("GMS 6");}
					if(ss.positionalMatch(r.sites.get(i), true)){r.sites.remove(i);}
				}
			}
		}
		
		
		final SiteScore ss=r.topSite();
		assert(ss==r.topSite());
		
//		assert(ss.slowScore>0) : ss.slowScore+", "+best+", "+r.mapScore;
		
		r.start=ss.start;
		r.stop=ss.stop;
		r.chrom=ss.chrom;
		r.setStrand(ss.strand);
		r.match=ss.match;
		r.gaps=ss.gaps;
		r.mapScore=ss.slowScore;
		r.setPerfect(ss.perfect());
		r.setRescued(ss.rescued());
		
		assert(Read.CHECKSITES(r, basesM));
		assert(checkTopSite(r));
		
//		assert(false) : r.numericID+", "+ss.slowScore+", "+r.mapScore;
	}
	
	
	protected final int genMatchStringForSite(final long id, final SiteScore ss, final byte[] basesP, final byte[] basesM, final int maxImperfectSwScore, final int maxSwScore, final Read mate){
		final byte[] bases=ss.plus() ? basesP : basesM;
		assert(Read.CHECKSITE(ss, bases, id));
		assert(msa!=null);
		
		
		final int minMsaLimit;
		if(PAIRED){
//			minMsaLimit=-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE*maxSwScore);
			minMsaLimit=-1+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore);
//			minMsaLimit=0;
		}else{
			minMsaLimit=-1+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO*maxSwScore);
//			minMsaLimit=0;
		}
		
		if(GEN_MATCH_FAST){
			
			assert(!(SLOW_ALIGN || AbstractIndex.USE_EXTENDED_SCORE) || AbstractIndex.GENERATE_BASE_SCORES_FROM_QUALITY || 
					(ss.slowScore==maxSwScore) == ss.perfect()) : 
				bases.length+", "+ss.toText()+", "+maxSwScore+", "+ss.slowScore+", "+ss.perfect()+", "+ss.semiperfect();
			
			//TODO: This WAS disabled because I saw a read marked perfect with a sub in it, probably with quality 0 at that point.
			if((SLOW_ALIGN || AbstractIndex.USE_EXTENDED_SCORE) && ss.perfect()){
				assert(ss.stop-ss.start==(bases.length-1));
				ss.match=makePerfectMatchString(bases.length);
				assert(ss.isPerfect(bases)) : id+", "+ss; //TODO: Slow assertion
			}else{
				int oldScore=ss.slowScore;
				assert(ss.gaps==null || ss.gaps[0]==ss.start && ss.gaps[ss.gaps.length-1]==ss.stop);
				int padding=(ss.perfect || ss.semiperfect ? 0 : Tools.max(SLOW_ALIGN_PADDING, 6));
				
				if(verbose){System.err.println("Attempting to realign read:\n"+id+", "+ss+"\npadding="+padding+"\nrescued="+ss.rescued());}
				
				TranslateColorspaceRead.realign_new(ss, bases, msa, padding, true, minMsaLimit, MAX_INDEL<1, id); //Also generates the match string
				ss.gaps=GapTools.fixGaps(ss.start, ss.stop, ss.gaps, Shared.MINGAP);
				
				if(verbose){System.err.println("Realigned read:\n"+id+", "+ss+"\npadding="+padding+"\nrescued="+ss.rescued());}
				assert(Read.CHECKSITE(ss, bases, id));
				
				if(ss.slowScore<oldScore){
					if(verbose){System.err.println("---- A ----");}
					if(verbose){
						System.err.print("Read "+id+": "+ss.start+","+ss.stop+": "+oldScore+">"+ss.slowScore);
					}
					
					int extra=(MAX_INDEL>0 ? 80 : 20)+SLOW_ALIGN_PADDING;
					int expectedLen=GapTools.calcGrefLen(ss.start, ss.stop, ss.gaps); //TODO Gaps should be correct here!!!
					int remaining=(msa.maxColumns-expectedLen-2);
					extra=Tools.max(0, Tools.min(remaining/2, extra));
					TranslateColorspaceRead.realign_new(ss, bases, msa, extra, true, minMsaLimit, false, id);
					ss.gaps=GapTools.fixGaps(ss.start, ss.stop, ss.gaps, Shared.MINGAP);
					assert(Read.CHECKSITE(ss, bases, id));
					
					if(verbose){
						System.err.println("\n-> "+ss.start+","+ss.stop+","+ss.slowScore+
								/*(r.originalSite==null ? "" : "\t*"+r.originalSite)+*/"\t(extra = "+extra+")");
					}
				}
				if(verbose){System.err.println("---- B ----");}
				assert(Read.CHECKSITE(ss, bases, id));
				
				if(verbose){
					System.err.println("---- D3 ----");
					System.err.println(ss);
					System.err.println("Checking perfect status: ss.perfect="+ss.perfect()+", ss.semi="+ss.semiperfect()+
							", maxSwScore="+maxSwScore+", ss.slowScore="+ss.slowScore);
				}
				ss.setPerfectFlag(maxSwScore, bases);
				if(verbose){
					System.err.println("---- E ----");
					System.err.println("Checking perfect status: ss.perfect="+ss.perfect()+", ss.semi="+ss.semiperfect()+
							", maxSwScore="+maxSwScore+", ss.slowScore="+ss.slowScore);
				}
				
				assert(Read.CHECKSITE(ss, bases, id));
			}
		}else{
			if(verbose){System.err.println("---- F ----");}
			ChromosomeArray cha=Data.getChromosome(ss.chrom);
			
			if(ss.perfect()){
				ss.match=makePerfectMatchString(bases.length);
			}else{
				assert(false) : "TODO: This does not take strand into account";
				if(ss.slowScore>=maxImperfectSwScore){
					//TODO
				}

				if(msa!=null){
					assert(false) : "0 is not good here; try a non-indel match string.";
					int[] max=msa.fillLimited(bases, cha.array, ss.start, ss.stop, 0, ss.gaps);
					//					System.err.print("*");
					ss.match=msa.traceback(bases, cha.array, ss.start, ss.stop, max[0], max[1], max[2], ss.gaps!=null);
				}
			}
		}
		if(verbose){System.err.println("---- G ----");}

		assert(Read.CHECKSITE(ss, bases, id));
		return ss.slowScore;
	}

	
	
	/** Returns the number of additional bases away that should be searched for slow align. 
	 * This should probably be called between quickMap and slowAlign, only on 
	 * sites where stop-start<=bases.length-1 */
	final void findTipDeletions(final Read r, final byte[] basesP, final byte[] basesM, final int maxSwScore, final int maxImperfectScore){

		boolean findRight=r.quality==null || (r.minQualityLastNBases(TIP_DELETION_MAX_TIPLEN)>=TIP_DELETION_MIN_QUALITY &&
				r.avgQualityLastNBases(TIP_DELETION_MAX_TIPLEN)>=TIP_DELETION_AVG_QUALITY);
		boolean findLeft=r.quality==null || (r.minQualityFirstNBases(TIP_DELETION_MAX_TIPLEN)>=TIP_DELETION_MIN_QUALITY &&
				r.avgQualityFirstNBases(TIP_DELETION_MAX_TIPLEN)>=TIP_DELETION_AVG_QUALITY);
		if(!findRight && !findLeft){
//			System.err.print(".");
			return;
		}
//		System.err.print("*");
		
		for(SiteScore ss : r.sites){
			final byte[] bases=(ss.strand==Gene.PLUS ? basesP : basesM);
			if(!ss.semiperfect && ss.slowScore<maxImperfectScore){
				boolean changed=findTipDeletions(ss, bases, maxImperfectScore, findRight, findLeft);
				if(changed){
					ss.match=null;
					ss.slowScore=msa.scoreNoIndels(bases, ss.chrom, ss.start);
					assert(!ss.perfect);
					if(ss.slowScore==maxSwScore){
						ss.stop=ss.start+bases.length-1;
						ss.perfect=ss.semiperfect=true;
					}else{
						ss.perfect=false;
						ss.setPerfect(bases, true);
					}
				}
			}
		}
	}

	final boolean findTipDeletions(SiteScore ss, final byte[] bases, final int maxImperfectScore, boolean lookRight, boolean lookLeft){
		if(ss.slowScore>=maxImperfectScore /*&& ss.stop-ss.start<=basesP.length-1*/){return false;}
		assert(lookRight || lookLeft);
		assert(TIP_DELETION_MAX_TIPLEN>2);
		if(bases.length<=2*TIP_DELETION_MAX_TIPLEN){return false;}
		assert(TIP_DELETION_MAX_TIPLEN<bases.length);
		assert(TIP_DELETION_SEARCH_RANGE>0);
		
		int maxSearch=TIP_DELETION_SEARCH_RANGE;
		maxSearch=Tools.min(maxSearch, ALIGN_COLUMNS-(SLOW_RESCUE_PADDING+8+Tools.max(bases.length, ss.stop-ss.start)));
		if(maxSearch<1){return false;}
		
		boolean changed=false;
		
		if(lookRight){
			int x=findTipDeletionsRight(bases, ss.chrom, ss.stop, maxSearch, TIP_DELETION_MAX_TIPLEN);
			if(x>0){
				assert(x+ss.stop-ss.start<ALIGN_COLUMNS);
				ss.stop+=x;
				changed=true;
				maxSearch=Tools.min(maxSearch, ALIGN_COLUMNS-(SLOW_RESCUE_PADDING+8+Tools.max(bases.length, ss.stop-ss.start)));
				if(maxSearch<1){return changed;}
			}
		}

		if(lookLeft){
			int y=findTipDeletionsLeft(bases, ss.chrom, ss.start, maxSearch, TIP_DELETION_MAX_TIPLEN);
			if(y>0){
				assert(y+ss.stop-ss.start<ALIGN_COLUMNS);
				ss.start-=y;
				changed=true;
			}
		}
		return changed;
	}
	
	
	final void rescue(Read anchor, Read loose, byte[] basesP, byte[] basesM, int searchDist){
		
		//Lists should be sorted at this point, and have a paired score if they are paired.
		
		if(anchor.sites==null || anchor.sites.isEmpty()){return;}
		if(loose.sites==null){
			loose.sites=new ArrayList<SiteScore>(anchor.sites.size());
		}
		
		final int maxLooseSwScore=msa.maxQuality(basesP.length);
		final int maxAnchorSwScore=msa.maxQuality(anchor.bases.length);
		final int maxImperfectScore=msa.maxImperfectScore(basesP.length);
		
		final int bestLooseScore=loose.sites.isEmpty() ? 0 : loose.topSite().slowScore;
		final int bestAnchorScore=anchor.topSite().slowScore;
		
		if(bestLooseScore==maxLooseSwScore && bestAnchorScore==maxAnchorSwScore
				&& anchor.topSite().pairedScore>0){return;}

		int rescueScoreLimit=(int)(0.95f*bestAnchorScore);
//		int retainScoreLimit=(int)(bestLooseScore>0 ? 0.58f*bestLooseScore : 0.58f*maxLooseSwScore);
		int retainScoreLimit=Tools.max((int)(0.68f*bestLooseScore), (int)(0.4f*maxLooseSwScore));
		int retainScoreLimit2=Tools.max((int)(0.95f*bestLooseScore), (int)(0.55f*maxLooseSwScore));
		final int maxMismatches=PERFECTMODE ? 0 : (bestLooseScore>maxImperfectScore) ? 5 : (int)(0.60f*basesP.length-1); //Higher number is more lenient
		assert(PERFECTMODE || maxMismatches>1 || loose.bases.length<16) : loose; //Added the <16 qualifier when a 4bp read failed this assertion
		
		final boolean findTipDeletions=FIND_TIP_DELETIONS && bestLooseScore<maxImperfectScore;
		
		//Data for finding tip deletions
		final boolean findRight=findTipDeletions && loose.quality==null || (loose.minQualityLastNBases(TIP_DELETION_MAX_TIPLEN)>=TIP_DELETION_MIN_QUALITY
				&& loose.avgQualityLastNBases(TIP_DELETION_MAX_TIPLEN)>=TIP_DELETION_AVG_QUALITY);
		final boolean findLeft=findTipDeletions && loose.quality==null || (loose.minQualityFirstNBases(TIP_DELETION_MAX_TIPLEN)>=TIP_DELETION_MIN_QUALITY
				&& loose.avgQualityFirstNBases(TIP_DELETION_MAX_TIPLEN)>=TIP_DELETION_AVG_QUALITY);
		
//		int searchIntoAnchor=Tools.max(20, Tools.min(anchor.bases.length, loose.bases.length));
		for(SiteScore ssa : anchor.sites){
			if(ssa.slowScore<rescueScoreLimit){break;}
			if(ssa.pairedScore==0 && !ssa.rescued){
//				int searchIntoAnchor=ssa.stop-ssa.start-1+(anchor.bases.length/2); //Allows rescue of fragments half the length of a read
				int searchIntoAnchor=ssa.stop-ssa.start-1+(anchor.bases.length*11/16); //Allows rescue of fragments 68% the length of a read
				int loc;
				int idealStart;
				byte[] bases;
				byte strand=(SAME_STRAND_PAIRS ? ssa.strand : (byte)(ssa.strand^1));
				boolean searchRight=(SAME_STRAND_PAIRS ? strand==Gene.PLUS : strand==Gene.MINUS);
				assert(strand==0 || strand==1);
				
				if(SAME_STRAND_PAIRS){
					if(ssa.strand==Gene.MINUS){
						bases=basesM;
						loc=ssa.start+searchIntoAnchor;
						idealStart=ssa.start-AVERAGE_PAIR_DIST;
					}else{
						bases=basesP;
						loc=ssa.stop-searchIntoAnchor;
						idealStart=ssa.stop+AVERAGE_PAIR_DIST;
					}
				}else{
					if(ssa.strand==Gene.PLUS){
						bases=basesM;//opposite strand
						loc=ssa.stop-searchIntoAnchor;
						idealStart=ssa.stop+AVERAGE_PAIR_DIST;
					}else{
						bases=basesP;//opposite strand
						loc=ssa.start+searchIntoAnchor;
						idealStart=ssa.start-AVERAGE_PAIR_DIST;
					}
				}
//				loc-=20; //Search for overlapping read ends
				SiteScore ss=quickRescue(bases, ssa.chrom, strand, loc, searchDist+searchIntoAnchor, searchRight, 
						idealStart, maxMismatches, POINTS_MATCH, POINTS_MATCH2);
				
				if(ss!=null && ss.isInBounds()){
					int mismatches=ss.slowScore;
					ss.slowScore=0;
					if(mismatches<=maxMismatches){
						slowRescue(bases, ss, maxLooseSwScore, maxImperfectScore, findRight, findLeft);
						if(ss.score>retainScoreLimit && ss.isInBounds()){
							if(ss.score>retainScoreLimit2){//Set them as paired to make them more resistant to being discarded
								ss.pairedScore=Tools.max(ss.pairedScore, ss.slowScore+ssa.slowScore/4);
								ssa.pairedScore=Tools.max(ssa.pairedScore, ssa.slowScore+ss.slowScore/4);
								assert(ss.pairedScore>0);
								assert(ssa.pairedScore>0);
							}
							loose.sites.add(ss);
						}
					}
				}
			}else{
				assert(ssa.pairedScore>0);
				assert(ssa.pairedScore>ssa.quickScore || ssa.pairedScore>ssa.slowScore) : ssa.toText();
			}
		}
	}
	
	
	final void slowRescue(final byte[] bases, SiteScore ss, final int maxScore, final int maxImperfectScore, 
			boolean findTipDeletionsRight, boolean findTipDeletionsLeft){
		
		int swscoreNoIndel=msa.scoreNoIndels(bases, ss.chrom, ss.start);
		final int oldStart=ss.start;
		
		if(swscoreNoIndel<maxImperfectScore && MAX_INDEL>0){
			ss.slowScore=swscoreNoIndel;
			if(findTipDeletionsRight || findTipDeletionsLeft){
				boolean changed=findTipDeletions(ss, bases, maxImperfectScore, findTipDeletionsRight, findTipDeletionsLeft);
				if(changed){
					ss.match=null;
					swscoreNoIndel=msa.scoreNoIndels(bases, ss.chrom, ss.start);
				}
			}
			
			final int minMsaLimit=-CLEARZONE1e+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxScore);				
			
			final int minscore=Tools.max(swscoreNoIndel, minMsaLimit);
			final int[] swscoreArray=msa.fillAndScoreLimited(bases, ss.chrom, ss.start, ss.stop, SLOW_RESCUE_PADDING, minscore, ss.gaps);
			
			if(swscoreArray!=null){
				ss.slowScore=ss.score=swscoreArray[0];
				ss.start=swscoreArray[1];
				ss.stop=swscoreArray[2];
				
				if(QUICK_MATCH_STRINGS && swscoreArray!=null && swscoreArray.length==6 && swscoreArray[0]>=minscore && (PRINT_SECONDARY_ALIGNMENTS || USE_SS_MATCH_FOR_PRIMARY)){
					assert(swscoreArray.length==6) : swscoreArray.length;
					assert(swscoreArray[0]>=minscore) : "\n"+Arrays.toString(swscoreArray)+"\n"+minscore;
					ss.match=msa.traceback(bases, Data.getChromosome(ss.chrom).array, ss.start-SLOW_RESCUE_PADDING, ss.stop+SLOW_RESCUE_PADDING, 
							swscoreArray[3], swscoreArray[4], swscoreArray[5], ss.gaps!=null);
					ss.fixXY(bases, true, msa);
				}else{ss.match=null;}
				
			}else{
				ss.slowScore=ss.score=swscoreNoIndel;
				ss.start=oldStart;
				ss.stop=ss.start+bases.length-1;
			}
		}else{
			ss.slowScore=ss.score=swscoreNoIndel;
			ss.stop=ss.start+bases.length-1;
		}
		ss.pairedScore=ss.score+1;
		assert(ss.slowScore<=maxScore);
		ss.perfect=(ss.slowScore==maxScore);
		if(ss.perfect){ss.semiperfect=true;}
		else{ss.setPerfect(bases);}
	}
	
	
	protected static final void capSiteList(Read r, int cap, boolean printSecondary){
		if(r==null || r.sites==null || cap<0){return;}
		if(cap==0){r.sites=null;}
		else{
			for(int i=r.sites.size()-1; i>=cap; i--){r.sites.remove(i);}
		}
		if(!printSecondary || r.numSites()<2){return;}
		int max=r.topSite().slowScore;
		int min=Tools.min(max-500, (int)(max*.95f));
		for(int i=r.sites.size()-1; i>0; i--){
			if(r.sites.get(i).slowScore<min){r.sites.remove(i);}
		}
//		assert(false) : r.mapScore+", "+max+", "+cap+", "+r.list;
//		assert(r.list.size()<2) : "\n"+max+", "+min+", "+r.list+"\n";
	}
	
	protected final int removeDuplicateBestSites(Read r){
		int x=0;
		if(r.numSites()<2){return 0;}

		//Remove duplicate best sites that may exist as a result of realignment.
		final SiteScore ss1=r.topSite();
		for(int i=r.sites.size()-1; i>0; i--){
			SiteScore ss2=r.sites.get(i);
			if(ss1.chrom==ss2.chrom && ss1.strand==ss2.strand && ss1.start==ss2.start && ss1.stop==ss2.stop){
				if(!Shared.anomaly){
//					Shared.anomaly=true;
//					System.err.println("Ignoring anomalous duplicate site: "+"\n"+r.toText(false)+(r.mate==null ? "" : "\n"+r.mate.toText(false))+"\n");
					System.err.println("Ignoring anomalous duplicate site for rid="+r.numericID);
//					new Exception().printStackTrace(System.err);
				}
				r.sites.remove(i);
				x++;
			}else{break;}
		}
		return x;
	}
	
	protected final void removeUnmapped(ArrayList<Read> list){
		for(int i=0; i<list.size(); i++){
			Read r=list.get(i);
			if(r.numSites()==0){
				if(r.mate==null || r.mate.numSites()==0){
					list.set(i, null);
				}
			}
		}
	}
	
	protected final void removeBlacklisted(ArrayList<Read> list){
		for(int i=0; i<list.size(); i++){
			Read r=list.get(i);
			if(Blacklist.inBlacklist(r)){
				list.set(i, null);
			}
		}
	}
	
	protected final void removeMapped(ArrayList<Read> list){
		for(int i=0; i<list.size(); i++){
			Read r=list.get(i);
			if(r.numSites()==0){
				if(r.mate==null || r.mate.numSites()==0){
					list.set(i, null);
				}
			}
		}
	}
	
	
	public abstract int trimList(ArrayList<SiteScore> list, boolean retainPaired, int maxScore, boolean specialCasePerfect, int minSitesToRetain, int maxSitesToRetain);
	
	public final int trimListAdvanced(ArrayList<SiteScore> list, boolean retainPaired, boolean retainSemiperfect, int maxScore, boolean specialCasePerfect, 
			int minSitesToRetain, int maxSitesToRetain, boolean indexUsesExtendedScore, float thresh){
		if(list==null || list.size()==0){return -99999;}
		if(list.size()==1){return list.get(0).score;}
		
		final int highestScore;
		if(indexUsesExtendedScore){

			highestScore=Tools.trimSiteList(list, .6f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
			if(highestScore==maxScore && specialCasePerfect){
				Tools.trimSiteList(list, .94f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				if(list.size()>8){Tools.trimSiteList(list, .99f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
				return highestScore;
			}

		}else{
			highestScore=Tools.trimSiteList(list, .4f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
			thresh=thresh*0.5f;
		}

		int lim, lastScore=list.get(0).score;
		long area=lastScore;
		for(lim=1; lim<list.size(); lim++){
			SiteScore ss=list.get(lim);
			lastScore=ss.score;
			area=area+lastScore;
			if(lastScore<0 || (lastScore/(float)area)<thresh){break;}
		}
		lim=Tools.max(minSitesToRetain, lim);
		Tools.trimSitesBelowCutoff(list, lastScore, retainPaired, true, minSitesToRetain, maxSitesToRetain);
		
		return highestScore;
	}
	
	
	public abstract void scoreSlow(final ArrayList<SiteScore> list, final byte[] basesP, final byte[] basesM, 
			final int maxSwScore, final int maxImperfectSwScore);
	
	/** This is only for saving ambiguous xy which is now irrelevant */
	public final boolean processAmbiguous(ArrayList<SiteScore> list, boolean primary, boolean removeAmbiguous, int clearzone, boolean save_xy){
		if(!save_xy){return true;}
		assert(false) : "Needs to be redone with contig names.";

		assert(list.size()>1);
		boolean ambiguous=true;
//		if(save_xy && minChrom<=24 && maxChrom>=24){
//			int best=list.get(0).score;
//			
//			//Remove everything outside of the clearzone
//			for(int i=list.size()-1; i>0; i--){
//				assert(best>=list.get(i).score);
//				if(best-list.get(i).score>clearzone){
////					assert(i>1); //No longer true because of clearzone/clearzone2
//					list.remove(i);
//				}else{
////					assert(i>0); //Maybe no longer true because of clearzone/clearzone2
//					break;
//				}
//			}
//			
//			
//			assert(list.size()>1);
//			int Xcount=0;
//			int Ycount=0;
//			for(SiteScore ss : list){
//				assert(ss.score-list.get(0).score<=clearzone);
//				if(ss.chrom==23){
//					Xcount++;
//				}else if(ss.chrom==24){
//					Ycount++;
//				}
//			}
//			if(Xcount>1 || Ycount>2 || (Xcount+Ycount)<list.size()){
//				ambiguous=true;
//			}else{
//				ambiguous=false;
//				for(int i=list.size()-1; i>0; i--){list.remove(i);}
//				assert(list.size()==1);
//			}
//		}
		assert(list.size()>=1);

		if(ambiguous){
			assert(list.size()>1);
			if(removeAmbiguous){
				list.clear();
			}
		}
		
		return ambiguous;
	}
	
	
	public abstract void calcStatistics1(final Read r, final int maxSwScore, final int maxPossibleQuickScore);
	
	
	public abstract void calcStatistics2(final Read r, final int maxSwScore, final int maxPossibleQuickScore);
	
//	/** Assumes list is sorted */
//	public abstract void genMatchString(final Read r, final byte[] basesP, final byte[] basesM, final int maxImperfectSwScore, final int maxSwScore, boolean setSSScore, boolean recur);
	
	public abstract void processRead(Read r, final byte[] basesM);
	
	@Deprecated
	protected final boolean applyClearzone3_old(Read r, int CLEARZONE3, float INV_CLEARZONE3){
		
		assert(!r.paired()); //This is currently for unpaired reads
		if(!r.mapped() || r.ambiguous() || r.discarded() || r.numSites()<2){return false;}

		final int score1=r.topSite().slowScore;
		final int score2=r.sites.get(1).slowScore;
		final int score3=(r.sites.size()>2 ? r.sites.get(2).slowScore : -1);
		int dif=score1-score2;
		
		assert(r.mapScore==score1) : r.mapScore+", "+r.topSite().toText();

		assert(score1==r.mapScore);
		assert(score1>=score2) : "\n"+r.topSite().toText()+"\t<\t"+r.sites.get(1).toText()+"\n"+r.toText(false)+"\n";
		if(dif>=CLEARZONE3){return false;}

//		final int dif2=40+(CLEARZONE3-dif)/3;
//		final int dif2=(CLEARZONE3-dif)/2;
		int dif2=(CLEARZONE3-dif);
		
		float f=dif2*INV_CLEARZONE3;
		
		int sub=(dif2+2*(int)(f*dif2));
		
		if(score3!=-1){
			assert(score1>=score3);
			dif=score1-score3;
			assert(score1>=score3);
			if(dif<CLEARZONE3){
				dif2=(CLEARZONE3-dif);
				
				f=dif2*INV_CLEARZONE3;
				sub=sub+(dif2+2*(int)(f*dif2))/4;
				
//				sub=sub+(dif2)/2;
			}
		}
		
		for(SiteScore ss : r.sites){
			ss.score-=sub;
			ss.slowScore-=sub;
		}
		r.mapScore-=sub;
		return sub>0;
	}
	
	
	protected final boolean applyClearzone3(Read r, int CLEARZONE3, float INV_CLEARZONE3){
		
		assert(!r.paired()); //This is currently for unpaired reads
		final ArrayList<SiteScore> list=r.sites;
		if(!r.mapped() || r.ambiguous() || r.discarded() || list==null || list.size()<2){return false;}

		final int score1=list.get(0).slowScore;
		assert(r.mapScore==score1) : r.mapScore+", "+list.get(0).toText()+"\n"+r;
		
		float sub=0;
		final int max=Tools.min(CZ3_MULTS.length, list.size());
		for(int i=1; i<max; i++){
			final SiteScore ss2=list.get(i);
			assert(ss2!=null) : r;
			if(ss2!=null){
				if(i>2 && ss2.slowScore<list.get(i-1).slowScore){break;}
//				final int score2=list.get(i).slowScore;
//				final int dif=score1-score2;
//				if(dif>=CLEARZONE3){break;}
//				int dif2=(CLEARZONE3-dif);
//				float f=dif2*INV_CLEARZONE3;
//				sub+=(dif2+2*(f*dif2))*CZ3_MULTS[i];
				float f=calcCZ3_fraction(score1, ss2.slowScore, CLEARZONE3, INV_CLEARZONE3);
				if(f<=0){break;}
				sub+=(f*CZ3_MULTS[i]);
			}
		}
		assert(sub>=0);
		if(sub<=0){return false;}
		
		float sub2;
//		float asymptote=8f+0.0267f*r.bases.length;
		float asymptote=4f+0.03f*r.bases.length;
		sub=sub*1.8f;
		sub2=CLEARZONE3*((asymptote*sub)/(sub+asymptote));
//		sub2=CLEARZONE3*sub;
//		System.out.println("sub="+sub+", sub2="+sub2+", CLEARZONE3="+CLEARZONE3+", (5*sub)="+(5*sub)+", (sub+5*CLEARZONE3)="+(sub+5*CLEARZONE3));
		int subi=(int)(sub2+0.5f);
		if(subi>=r.mapScore-300){
			subi=r.mapScore-300;
		}
		if(subi<=0){return false;}
		
		for(SiteScore ss : list){
			ss.score-=subi;
			ss.slowScore-=subi;
		}
		r.mapScore-=subi;
		assert(r.mapScore>200);
		return true;
	}
	
	
//	protected float calcCZ3(int score1, int score2, int CLEARZONE3, float INV_CLEARZONE3){
//		
//		int dif=score1-score2;
//		if(dif>=CLEARZONE3){return 0;}
//		//Now dif is between 0 and CZ3
//
////		final int dif2=40+(CLEARZONE3-dif)/3;
////		final int dif2=(CLEARZONE3-dif)/2;
//		int dif2=(CLEARZONE3-dif); //dif2 is higher if the scores are closer.
//		
//		float f=dif2*INV_CLEARZONE3; //f ranges linearly from 1 (if the scores are identical) to 0 (when score2 is maximally below score1) 
//		
//		float f2=f*f;
//		float f7=(float)Math.pow(f, .7);
//		
////		return (dif2+2f*f*dif2+2f*Tools.min(f2,0.5f)*dif2);		
//		return (CLEARZONE3*f7+2f*f*dif2+2f*Tools.min(f2,0.5f)*dif2);
//	}
	
	
	protected float calcCZ3_fraction(int score1, int score2, int CLEARZONE3, float INV_CLEARZONE3){
		
		int dif=score1-score2;
		if(dif>=CLEARZONE3){return 0;}
		//Now dif is between 0 and CZ3

//		final int dif2=40+(CLEARZONE3-dif)/3;
//		final int dif2=(CLEARZONE3-dif)/2;
		int dif2=(CLEARZONE3-dif); //dif2 is higher if the scores are closer.
		
		float f=dif2*INV_CLEARZONE3; //f ranges linearly from 1 (if the scores are identical) to 0 (when score2 is maximally below score1) 
		
		float f2=f*f;
//		float f7=(float)Math.pow(f, .7);
		
//		return (dif2+2f*f*dif2+2f*Tools.min(f2,0.5f)*dif2);		
		return f+2f*f2+2f*f2*f;
	}
	
	/** Returns number of perfect pairs */
	public abstract int pairSiteScoresInitial(Read r, Read r2, boolean trim);


	
	

	protected static void pairSiteScoresFinal(Read r, Read r2, boolean trim, boolean setScore, int MAX_PAIR_DIST, int AVERAGE_PAIR_DIST, 
			boolean SAME_STRAND_PAIRS, boolean REQUIRE_CORRECT_STRANDS_PAIRS, int maxTrimSitesToRetain){
		
		if(r.sites!=null){
			for(SiteScore ss : r.sites){ss.pairedScore=0;}
		}
		if(r2.sites!=null){
			for(SiteScore ss : r2.sites){ss.pairedScore=0;}
		}
		
		if(r.numSites()<1 || r2.numSites()<1){return;}
		
		SiteScore.PCOMP.sort(r.sites);
		SiteScore.PCOMP.sort(r2.sites);

		int maxPairedScore1=-1;
		int maxPairedScore2=-1;
		
		
//		if(verbose){
//			System.out.println(r.list.size()+", "+r2.list.size());
//			System.out.println();
//			for(SiteScore ss : r.list){
//				System.out.println(ss.toText());
//			}
//			System.out.println();
//			for(SiteScore ss : r2.list){
//				System.out.println(ss.toText());
//			}
//			System.out.println();
//		}
		
		final float mult1=Tools.min(1/2f, Tools.max(1/4f, (r.bases.length/(4f*r2.bases.length))));
		final float mult2=Tools.min(1/2f, Tools.max(1/4f, (r2.bases.length/(4f*r.bases.length))));
		
		final int ilimit=r.sites.size()-1;
		final int jlimit=r2.sites.size()-1;
		
		final int outerDistLimit=(Tools.max(r.bases.length, r2.bases.length)*OUTER_DIST_MULT)/OUTER_DIST_DIV; //Minimum pairing distance
		final int expectedFragLength=AVERAGE_PAIR_DIST+r.bases.length+r2.bases.length;
		
		if(verboseS){
			System.err.println("**************************   PAIRING   ********************************");
			System.err.println("outerDistLimit="+outerDistLimit+", MAX_PAIR_DIST="+MAX_PAIR_DIST);
		}
		
		for(int i=0, j=0; i<=ilimit && j<=jlimit; i++){
			SiteScore ss1=r.sites.get(i);
			SiteScore ss2=r2.sites.get(j);
			
			while(j<jlimit && (ss2.chrom<ss1.chrom || (ss2.chrom==ss1.chrom && ss1.start-ss2.stop>MAX_PAIR_DIST))){
				j++;
//				if(verbose){System.err.println("a.Incrementing j->"+j);}
				ss2=r2.sites.get(j);
			}

			for(int k=j; k<=jlimit; k++){
				ss2=r2.sites.get(k);
				
				if(verboseS){
					System.err.println("Considering sites:\n"+ss1+"\n"+ss2);
				}

				if(ss2.chrom>ss1.chrom){break;}
				//				if(verbose){System.err.println("Same chrom");}
				if(ss2.start-ss1.stop>MAX_PAIR_DIST){break;}

				final int innerdist;
				final int outerdist;
				
				//assert(!SAME_STRAND_PAIRS) : "TODO";
				
				if(REQUIRE_CORRECT_STRANDS_PAIRS){
					if(ss1.strand!=ss2.strand){
						if(ss1.strand==Gene.PLUS){
							innerdist=ss2.start-ss1.stop;
							outerdist=ss2.stop-ss1.start;
						}else{
							innerdist=ss1.start-ss2.stop;
							outerdist=ss1.stop-ss2.start;
						}
					}else{
						if(ss1.start<=ss2.start){
							innerdist=ss2.start-ss1.stop;
							outerdist=ss2.stop-ss1.start;
						}else{
							innerdist=ss1.start-ss2.stop;
							outerdist=ss1.stop-ss2.start;
						}
					}
				}else{
					if(ss1.start<=ss2.start){
						innerdist=ss2.start-ss1.stop;
						outerdist=ss2.stop-ss1.start;
					}else{
						innerdist=ss1.start-ss2.stop;
						outerdist=ss1.stop-ss2.start;
					}
				}
				
				if(verboseS){
					System.err.println("innerdist="+innerdist+", outerdist="+outerdist);
				}
				
//				if(ss1.start<=ss2.start){
//					innerdist=ss2.start-ss1.stop;
//					outerdist=ss2.stop-ss1.start;
//				}else{
//					innerdist=ss1.start-ss2.stop;
//					outerdist=ss1.stop-ss2.start;
//				}
				assert(outerdist>=innerdist) : "outerdist<innerdist:\n"+innerdist+", "+outerdist+", "+ss1+", "+ss2;
				
				if(outerdist>=outerDistLimit && innerdist<=MAX_PAIR_DIST){

					boolean strandOK=((ss1.strand==ss2.strand)==SAME_STRAND_PAIRS);
					//					if(verbose){System.err.println("strandOK="+strandOK);}

					if(strandOK || !REQUIRE_CORRECT_STRANDS_PAIRS){

						int deviation=absdif(AVERAGE_PAIR_DIST, innerdist);

						final int pairedScore1;
						final int pairedScore2;
						if(strandOK){
							//							pairedScore1=ss1.score+(int)(ss2.score*mult1);
							//							pairedScore2=ss2.score+(int)(ss1.score*mult2);

							pairedScore1=ss1.score+1+
									Tools.max(1, (int)(ss2.score*mult1)-(((deviation)*ss2.score)/Tools.max(100,(10*expectedFragLength+100))));
							pairedScore2=ss2.score+1+
									Tools.max(1, (int)(ss1.score*mult2)-(((deviation)*ss1.score)/Tools.max(100,(10*expectedFragLength+100))));
							
							
						}else{//e.g. a junction
							pairedScore1=ss1.score+ss2.score/16;
							pairedScore2=ss2.score+ss1.score/16;
						}
						
						if(verboseS){
							System.err.println("strandOK="+strandOK+"\tpairedScore1="+pairedScore1+", pairedScore2="+pairedScore2);
							System.err.println("             \tscore1="+ss1.score+", score2="+ss2.score);
						}
						
						ss1.pairedScore=Tools.max(ss1.pairedScore, pairedScore1);
						ss2.pairedScore=Tools.max(ss2.pairedScore, pairedScore2);
						maxPairedScore1=Tools.max(ss1.score, maxPairedScore1);
						maxPairedScore2=Tools.max(ss2.score, maxPairedScore2);
						//						if(verbose){System.err.println("Paired:\nss1="+ss1.toText()+", ss2="+ss2.toText());}
					}
				}else{
					//					if(verbose){System.err.println("Out of range");}
				}
			}
			//			if(verbose){System.err.println("\nss1="+ss1.toText()+", ss2="+ss2.toText());}
			
		}
		
		if(setScore){
			for(SiteScore ss : r.sites){
				if(ss.pairedScore>ss.score){ss.score=ss.pairedScore;}
				else{assert(ss.pairedScore==0);}
			}
			for(SiteScore ss : r2.sites){
				if(ss.pairedScore>ss.score){ss.score=ss.pairedScore;}
				else{assert(ss.pairedScore==0);}
			}
		}
		
		if(trim){
//			Tools.trimSitesBelowCutoffInplace(r.list, (int)(maxPairedScore1*.95f), false);
//			Tools.trimSitesBelowCutoffInplace(r2.list, (int)(maxPairedScore2*.95f), false);
			Tools.trimSitesBelowCutoff(r.sites, (int)(maxPairedScore1*.95f), false, true, 1, maxTrimSitesToRetain);
			Tools.trimSitesBelowCutoff(r2.sites, (int)(maxPairedScore2*.95f), false, true, 1, maxTrimSitesToRetain);
		}
	}
	
	protected final boolean canPair(SiteScore ss1, SiteScore ss2, int len1, int len2, 
			boolean REQUIRE_CORRECT_STRANDS_PAIRS, boolean SAME_STRAND_PAIRS, int MAX_PAIR_DIST){
		if(ss1.chrom!=ss2.chrom){return false;}
		if(REQUIRE_CORRECT_STRANDS_PAIRS){
			boolean strandOK=((ss1.strand==ss2.strand)==SAME_STRAND_PAIRS);
			if(!strandOK){return false;}
		}
//		int dist=0;
//
//		if(ss1.start<=ss2.start){
//			dist=ss2.start-ss1.stop;
//		}else if(ss1.start>ss2.start){
//			dist=ss1.start-ss2.stop;
//		}
//
//		return (dist>=MIN_PAIR_DIST && dist<=MAX_PAIR_DIST);
		
//		final int outerDistLimit=MIN_PAIR_DIST+len1+len2;
//		final int outerDistLimit=(Tools.max(len1, len2)*(OUTER_DIST_MULT2))/OUTER_DIST_DIV;
		final int outerDistLimit=(Tools.max(len1, len2)*(OUTER_DIST_MULT))/OUTER_DIST_DIV;
		int innerdist=0;
		int outerdist=0;
		
		if(verboseS){
			System.err.println("canPair: outerDistLimit="+outerDistLimit);
		}
		
//		if(ss1.start<=ss2.start){
//			innerdist=ss2.start-ss1.stop;
//			outerdist=ss2.stop-ss1.start;
//		}else if(ss1.start>ss2.start){
//			innerdist=ss1.start-ss2.stop;
//			outerdist=ss1.stop-ss2.start;
//		}
//		assert(outerdist>=innerdist);
		
		//assert(!SAME_STRAND_PAIRS) : "TODO";
		
		if(REQUIRE_CORRECT_STRANDS_PAIRS){
			if(ss1.strand!=ss2.strand){
				if(ss1.strand==Gene.PLUS){
					innerdist=ss2.start-ss1.stop;
					outerdist=ss2.stop-ss1.start;
				}else{
					innerdist=ss1.start-ss2.stop;
					outerdist=ss1.stop-ss2.start;
				}
			}else{
				if(ss1.start<=ss2.start){
					innerdist=ss2.start-ss1.stop;
					outerdist=ss2.stop-ss1.start;
				}else{
					innerdist=ss1.start-ss2.stop;
					outerdist=ss1.stop-ss2.start;
				}
			}
		}else{
			if(ss1.start<=ss2.start){
				innerdist=ss2.start-ss1.stop;
				outerdist=ss2.stop-ss1.start;
			}else{
				innerdist=ss1.start-ss2.stop;
				outerdist=ss1.stop-ss2.start;
			}
		}
		
		return (outerdist>=outerDistLimit && innerdist<=MAX_PAIR_DIST);
	}
	
	
//	/** Returns the number of additional bases away that should be searched for slow align. 
//	 * This should probably be called between quickMap and slowAlign, only on 
//	 * sites where stop-start<=bases.length-1 */
//	public abstract void findTipDeletions(final Read r, final byte[] basesP, final byte[] basesM, final int maxSwScore, final int maxImperfectScore);
//
//	public abstract boolean findTipDeletions(SiteScore ss, final byte[] bases, final int maxImperfectScore, boolean lookRight, boolean lookLeft);
	
	
	/** Returns the number of additional bases away that should be searched for slow align. 
	 * This should probably be called between quickMap and slowAlign, only on 
	 * sites where stop-start<=bases.length-1 */
	protected final int findTipDeletionsRight(final byte[] bases, final int chrom,
			int originalStop, int searchDist, int tiplen){
		ChromosomeArray cha=Data.getChromosome(chrom);
		byte[] ref=cha.array;
		if(originalStop<cha.minIndex+tiplen-1){return 0;} //fail
		
		int minMismatches=tiplen;
		int bestStart=originalStop;
		
		final int tipCoord=bases.length-1;
		
		int lastMismatch=0;
		int originalMismatches=0;
		int contig=0;
		for(int i=0; i<tiplen && contig<5; i++){
			if(bases[tipCoord-i]!=ref[originalStop-i]){
				originalMismatches++;
				lastMismatch=i;
				contig=0;
			}else{
				contig++;
			}
		}
//		System.err.print("!");
		if(originalMismatches<3){return 0;}
		minMismatches=originalMismatches;
		tiplen=lastMismatch+1;
		assert(tiplen>1);
		if(tiplen<4){return 0;}
//		System.err.println("Tiplen="+tiplen+", mismatches="+originalMismatches);
//		System.err.print("* ");
		
		searchDist=Tools.min(searchDist, 30*originalMismatches);
		int lastIndexToStart=Tools.min(ref.length-1, originalStop+searchDist);
		for(int start=originalStop+1; start<=lastIndexToStart && minMismatches>0; start++){
//			System.err.print("_");
			int mismatches=0;
			for(int j=0; j<tiplen && mismatches<minMismatches; j++){
				if(bases[tipCoord-j]!=ref[start-j]){
					mismatches++;
				}
			}
//			System.err.print(mismatches+" ");
			if(mismatches<minMismatches){
				bestStart=start;
				minMismatches=mismatches;
			}
		}
//		System.err.println("\noriginalStop+1:"+(originalStop+1)+"\nlastIndexToStart:"+(lastIndexToStart)+"\ntiplen: "+tiplen+"\noriginalMismatches: "+originalMismatches+"\nminMismatches: "+minMismatches+"\n");
		if(minMismatches>2 || originalMismatches-minMismatches<2){
			return 0;
		}
//		System.err.println(" $$$ ");
		return bestStart-originalStop;
	}
	
	
	/** Returns the number of additional bases away that should be searched for slow align. 
	 * This should probably be called between quickMap and slowAlign, only on 
	 * sites where stop-start<=bases.length-1 */
	protected final int findTipDeletionsLeft(final byte[] bases, final int chrom,
			final int originalStart, int searchDist, int tiplen){
		ChromosomeArray cha=Data.getChromosome(chrom);
		byte[] ref=cha.array;
		if(originalStart+tiplen>=ref.length){return 0;} //fail
		
		if(cha.minIndex>=originalStart){return 0;} //fail
		
		int minMismatches=tiplen;
		int bestStart=originalStart;
		
		int lastMismatch=0;
		int originalMismatches=0;
		int contig=0;
		for(int i=0; i<tiplen && contig<5; i++){
			if(bases[i]!=ref[originalStart+i]){
				originalMismatches++;
				lastMismatch=i;
				contig=0;
			}else{
				contig++;
			}
		}
//		System.err.print("!");
		if(originalMismatches<3){return 0;}
		minMismatches=originalMismatches;
		tiplen=lastMismatch+1;
		assert(tiplen>1);
		if(tiplen<4){return 0;}
//		System.err.println("Tiplen="+tiplen+", mismatches="+originalMismatches);
//		System.err.print("* ");
		
		searchDist=Tools.min(searchDist, 16+16*originalMismatches+8*tiplen);
		int lastIndexToStart=Tools.max(cha.minIndex, originalStart-searchDist);
		for(int start=originalStart-1; start>=lastIndexToStart && minMismatches>0; start--){
//			System.err.print("_");
			int mismatches=0;
			for(int j=0; j<tiplen && mismatches<minMismatches; j++){
				if(bases[j]!=ref[start+j]){
					mismatches++;
				}
			}
//			System.err.print(mismatches+" ");
			if(mismatches<minMismatches){
				bestStart=start;
				minMismatches=mismatches;
			}
		}
//		System.err.println("\noriginalStop+1:"+(originalStop+1)+"\nlastIndexToStart:"+(lastIndexToStart)+"\ntiplen: "+tiplen+"\noriginalMismatches: "+originalMismatches+"\nminMismatches: "+minMismatches+"\n");
		if(minMismatches>2 || originalMismatches-minMismatches<2){
			return 0;
		}
//		System.err.println(" $$$ ");
		return originalStart-bestStart;
	}
	
	
//	public abstract void rescue(Read anchor, Read loose, byte[] basesP, byte[] basesM, int searchDist);
	
	
//	public abstract void slowRescue(final byte[] bases, SiteScore ss, final int maxScore, final int maxImperfectScore, 
//			boolean findTipDeletionsRight, boolean findTipDeletionsLeft);
	
	
	/** Assumes bases/colors are already on the correct strand */
	public final SiteScore quickRescue(final byte[] bases, final int chrom, final byte strand, final int loc, final int searchDist,
			final boolean searchRight, final int idealStart, final int maxAllowedMismatches, int POINTS_MATCH, int POINTS_MATCH2){
		if(bases==null || bases.length<10){return null;}
		ChromosomeArray cha=Data.getChromosome(chrom);
		byte[] ref=cha.array;
		
		int lowerBound, upperBound;
		if(searchRight){
			lowerBound=Tools.max(cha.minIndex, loc);
			upperBound=Tools.min(ref.length-bases.length, loc+searchDist);
		}else{
			lowerBound=Tools.max(cha.minIndex, loc-searchDist);
			upperBound=Tools.min(ref.length-bases.length, loc);
		}

//		int minMismatches=(int)(bases.length*.6f); //Default: .75f.  Lower numbers are faster with lower quality.
		int minMismatches=maxAllowedMismatches+1;
		//For situations like RNASEQ with lots of deletions, a higher value of at least .75 should be used.
		
		int maxContigMatches=0;
		int bestScore=0;
		int bestStart=-1;
		int bestAbsdif=Integer.MAX_VALUE;
		
		if(searchRight){
			for(int start=lowerBound; start<=upperBound/* && minMismatches>0*/; start++){
				int mismatches=0;
				int contig=0;
				int currentContig=0;
				for(int j=0; j<bases.length && mismatches<=minMismatches; j++){
					final byte c=bases[j], r=ref[start+j];
					if(c!=r || c=='N'){
						mismatches++;
						contig=Tools.max(contig, currentContig);
						currentContig=0;
					}else{
						currentContig++;
					}
				}
				
				int score=(bases.length-mismatches)+contig;
				int absdif=absdif(start, idealStart);
				if(mismatches<=minMismatches && (score>bestScore || (score==bestScore && absdif<bestAbsdif))){
					bestStart=start;
					minMismatches=mismatches;
					maxContigMatches=contig;
					bestScore=score;
					bestAbsdif=absdif;
					if(mismatches==0){upperBound=Tools.min(upperBound, idealStart+absdif);}
//					assert(upperBound>=start && lowerBound<=start);
//					assert(upperBound>=idealStart);
//					assert(lowerBound<=idealStart);
				}
			}
		}else{
			for(int start=upperBound; start>=lowerBound/* && minMismatches>0*/; start--){
				int mismatches=0;
				int contig=0;
				int currentContig=0;
				for(int j=0; j<bases.length && mismatches<=minMismatches; j++){
					final byte c=bases[j], r=ref[start+j];
					if(c!=r || c=='N'){
						mismatches++;
						contig=Tools.max(contig, currentContig);
						currentContig=0;
					}else{
						currentContig++;
					}
				}

				int score=(bases.length-mismatches)+contig;
				int absdif=absdif(start, idealStart);
				if(mismatches<=minMismatches && (score>bestScore || (score==bestScore && absdif<bestAbsdif))){
					bestStart=start;
					minMismatches=mismatches;
					maxContigMatches=contig;
					bestScore=score;
					bestAbsdif=absdif;
					if(mismatches==0){lowerBound=Tools.max(lowerBound, idealStart-absdif);}
//					assert(upperBound>=start && lowerBound<=start);
//					assert(upperBound>=idealStart);
//					assert(lowerBound<=idealStart);
				}
			}
		}
		
		if(bestStart<0){return null;}
		
		//These scores are dummies and will not quite match the normally generated scores.
		final int scoreOut;
		if(USE_AFFINE_SCORE){
			scoreOut=POINTS_MATCH+(POINTS_MATCH2*(bases.length-1-minMismatches));
		}else{
			scoreOut=maxContigMatches+(BASE_HIT_SCORE*(bases.length-minMismatches));
		}
		
		SiteScore ss=new SiteScore(chrom, strand, bestStart, bestStart+bases.length-1, 0, scoreOut);
		ss.setPerfect(bases);
		ss.rescued=true;
		ss.slowScore=minMismatches; //TODO: Clear this field later!
		return ss;
	}
	
	
	/** Assumes bases/colors are already on the correct strand */
	protected final int[] quickerRescue(final byte[] bases, final int chrom, int loc, final int searchDist){
		ChromosomeArray cha=Data.getChromosome(chrom);
		byte[] ref=cha.array;
		if(loc<cha.minIndex){loc=cha.minIndex;}
		
		int lastIndexToStart=loc+searchDist-1;
		final int limit=Tools.min(lastIndexToStart, ref.length-bases.length)+1;

		int minMismatches=bases.length;
		int bestStart=-1;
		for(int start=loc; start<limit && minMismatches>0; start++){
			int mismatches=0;
			for(int j=0; j<bases.length && mismatches<minMismatches; j++){
				if(bases[j]!=ref[start+j]){
					mismatches++;
				}
			}
			if(mismatches<minMismatches){
				bestStart=start;
				minMismatches=mismatches;
			}
		}
		
		return new int[] {bestStart, bestStart+bases.length-1, minMismatches};
	}
	
	
	public abstract void processReadPair(final Read r, final byte[] basesM1, final byte[] basesM2);
	
	/** TODO: Iterate through loop backwards when removing sites.
	 * @param r
	 * @param DONT_OUTPUT_UNMAPPED_READS
	 * @param SAM_OUT
	 * @param EXPECTED_LEN_LIMIT
	 * @return Number of sites removed
	 */
	protected final int removeOutOfBounds(Read r, boolean DONT_OUTPUT_UNMAPPED_READS, boolean SAM_OUT, int EXPECTED_LEN_LIMIT){
//		assert(false) : DONT_OUTPUT_UNMAPPED_READS+", "+SAM_OUT+", "+EXPECTED_LEN_LIMIT;
		ArrayList<SiteScore> ssl=r.sites;
		if(ssl==null){return 0;}
		int initial=ssl.size();
		for(int i=0; i<ssl.size(); i++){
			SiteScore ss=ssl.get(i);
//			System.out.println("Estimated greflen: "+GapTools.calcGrefLen(ss.start, ss.stop, ss.gaps));
			ChromosomeArray cha=Data.getChromosome(ss.chrom);
			int max=cha.maxIndex;
			if(ss.start<0 || ss.stop>max){
				ssl.remove(i);
				i--;
				ss=null;
			}else if(/*DONT_OUTPUT_UNMAPPED_READS && */SAM_OUT){
				if(!Data.isSingleScaffold(ss.chrom, ss.start, ss.stop)){
					//TODO: Attempt to trim instead of removing
					ssl.remove(i);
					i--;
					ss=null;
				}
			}
			if(ss!=null){
				int expectedLen=GapTools.calcGrefLen(ss);
				if(expectedLen>=EXPECTED_LEN_LIMIT){
					//TODO: Alternately, I could kill the site.
					ss.stop=ss.start+Tools.min(r.bases.length+40, EXPECTED_LEN_LIMIT);
					if(ss.gaps!=null){GapTools.fixGaps(ss);}
				}
			}
		}
		
//		System.out.println("Estimated greflen: "+GapTools.calcGrefLen(r.start, r.stop, r.gaps));
//		assert(false);
		
		return initial-ssl.size();
	}
	
	protected static final int forbidSelfMapping(ArrayList<SiteScore> ssl, SiteScore original){
//		assert(original!=null);
		if(ssl==null || ssl.isEmpty() || original==null){return 0;}
		int removed=0;
		for(int i=0; i<ssl.size(); i++){
			SiteScore ss=ssl.get(i);
			if(ss.overlaps(original, true)){
				ssl.set(i, null);
				removed++;
			}
		}
		if(removed>0){Tools.condenseStrict(ssl);}
		return removed;
	}
	

	/** Generate a score penalty based on the presence of errors near the read tips. */
	public static int calcTipScorePenalty(final Read r, final int maxScore, final int tiplen){
		if(!r.mapped() || r.match==null || r.bases.length<2*tiplen){return 0;}
		
		int points=0;
		final byte[] match=r.match;
		final byte[] bases=r.bases;
		final int last=r.bases.length-1;
		byte prev='m';
		for(int i=0, cpos=0; cpos<=tiplen; i++){
			byte b=match[i];
			if(b=='m'){
				cpos++;
			}else if(b=='D'){
				if(prev!='D'){points+=2*(tiplen+2-cpos);}
			}else if(b=='N' || b=='C'){
				points+=(tiplen+2-cpos);
				cpos++;
			}else{
				if(Character.isDigit(b)){
					r.match=Read.toLongMatchString(r.match);
					return calcTipScorePenalty(r, maxScore, tiplen);
				}
				assert(b=='I' || b=='S') : ((char)b)+"\n"+new String(match)+"\n"+new String(bases)+"\n";
				points+=2*(tiplen+2-cpos);
				cpos++;
			}
			prev=b;
		}
		
		prev='m';
		for(int i=match.length-1, cpos=0; cpos<=tiplen; i--){
			byte b=match[i];
			if(b=='m'){
				cpos++;
			}else if(b=='D'){
				if(prev!='D'){points+=2*(tiplen+2-cpos);}
			}else if(b=='N' || b=='C'){
				points+=(tiplen+2-cpos);
				cpos++;
			}else{
				assert(b=='I' || b=='S');
				points+=2*(tiplen+2-cpos);
				cpos++;
			}
			prev=b;
		}
		
		byte b=bases[0];
		//homopolymer tip penalty
		if(b!='N' && b==bases[1]){
			for(int i=2; i<=tiplen && bases[i]==b; i++){points++;}
		}
		
		//homopolymer tip penalty
		b=bases[last];
		if(b!='N' && b==bases[last-1]){
			for(int i=last-2; i>=(last-tiplen) && bases[i]==b; i--){points++;}
		}
		
		//Did not seem to help
//		int hits=r.list.get(0).hits;
//		float desired=Tools.min(6, bases.length/12f);
//		if(hits<desired){points+=20*(1-(hits/desired));}
		
		if(points<1){return 0;}
//		points=Tools.min(points, 40);
		
		float asymptote=80;
		float f=((asymptote*points)/(points+asymptote));
		
		int penalty=(int)(f*.0022f*maxScore);
		int maxPenalty=r.mapScore-maxScore/10;
		if(maxPenalty<=0){return 0;}
		return Tools.min(penalty, maxPenalty);
		
//		final int len=7;
//		int dist1=len+1, dist2=len+1;
//		for(int i=0; i<=len; i++){
//			if(r.match[i]!='m'){
//				dist1=i;
////				System.out.println("dist1="+dist1+"\n"+new String(r.match));
//				break;
//			}
//		}
//		for(int i=0, tip=r.match.length-1; i<=len; i++){
//			if(r.match[tip-i]!='m'){
//				dist2=i;
////				System.out.println("dist2="+dist2+"\n"+new String(r.match));
//				break;
//			}
//		}
//		float penalty=0;
//		if(dist1<len){
//			penalty+=(len-dist1+1)*.005f*maxScore;
//		}
//		if(dist2<len){
//			penalty+=(len-dist2+1)*.005f*maxScore;
//		}
//		return (int)penalty;
	}
	
	
	public static void applyScorePenalty(Read r, int penalty){
		if(penalty>0){
			r.mapScore-=penalty;
			for(SiteScore ss : r.sites){
				ss.score-=penalty;
				ss.slowScore-=penalty;
				ss.pairedScore-=penalty;
			}
		}
	}
	
	
	/** {group of correct hit (or -1), size of correct group, number of groups, 
	 * number of elements, correctScore, maxScore, size of top group, num correct, firstElementCorrect, 
	 * firstElementCorrectLoose, firstGroupCorrectLoose} */
	protected int[] calcCorrectness(Read r, int thresh){
		//assume sorted.
		ArrayList<SiteScore> ssl=r.sites;
		
		if(ssl==null || ssl.isEmpty()){
			return new int[] {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		}

		SiteScore original=r.originalSite;
		assert((original==null) != (r.synthetic()));
		if(original==null){
			original=ssl.get(0);
		}
		
		int group=0;
		int correctGroup=-1;
		int groupSize=0;
		int correctGroupSize=-1;
		int prevScore=Integer.MAX_VALUE;
		int sizeOfTopGroup=0;
		SiteScore correct=null;

		int firstElementCorrect=0;
		int firstElementCorrectLoose=0;
		int firstGroupCorrectLoose=0;
		
		int numCorrect=0;
		
		for(int i=0; i<ssl.size(); i++){
			SiteScore ss=ssl.get(i);
			if(ss.score==ssl.get(0).score){sizeOfTopGroup++;}
			
			if(prevScore!=ss.score){
				assert(prevScore>ss.score || (AMBIGUOUS_RANDOM && r.ambiguous()) || r.mate!=null) : "i="+i+", r="+r;
				
				if(correctGroup==group){
					correctGroupSize=groupSize;
				}
				
				group++;
				groupSize=0;
				prevScore=ss.score;
			}
			groupSize++;
			

//			boolean b=isCorrectHit(ss, original.chrom, original.strand, original.start, 1, thresh);
			boolean b=isCorrectHit(ss, original.chrom, original.strand, original.start, original.stop, thresh);
			boolean b2=isCorrectHitLoose(ss, original.chrom, original.strand, original.start, original.stop, thresh+20);
			if(b){
				if(i==0){firstElementCorrect=1;}
				numCorrect++;
				if(correct==null){
					correct=ss;
					correctGroup=group;
				}
			}
			if(b2){
				if(i==0){firstElementCorrectLoose=1;}
				if(group==0){firstGroupCorrectLoose=1;}
			}
		}
		if(correctGroup==group){
			correctGroupSize=groupSize;
		}

		assert(correctGroup!=0 && correctGroup<=group);
		assert(group<=ssl.size());
		assert(sizeOfTopGroup>0 && sizeOfTopGroup<=ssl.size());
		assert((correctGroup>0) == (correctGroupSize>0));
		return new int[] {correctGroup, correctGroupSize, group, ssl.size(), 
				correct==null ? 0 : correct.score, ssl.get(0).score, sizeOfTopGroup, numCorrect, firstElementCorrect, 
						firstElementCorrectLoose, firstGroupCorrectLoose};
	}
	
	
	public static final boolean isCorrectHit(SiteScore ss, int trueChrom, byte trueStrand, int trueStart, int trueStop, int thresh){
//		boolean b=(ss.chrom==trueChrom && ss.strand==trueStrand);
		if(ss.chrom!=trueChrom || ss.strand!=trueStrand){return false;}

		assert(ss.stop>ss.start) : ss.toText()+", "+trueStart+", "+trueStop;
		assert(trueStop>trueStart) : ss.toText()+", "+trueStart+", "+trueStop;

		return (absdif(ss.start, trueStart)<=thresh && absdif(ss.stop, trueStop)<=thresh);
//		return (absdif(ss.start, trueStart)<=thresh || absdif(ss.stop, trueStop)<=thresh);
		
//		if(absdif(ss.start, trueStart)<=thresh){return true;}
//		if(absdif(ss.stop, trueStop)<=thresh){return true;}
//		return false;
		
//		if(absdif(ss.start, trueStart)>thresh){return false;}
//		if(absdif(ss.stop, trueStop)>thresh){return false;}
//		return true;
	}
	
	
	public static final boolean isCorrectHitLoose(SiteScore ss, int trueChrom, byte trueStrand, int trueStart, int trueStop, int thresh){
//		boolean b=(ss.chrom==trueChrom && ss.strand==trueStrand);
		if(ss.chrom!=trueChrom || ss.strand!=trueStrand){return false;}

		assert(ss.stop>ss.start) : ss.toText()+", "+trueStart+", "+trueStop;
		assert(trueStop>trueStart) : ss.toText()+", "+trueStart+", "+trueStop;
		
		return (absdif(ss.start, trueStart)<=thresh || absdif(ss.stop, trueStop)<=thresh);
		
//		if(absdif(ss.start, trueStart)<=thresh){return true;}
//		if(absdif(ss.stop, trueStop)<=thresh){return true;}
//		return false;
		
//		if(absdif(ss.start, trueStart)>thresh){return false;}
//		if(absdif(ss.stop, trueStop)>thresh){return false;}
//		return true;
	}
	
	protected static final byte[] makePerfectMatchString(int len){
		byte[] r=new byte[len];
		Arrays.fill(r, (byte)'m');
		return r;
	}
	
	protected static final int absdif(int a, int b){
		return a>b ? a-b : b-a;
	}
	
	/** Returns maximum read length supported by this mapper */
	public abstract int maxReadLength();
	
	/** Ensure top site is congruent with read */
	protected static final boolean checkTopSite(Read r){
		if(!r.mapped()){return true;}
		if(r.numSites()==0){return false;}
		SiteScore ss=r.topSite();
		if(ss==null){return false;}
		boolean b=(ss.start==r.start) && (ss.stop==r.stop) && (ss.strand==r.strand()) && (ss.chrom==r.chrom) && (ss.match==r.match);
		assert(b) : "\nread="+r+"\nmate="+r.mate+"\nss="+ss+"\n"+(ss==null ? "ss is null" : 
			((ss.start==r.start)+", "+(ss.stop==r.stop)+", "+(ss.strand==r.strand())+", "+(ss.chrom==r.chrom)+", "+(ss.match==r.match))+"\nlist="+r.sites);
		return b;
	}
	
	
	protected static final int removeLongIndels(ArrayList<SiteScore> list, int maxlen){
		if(list==null || list.size()<1){return 0;}
		int removed=0;
		for(int i=list.size()-1; i>=0; i--){
			SiteScore ss=list.get(i);
			if(hasLongIndel(ss.match, maxlen)){
				list.remove(i);
				removed++;
			}
		}
		return removed;
	}
	
	protected static final boolean hasLongIndel(byte[] match, int maxlen){
		if(match==null || match.length<maxlen){return false;}
		byte prev='0';
		int len=0;
		for(byte b : match){
			if(b=='D' || b=='I'){
				if(b==prev){len++;}
				else{len=1;}
				if(len>maxlen){return true;}
			}else{
				len=0;
			}
			prev=b;
		}
		return false;
	}
	
	/** TODO */
	final void processReadSplit(Read r, byte[] basesM, int minlen, int maxlen){
		assert(minlen>=KEYLEN && maxlen>=minlen) : KEYLEN+", "+maxlen+", "+minlen;
		int len=r.bases==null ? 0 : r.bases.length;
		if(len<=maxlen){
			processRead(r, basesM);
			return;
		}
		ArrayList<Read> subreads=r.split(minlen, maxlen);
	}
	
	public final synchronized boolean finished(){return finished;}
	
	public final synchronized boolean working(){return !finished;}
	
	final synchronized void finish(){
		assert(!finished);
		finished=true;
		notifyAll();
	}
	
	private boolean finished=false;
	
	private static final float[] CZ3_MULTS=new float[] {0f, 1f, .75f, 0.5f, 0.25f, 0.125f, 0.0625f};
	
	/*--------------------------------------------------------------*/
	
	/** Input read source. */
	protected final ConcurrentReadStreamInterface cris;

	
	/** All reads go here. <br>  
	 * If outputunmapped=false, omit unmapped single reads and double-unmapped paired reads. */
	protected final RTextOutputStream3 outStream;
	/** All mapped reads (and half-mapped pairs) go here except reads that only map to the blacklist. */
	protected final RTextOutputStream3 outStreamMapped;
	/** All unmapped reads (and double-unmapped pairs) go here. */
	protected final RTextOutputStream3 outStreamUnmapped;
	/** All reads (and half-mapped pairs) that map best to the blacklist go here. */
	protected final RTextOutputStream3 outStreamBlack;
	
	
	/*--------------------------------------------------------------*/
	
	
	public final String MSA_TYPE;
	final MSA msa;
	final TranslateColorspaceRead tcr;
	public final ReadStats readstats;
	public final int POINTS_MATCH, POINTS_MATCH2;
	public final int KEYLEN;
	
	protected final boolean PERFECTMODE; //Only look for perfect matches
	protected final boolean SEMIPERFECTMODE; //Only look for perfect and semiperfect matches
	protected final boolean FORBID_SELF_MAPPING; //Do not allow reads to map to their official origin.  Allows you to find next-best matches (when supported)
	protected final boolean RCOMP_MATE; //Reverse-complement mate prior to mapping
	/** True if this thread should generate a match string for the best match */
	protected final boolean MAKE_MATCH_STRING;
	
	protected final boolean DONT_OUTPUT_UNMAPPED_READS;
	protected final boolean DONT_OUTPUT_BLACKLISTED_READS;
	protected final boolean PRINT_SECONDARY_ALIGNMENTS;
	protected final boolean QUICK_MATCH_STRINGS;
	protected final boolean USE_SS_MATCH_FOR_PRIMARY=true;

	protected final int MAX_SITESCORES_TO_PRINT;
	
	/** Scores below the (max possible alignment score)*(MINIMUM_ALIGNMENT_SCORE_RATIO) will be discarded.
	 * Default: 0.4 for synthetic data. */
	protected final float MINIMUM_ALIGNMENT_SCORE_RATIO;
	protected final float MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE;
	protected final float MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED;

	protected final float keyDensity;
	protected final float maxKeyDensity;
	protected final float minKeyDensity;
	protected final int maxDesiredKeys;
	
	/*--------------------------------------------------------------*/

	final int CLEARZONE1e;
	
	/*--------------------------------------------------------------*/
	
	final int MIN_APPROX_HITS_TO_KEEP;
	final boolean USE_EXTENDED_SCORE;
	public static final boolean GENERATE_BASE_SCORES_FROM_QUALITY=AbstractIndex.GENERATE_BASE_SCORES_FROM_QUALITY;
	final int BASE_HIT_SCORE;
	final int BASE_KEY_HIT_SCORE;
	final boolean USE_AFFINE_SCORE;
	final int EXPECTED_LEN_LIMIT;
	final int MAX_INDEL;
	
	final boolean TRIM_LIST;
	final int TIP_DELETION_SEARCH_RANGE;
	final boolean FIND_TIP_DELETIONS;
	final int ALIGN_COLUMNS;
	
	/*--------------------------------------------------------------*/
	
	
	/** Deprecated.  Must be set to false.  Reads and index are in SOLiD colorspace. */
	protected final boolean colorspace;
	/** Use dynamic programming slow-alignment phase to increase quality.  Program may not run anymore if this is disabled. */
	protected final boolean SLOW_ALIGN;
	/** Produce local alignments instead of global alignments */
	protected final boolean LOCAL_ALIGN;
	/** Discard reads with ambiguous alignments (consider them unmapped). */
	protected final boolean AMBIGUOUS_TOSS;
	/** Choose a random site for reads with ambiguous alignments. */
	protected final boolean AMBIGUOUS_RANDOM;
	/** Output all sites for reads with ambiguous alignments. */
	protected final boolean AMBIGUOUS_ALL;		
	/** Quality-trim left side of reads before mapping. */
	protected final boolean TRIM_LEFT;
	/** Quality-trim right side of reads before mapping. */
	protected final boolean TRIM_RIGHT;
	/** Undo quality trimming after mapping. */
	protected final boolean UNTRIM;
	/** Trim until 2 consecutive bases are encountered with at least this quality. */
	protected final byte TRIM_QUAL;
	/** Don't trim reads to be shorter than this */
	protected final int TRIM_MIN_LENGTH=30;
	/** Distance cutoff for classifying a read as loosely correct */
	protected final int THRESH;
	/** Semi-deprecated.  Minimum chrom to index or load. */
	protected final int minChrom;
	/** Semi-deprecated.  Maximum chrom to index or load. */
	protected final int maxChrom;
	/** Disallow sites that do not have at least k consecutive matching bases. */
	protected final int KFILTER;
	
	
	/** When reads are not in valid pairing orientation, eliminate (mark unmapped) the lower-scoring read. */
	protected final boolean KILL_BAD_PAIRS;
	/** For human genome, map ambiguous reads in the PAR to the X chromosome. */
	protected final boolean SAVE_AMBIGUOUS_XY;
	/** Deprecated.  Must be set to true. */
	protected final boolean GEN_MATCH_FAST=true;
	/** For colorspace reads, translate to base space before outputting them. */
	protected final boolean translateToBaseSpace;
	
	/** Padding for dynamic-programming slow alignment. */
	protected final int SLOW_ALIGN_PADDING;
	/** Padding for dynamic-programming slow alignment for rescued reads (which typically may need more padding). */
	protected final int SLOW_RESCUE_PADDING;
	/** If a site is unpaired, search nearby for a possible site for the other read. */
	protected final boolean DO_RESCUE;
	/** Forbid alignments with indels longer than MAX_INDEL */
	protected final boolean STRICT_MAX_INDEL;
	/** Bandwidth of banded MSA */
	protected final int BANDWIDTH;
	
	protected final boolean PAIRED;
	protected final boolean REQUIRE_CORRECT_STRANDS_PAIRS;
	protected final boolean SAME_STRAND_PAIRS;
	
	/*--------------------------------------------------------------*/

	protected int AVERAGE_PAIR_DIST=100;

	/** Extra padding for when slow alignment fails. */
	protected int EXTRA_PADDING=10;

	protected final boolean GENERATE_KEY_SCORES_FROM_QUALITY;
	
	/*--------------------------------------------------------------*/
	
	protected static boolean CALC_STATISTICS=true;
	protected static int MIN_PAIR_DIST=-160;
	protected static int MAX_PAIR_DIST=32000;
	/** IMPORTANT!!!!  This option causes non-deterministic output. */
	protected static final boolean DYNAMIC_INSERT_LENGTH=true;
	/** Counts undefined bases. */
	protected static final boolean DISCARD_MOSTLY_UNDEFINED_READS=true;

	protected static final byte TIP_DELETION_MIN_QUALITY=6;
	protected static final byte TIP_DELETION_AVG_QUALITY=14;
	protected static final int TIP_DELETION_MAX_TIPLEN=8;

	protected static final int OUTER_DIST_MULT=14;
//	protected static final int OUTER_DIST_MULT2=OUTER_DIST_MULT-1;
	protected static final int OUTER_DIST_DIV=32;
	
	protected static long SKIP_INITIAL=0;
	
	protected static boolean OUTPUT_PAIRED_ONLY=false;
	
//	static{if(OUTER_DIST_MULT2<1){throw new RuntimeException();}}
	
	/*--------------------------------------------------------------*/
	
	public int totalNumCorrect1=0;
	public int totalNumIncorrect1=0;
	public int totalNumIncorrectPrior1=0;
	public int totalNumCapturedAllCorrect1=0;
	public int totalNumCapturedAllCorrectTop1=0;
	public int totalNumCapturedAllCorrectOnly1=0;
	
	public int totalNumCorrect2=0;
	public int totalNumIncorrect2=0;
	public int totalNumIncorrectPrior2=0;
	public int totalNumCapturedAllCorrect2=0;
	public int totalNumCapturedAllCorrectTop2=0;
	public int totalNumCapturedAllCorrectOnly2=0;
	
	/*--------------------------------------------------------------*/

	public boolean verbose=false;
	public static final boolean verboseS=false;
	
	public long readsUsed=0;
	public long readsUsed2=0;
	public long numMated=0;
	public long badPairs=0;
	public long innerLengthSum=0;
	public long outerLengthSum=0;
	public long insertSizeSum=0;
	public long keysUsed=0;
	public long basesUsed=0; //basesUsed and basesAtQuickmap are identical
	public long basesAtQuickmap=0; //basesUsed and basesAtQuickmap are identical
	public long syntheticReads=0;

	public int mapped1=0;
	public int mappedRetained1=0;
	public int rescuedP1=0;
	public int rescuedM1=0;
	public int truePositiveP1=0;
	public int truePositiveM1=0;
	public int falsePositive1=0;
	public int totalCorrectSites1=0;

	public int firstSiteCorrectP1=0;
	public int firstSiteCorrectM1=0;
	public int firstSiteIncorrect1=0;
	public int firstSiteCorrectLoose1=0;
	public int firstSiteIncorrectLoose1=0;
	public int firstSiteCorrectPaired1=0;
	public int firstSiteCorrectSolo1=0;
	public int firstSiteCorrectRescued1=0;

	public long matchCountS1=0;
	public long matchCountI1=0;
	public long matchCountD1=0;
	public long matchCountM1=0;
	public long matchCountN1=0;
	
	
	public int perfectHit1=0; //Highest quick score is max quick score
	public int uniqueHit1=0; //Only one hit has highest score
	public int correctUniqueHit1=0; //unique highest hit on answer site 
	public int correctMultiHit1=0;  //non-unique highest hit on answer site 
	public int correctLowHit1=0;  //hit on answer site, but not highest scorer 
	public int noHit1=0;
	
	/** Number of perfect hit sites found */
	public int perfectHitCount1=0;
	/** Number of sites found that are perfect except for no-ref */
	public int semiPerfectHitCount1=0;
	
	
	public int perfectMatch1=0; //Highest slow score is max slow score
	public int semiperfectMatch1=0;
	
	public int ambiguousBestAlignment1=0;

	public long initialSiteSum1=0;
	public long postTrimSiteSum1=0;
	public long postRescueSiteSum1=0;
	public long siteSum1=0;
	public long topSiteSum1=0;
	
	public long lowQualityReadsDiscarded1=0;
	
	public int mapped2=0;
	public int mappedRetained2=0;
	public int rescuedP2=0;
	public int rescuedM2=0;
	public int truePositiveP2=0;
	public int truePositiveM2=0;
	public int falsePositive2=0;
	public int totalCorrectSites2=0;

	public int firstSiteCorrectP2=0;
	public int firstSiteCorrectM2=0;
	public int firstSiteIncorrect2=0;
	public int firstSiteCorrectLoose2=0;
	public int firstSiteIncorrectLoose2=0;
	public int firstSiteCorrectPaired2=0;
	public int firstSiteCorrectSolo2=0;
	public int firstSiteCorrectRescued2=0;
	
	public long matchCountS2=0;
	public long matchCountI2=0;
	public long matchCountD2=0;
	public long matchCountM2=0;
	public long matchCountN2=0;
	
	public int perfectHit2=0; //Highest quick score is max quick score
	public int uniqueHit2=0; //Only one hit has highest score
	public int correctUniqueHit2=0; //unique highest hit on answer site 
	public int correctMultiHit2=0;  //non-unique highest hit on answer site 
	public int correctLowHit2=0;  //hit on answer site, but not highest scorer 
	public int noHit2=0;
	
	/** Number of perfect hit sites found */
	public int perfectHitCount2=0;
	/** Number of sites found that are perfect except for no-ref */
	public int semiPerfectHitCount2=0;

	public int perfectMatch2=0; //Highest slow score is max slow score
	public int semiperfectMatch2=0;
	
	public int ambiguousBestAlignment2=0;

	public long initialSiteSum2=0;
	public long postTrimSiteSum2=0;
	public long postRescueSiteSum2=0;
	public long siteSum2=0;
	public long topSiteSum2=0;
	
	public long lowQualityReadsDiscarded2=0;
	
	/*--------------------------------------------------------------*/
	
	int idmodulo;
}
