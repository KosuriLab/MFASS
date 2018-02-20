package align2;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLongArray;

import stream.SiteScore;

import dna.Data;

public final class Tools {
	
	/** Return true if the user seems confused */
	public static boolean parseHelp(String[] args){
		if(args==null || args.length==0 || (args.length==1 && args[0]==null)){return true;}
		if(args.length>1){return false;}
		final String s=args[0].toLowerCase();
		return s.equals("-h") || s.equals("-help") || s.equals("--help") 
				|| s.equals("-version") || s.equals("--version") || s.equals("?") || s.equals("-?") || (s.equals("help") && !new File(s).exists());
	}
	
	/** Checks for permission to overwrite files, and output name collisions. */
	public static boolean testOutputFiles(boolean overwrite, boolean allowDuplicates, String...args){
		if(args==null || args.length==0){return true;}
		HashSet<String> set=new HashSet<String>(args.length*2);
		int terms=0;
		for(String s : args){
			if(s!=null){
				if(isOutputFileName(s)){
					terms++;
					
					if(!overwrite && new File(s).exists()){
						assert(overwrite) : "File "+s+" exists and overwrite=false";
						return false;
					}
					
					if(!allowDuplicates && set.contains(s)){
						assert(false) : "Duplicate file "+s+" was specified for multiple output streams.";
						return false;
					}
					
					set.add(s);
				}		
			}
		}
		return true;
	}
	
	public static final boolean canWrite(String s, boolean overwrite){
		if(isNullFileName(s) || isSpecialOutputName(s)){return true;}
		File f=new File(s);
		if(f.exists()){return overwrite && f.canWrite();}
		return true;
	}
	
//	public static final boolean outputDestinationExists(String s){
//		if(isNullFileName(s)){return false;}
//		if(isSpecialOutputName(s)){return false;}
//		File f=new File(s);
//		return f.exists();
//	}
	
	public static final boolean isOutputFileName(String s){
		return !(isNullFileName(s) || isSpecialOutputName(s)); 
	}
	
	public static final boolean isNullFileName(String s){
		if(s==null || s.equalsIgnoreCase("null") || s.equalsIgnoreCase("none")){return true;}
		for(int i=0; i<s.length(); i++){
			if(!Character.isWhitespace(s.charAt(i))){return false;}
		}
		return true;
	}
	
	public static final boolean isSpecialOutputName(String s){
		if(s==null){return false;}
		s=s.toLowerCase();
		return s.equals("stdout") || s.equals("stderr") || s.equals("standardout") || s.equals("standarderr")
				|| s.equals("/dev/null") || s.startsWith("stdout.") || s.startsWith("stderr.");
	}
	
	public static final boolean isSpecialInputName(String s){
		if(s==null){return false;}
		s=s.toLowerCase();
		return s.equals("stdin") || s.equals("standardin") || s.startsWith("stdin.");
	}
	
	public static final boolean canRead(String s){
		if(s==null){return false;}
		if(isSpecialInputName(s)){return true;}
		File f=new File(s);
		return f.canRead();
	}
	
	/** Removes null elements by shrinking the list.  May change list order. */
	public static final <X> int condense(ArrayList<X> list){
		if(list==null || list.size()==0){return 0;}
		int removed=0;
		
		for(int i=list.size()-1; i>0; i--){
			if(list.get(i)==null){
				removed++;
				X last=list.get(list.size()-1);
				list.set(i, last);
				list.remove(list.size()-1);
			}
		}
		return removed;
	}
	
	/** Removes null elements by shrinking the list.  Will not change list order. */
	public static final <X> int condenseStrict(ArrayList<X> list){
		if(list==null || list.size()==0){return 0;}
		int removed=0;
		
		int insertPos=0;
		for(int i=0; i<list.size(); i++){
			X x=list.get(i);
			if(x!=null){
				if(insertPos!=i){
					assert(insertPos<i);
					while(list.get(insertPos)!=null){insertPos++;}
					assert(insertPos<i && list.get(insertPos)==null) : insertPos+", "+i; //slow, temporary
					list.set(i, null);
					list.set(insertPos, x);
				}
				insertPos++;
			}else{
				removed++;
			}
		}
		for(int i=0; i<removed; i++){
			X x=list.remove(list.size()-1);
			assert(x==null);
		}
		return removed;
	}
	
	/** Creates a new list without null elements. */
	public static final <X> ArrayList<X> condenseNew(ArrayList<X> list){
		ArrayList<X> temp=new ArrayList<X>(list.size());
		for(X x : list){
			if(x!=null){temp.add(x);}
		}
		return temp;
	}
	
	//This should also be correct.  I'm not sure which is faster.
//	/** Removes null elements by shrinking the list.  Will not change list order. */
//	public static final <X> int condenseStrict(ArrayList<X> list){
//		if(list==null || list.size()==0){return 0;}
//		int removed=0;
//		int last=0;
//		
//		for(int i=0; i<list.size(); i++){
//			X x=list.get(i);
//			if(x==null){
//				removed++;
//			}else{
//				while(last<i && list.get(last)!=null){last++;}
//				assert(last==i || list.get(last)==null);
//				if(last!=i){
//					assert(last<i);
//					list.set(last, x);
//					list.set(i, null);
//				}
//			}
//		}
//		for(int i=0; i<removed; i++){
//			X x=list.remove(list.size()-1);
//			assert(x==null);
//		}
//		return removed;
//	}
	
	
	public static final int calcMedianDistance(int[] array){
		if(array==null || array.length<2){return 500000000;}
		int[] dif=new int[array.length-1];
		for(int i=0; i<array.length; i++){
			dif[i]=(array[i+1]-array[i]);
		}
		Arrays.sort(dif);
		return dif[dif.length/2];
	}
	
	
	

//	public static final int trimSiteList(ArrayList<SiteScore> ssl, float fractionOfMax, boolean retainPaired){
////		assert(false);
//		if(ssl==null || ssl.size()==0){return -999999;}
//		if(ssl.size()==1){return ssl.get(0).score;}
//		int maxScore=-999999;
//		for(SiteScore ss : ssl){
//			maxScore=Tools.max(maxScore, ss.score);
//		}
//		
//		int cutoff=(int) (maxScore*fractionOfMax);
//		trimSitesBelowCutoff(ssl, cutoff, retainPaired);
////		trimSitesBelowCutoffInplace(ssl, cutoff);
//		return maxScore;
//	}
	
	/** minSitesToRetain should be set to 1 if the list is not sorted by score (for efficiency of removal).  Otherwise, it can be higher. */
	public static final int trimSiteList(ArrayList<SiteScore> ssl, float fractionOfMax, boolean retainPaired, boolean retainSemiperfect, 
			int minSitesToRetain, int maxSitesToRetain){
//		assert(false);
		if(ssl==null || ssl.size()==0){return -999999;}
		if(ssl.size()==1){return ssl.get(0).score;}
		int maxScore=-999999;
		
		if(minSitesToRetain>1 && minSitesToRetain<ssl.size()){
			assert(inOrder(ssl));
			maxScore=ssl.get(0).score;
		}else{
			for(SiteScore ss : ssl){
				maxScore=Tools.max(maxScore, ss.score);
			}
		}
		
		int cutoff=(int) (maxScore*fractionOfMax);
		trimSitesBelowCutoff(ssl, cutoff, retainPaired, retainSemiperfect, minSitesToRetain, maxSitesToRetain);
		return maxScore;
	}
	
	/** minSitesToRetain should be set to 1 if the list is not sorted by score.  Otherwise, it can be higher. */
	public static final void trimSiteListByMax(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired, boolean retainSemiperfect, 
			int minSitesToRetain, int maxSitesToRetain){
//		assert(false);
		if(ssl==null || ssl.size()==0){return;}
		if(ssl.size()==1){return;}
		
		trimSitesBelowCutoff(ssl, cutoff, retainPaired, retainSemiperfect, minSitesToRetain, maxSitesToRetain);
	}
	
	public static final <X extends Comparable<? super X>> boolean inOrder(ArrayList<X> list){
		if(list==null || list.size()<2){return true;}
		for(int i=1; i<list.size(); i++){
			X xa=list.get(i-1);
			X xb=list.get(i);
			if(xa.compareTo(xb)>0){return false;}
		}
		return true;
	}
	

	
	public static final int mergeDuplicateSites(ArrayList<SiteScore> list, boolean doAssertions, boolean mergeDifferentGaps){
		if(list==null || list.size()<2){return 0;}
		Collections.sort(list, SiteScore.PCOMP);
		
		int removed=0;
		
		SiteScore a=list.get(0);
		for(int i=1; i<list.size(); i++){
			SiteScore b=list.get(i);
			if(a.positionalMatch(b, true)){
				
				if(doAssertions){
					if(!(a.perfect==b.perfect || 
							(a.perfect && (a.score>b.score || a.slowScore>b.slowScore)))){
						throw new RuntimeException("\n"+SiteScore.header()+"\n"+a.toText()+"\n"+b.toText()+"\n");
					}

					assert(a.perfect==b.perfect || 
							(a.perfect && (a.score>b.score || a.slowScore>b.slowScore))) : 
								"\n"+SiteScore.header()+"\n"+a.toText()+"\n"+b.toText()+"\n";
				}
				
				a.score=max(a.score, b.score);
				a.slowScore=max(a.slowScore, b.slowScore);
				a.pairedScore=max(a.pairedScore, b.pairedScore);
				a.perfect=(a.perfect || b.perfect);
				if(a.pairedScore>0 && a.pairedScore<=a.score){a.pairedScore=a.score+1;}
				
				removed++;
				list.set(i, null);
			}else if(mergeDifferentGaps && a.positionalMatch(b, false)){ //Same outermost boundaries, different gaps
				
				SiteScore better=null;
				if(a.score!=b.score){
					better=(a.score>b.score ? a : b);
				}else if(a.slowScore!=b.slowScore){
					better=(a.slowScore>b.slowScore ? a : b);
				}else if(a.pairedScore!=b.pairedScore){
					better=(a.pairedScore>b.pairedScore ? a : b);
				}else{
					better=a;
				}
				
				a.score=max(a.score, b.score);
				a.slowScore=max(a.slowScore, b.slowScore);
				a.pairedScore=max(a.pairedScore, b.pairedScore);
				a.perfect=(a.perfect || b.perfect);
				if(a.pairedScore>0 && a.pairedScore<=a.score){a.pairedScore=a.score+1;}
				a.gaps=better.gaps;
				
				removed++;
				list.set(i, null);
			}
			else{
				a=b;
			}
		}

//		if(removed>0){condense(list);}
		if(removed>0){condenseStrict(list);}
		return removed;
	}
	

	
	public static final int subsumeOverlappingSites(ArrayList<SiteScore> list, boolean subsumeIfOnlyStartMatches, boolean subsumeInexact){
		if(list==null || list.size()<2){return 0;}
		Collections.sort(list, SiteScore.PCOMP);
		
		int removed=0;
		
		
		for(int i=0; i<list.size(); i++){
			SiteScore a=list.get(i);
			
			assert(a==null || !a.perfect || a.semiperfect);
			
			boolean overlappingA=true;
			if(a!=null){
				for(int j=i+1; overlappingA && j<list.size(); j++){
					SiteScore b=list.get(j);
					assert(b==null || !b.perfect || b.semiperfect);
					if(b!=null){
						overlappingA=(a.chrom==b.chrom && b.start<a.stop && b.stop>a.start);
						if(overlappingA && a.strand==b.strand){
							
							SiteScore better=null;
							if(a.perfect!=b.perfect){
								better=a.perfect ? a : b;
							}if(a.semiperfect!=b.semiperfect){
								better=a.semiperfect ? a : b;
							}else if(a.score!=b.score){
								better=(a.score>b.score ? a : b);
							}else if(a.slowScore!=b.slowScore){
								better=(a.slowScore>b.slowScore ? a : b);
							}else if(a.pairedScore!=b.pairedScore){
								better=(a.pairedScore>b.pairedScore ? a : b);
							}else if(a.pairedScore!=b.pairedScore){
								better=(a.quickScore>b.quickScore ? a : b);
							}else{
								better=a;
							}
							
//							if((a.perfect && b.perfect) || (a.semiperfect && b.semiperfect)){
							if(a.semiperfect && b.semiperfect){
								if(a.start==b.start || a.stop==b.stop){
									list.set(i, better);
									list.set(j, null);
									removed++;
									a=better;
								}else{
									//retain both of them
								}
							}else if(a.perfect || b.perfect){
								list.set(i, better);
								list.set(j, null);
								removed++;
								a=better;
							}else if(a.semiperfect || b.semiperfect){
								if(a.start==b.start && a.stop==b.stop){
									list.set(i, better);
									list.set(j, null);
									removed++;
									a=better;
								}else{
									//retain both of them
								}
							}else if(subsumeInexact || (a.start==b.start && (subsumeIfOnlyStartMatches || a.stop==b.stop))){
								assert(!a.semiperfect && !a.perfect && !b.semiperfect && !b.perfect);
								a.start=min(a.start, b.start);
								a.stop=max(a.stop, b.stop);
								a.score=max(a.score, b.score);
								a.slowScore=max(a.slowScore, b.slowScore);
								a.pairedScore=max(a.pairedScore, b.pairedScore);
								a.quickScore=max(a.quickScore, b.quickScore);
								if(a.pairedScore>0 && a.pairedScore<=a.score){a.pairedScore=a.score+1;}
								a.gaps=better.gaps;//Warning!  Merging gaps would be better; this could cause out-of-bounds.
								//TODO: Test for a subsumption length limit.
								list.set(j, null);
								removed++;
							}
						}
					}
				}
			}
		}
		
//		if(removed>0){condense(list);}
		if(removed>0){condenseStrict(list);}
		return removed;
	}
	

	
	public static final int removeOverlappingSites(ArrayList<SiteScore> list, boolean requireAMatchingEnd){
		if(list==null || list.size()<2){return 0;}
		Collections.sort(list, SiteScore.PCOMP);
		
		int removed=0;
		
		
		for(int i=0; i<list.size(); i++){
			SiteScore a=list.get(i);
			boolean overlappingA=true;
			if(a!=null){
				for(int j=i+1; overlappingA && j<list.size(); j++){
					SiteScore b=list.get(j);
					if(b!=null){
						overlappingA=(a.chrom==b.chrom && b.start<a.stop && b.stop>a.start);
						if(overlappingA && a.strand==b.strand){
							
							SiteScore better=null;
							if(a.perfect!=b.perfect){
								better=a.perfect ? a : b;
							}else if(a.score!=b.score){
								better=(a.score>b.score ? a : b);
							}else if(a.slowScore!=b.slowScore){
								better=(a.slowScore>b.slowScore ? a : b);
							}else if(a.pairedScore!=b.pairedScore){
								better=(a.pairedScore>b.pairedScore ? a : b);
							}else if(a.pairedScore!=b.pairedScore){
								better=(a.quickScore>b.quickScore ? a : b);
							}else{
								better=a;
							}
							
							if(a.start==b.start && a.stop==b.stop){
								list.set(i, better);
								list.set(j, null);
								a=better;
								removed++;
							}else if(a.start==b.start || a.stop==b.stop){ //In this case they cannot both be perfect
								list.set(i, better);
								list.set(j, null);
								a=better;
								removed++;
							}else if(!requireAMatchingEnd && a.score!=b.score){
								list.set(i, better);
								list.set(j, null);
								a=better;
								removed++;
							}
						}
					}
				}
			}
		}
		
//		if(removed>0){condense(list);}
		if(removed>0){condenseStrict(list);}
		return removed;
	}
	

	
	/** Returns the number of sitescores in the list within "thresh" of the top score.  Assumes list is sorted descending. 
	 * This is used to determine whether a mapping is ambiguous. */
	public static final int countTopScores(ArrayList<SiteScore> list, int thresh){
		assert(thresh>=0) : thresh;
		if(list==null || list.isEmpty()){return 0;}
		int count=1;
		final SiteScore ss=list.get(0);
		final int limit=ss.score-thresh;
		
		for(int i=1; i<list.size(); i++){
			SiteScore ss2=list.get(i);
			if(ss2.score<limit){break;}
			if(ss.start!=ss2.start && ss.stop!=ss2.stop){ //Don't count mappings to the same location
				count++;
			}
		}
		return count;
	}
	

	
	/** Assumes list is sorted by NON-PAIRED score. 
	 * Returns number removed. */
	public static final int removeLowQualitySitesPaired(ArrayList<SiteScore> list, int maxSwScore, float multSingle, float multPaired){
		if(list==null || list.size()==0){return 0;}
		
		assert(multSingle>=multPaired);
		
		int initialSize=list.size();
		final int swScoreThresh=(int)(maxSwScore*multSingle); //Change low-quality alignments to no-hits.
		final int swScoreThreshPaired=(int)(maxSwScore*multPaired);
		if(list.get(0).score<swScoreThreshPaired){list.clear(); return initialSize;}
		
		for(int i=list.size()-1; i>=0; i--){
			SiteScore ss=list.get(i);
			assert(ss.score==ss.slowScore);
			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
			if(ss.pairedScore>0){
				assert(ss.pairedScore>ss.quickScore || ss.pairedScore>ss.slowScore) : ss;
				if(ss.slowScore<swScoreThreshPaired){list.remove(i);}
			}else{
				assert(ss.pairedScore==0) : ss.toText();
				if(ss.slowScore<swScoreThresh){list.remove(i);}
			}
		}
		
		return initialSize-list.size();
	}
	

	
//	/** Assumes list is sorted by NON-PAIRED score. 
//	 * Returns number removed. */
//	public static final int removeLowQualitySitesUnpaired(ArrayList<SiteScore> list, int maxSwScore, float multSingle){
//		if(list==null || list.size()==0){return 0;}
//		
//		int initialSize=list.size();
//		final int swScoreThresh=(int)(maxSwScore*multSingle); //Change low-quality alignments to no-hits.
//		if(list.get(0).score<swScoreThresh){list.clear(); return initialSize;}
//		
////		for(int i=list.size()-1; i>=0; i--){
//		for(int i=list.size()-1; i>1; i--){
//			SiteScore ss=list.get(i);
//			assert(ss.score==ss.slowScore);
//			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
//			assert(ss.pairedScore==0) : ss.toText();
//			if(ss.slowScore<swScoreThresh){list.remove(i);}
//		}
//		
//		return initialSize-list.size();
//	}

	
	/** Assumes list is sorted by NON-PAIRED score. 
	 * Returns number removed. */
	public static final int removeLowQualitySitesUnpaired(ArrayList<SiteScore> list, int thresh){
		if(list==null || list.size()==0){return 0;}
		
		int initialSize=list.size();
		if(list.get(0).score<thresh){list.clear(); return initialSize;}
		
//		for(int i=list.size()-1; i>=0; i--){
		for(int i=list.size()-1; i>1; i--){
			SiteScore ss=list.get(i);
			assert(ss.score==ss.slowScore);
			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
			assert(ss.pairedScore==0) : ss.toText();
			if(ss.slowScore<thresh){list.remove(i);}
		}
		
		return initialSize-list.size();
	}
	

	
	/** Assumes list is sorted by NON-PAIRED score. 
	 * Returns number removed. */
	public static final int removeLowQualitySitesPaired2(ArrayList<SiteScore> list, int maxSwScore, float multSingle, float multPaired, int expectedSites){
		if(list==null || list.size()==0){return 0;}
		
		assert(multSingle>=multPaired);
		
		int initialSize=list.size();
		final int swScoreThresh=(int)(maxSwScore*multSingle); //Change low-quality alignments to no-hits.
		final int swScoreThreshPaired=(int)(maxSwScore*multPaired);
		final int swScoreThresh2=(int)(maxSwScore*multSingle*1.2f);
		final int swScoreThreshPaired2=(int)(maxSwScore*multPaired*1.1f);
		if(list.get(0).score<swScoreThreshPaired){list.clear(); return initialSize;}
		final int nthBest=list.get(Tools.min(list.size(), expectedSites)-1).score-maxSwScore/64;
		
		for(int i=list.size()-1, min=expectedSites*2; i>min; i--){
			if(list.get(i).slowScore>=nthBest){break;}
			list.remove(i);
		}
		
		for(int i=list.size()-1; i>=0; i--){
			SiteScore ss=list.get(i);
			assert(ss.score==ss.slowScore);
			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
			if(ss.pairedScore>0){
				int thresh=(i>=expectedSites ? swScoreThreshPaired2 : swScoreThreshPaired);
				assert(ss.pairedScore>ss.quickScore || ss.pairedScore>ss.slowScore) : ss;
				if(ss.slowScore<thresh){list.remove(i);}
			}else{
				int thresh=(i>=expectedSites ? swScoreThresh2 : swScoreThresh);
				assert(ss.pairedScore==0) : ss.toText();
				if(ss.slowScore<thresh){list.remove(i);}
			}
		}
		
		return initialSize-list.size();
	}
	

	
	/** Assumes list is sorted by NON-PAIRED score. 
	 * Returns number removed. 
	 * This has a couple of changes (like potentially removing the second-best site) that make it applicable to SKIMMER not MAPPER.
	 * */
	public static final int removeLowQualitySitesUnpaired2(ArrayList<SiteScore> list, int maxSwScore, float multSingle, int expectedSites){
		if(list==null || list.size()==0){return 0;}
		
		for(int i=expectedSites/2; i<list.size(); i++){
			if(list.get(i).perfect){expectedSites++;}
		}
		
		int initialSize=list.size();
		final int swScoreThresh=(int)(maxSwScore*multSingle); //Change low-quality alignments to no-hits.
		final int swScoreThresh2=(int)(maxSwScore*multSingle*1.2f); //Change low-quality alignments to no-hits.
		if(list.get(0).score<swScoreThresh){list.clear(); return initialSize;}
		final int nthBest=list.get(Tools.min(list.size(), expectedSites)-1).score-maxSwScore/64;
		
		for(int i=list.size()-1, min=expectedSites*2; i>min; i--){
			if(list.get(i).slowScore>=nthBest){break;}
			list.remove(i);
		}
		
//		for(int i=list.size()-1; i>=0; i--){
		for(int i=list.size()-1; i>=1; i--){
			SiteScore ss=list.get(i);
			assert(ss.score==ss.slowScore);
			assert(i==0 || ss.slowScore<=list.get(i-1).slowScore) : "List is not sorted by singleton score!";
			assert(ss.pairedScore==0) : ss.toText();
			int thresh=(i>=expectedSites ? swScoreThresh2 : swScoreThresh);
			if(ss.slowScore<thresh){list.remove(i);}
		}
		
		return initialSize-list.size();
	}
	
	
//	public static final void trimSitesBelowCutoff(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired){
//		trimSitesBelowCutoff(ssl, cutoff, retainPaired, 1);
//	}
	
	
//	public static final void trimSitesBelowCutoff(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired, int minSitesToRetain){
////		assert(false);
//		assert(minSitesToRetain>=1);
//		if(ssl==null || ssl.size()<minSitesToRetain){return;}
//		
//		ArrayList<SiteScore> ssl2=new ArrayList<SiteScore>(ssl.size());
//		for(SiteScore ss : ssl){
//			if(ss.score>=cutoff || (retainPaired && ss.pairedScore>0)){
//				ssl2.add(ss);
//			}
//		}
//		
////		Collections.sort(ssl2);
////		System.err.println("Cutoff: "+cutoff);
////		for(SiteScore ss : ssl2){
////			System.err.print("("+ss.chrom+", "+ss.score+"), ");
////		}
////		System.err.println();
//		
//		if(ssl2.size()==ssl.size()){return;}
////		System.err.println("cutoff: "+cutoff+",\tsize: "+ssl.size()+" -> "+ssl2.size());
//		ssl.clear();
//		ssl.addAll(ssl2);
//	}
	
	
	public static final void trimSitesBelowCutoff(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired, boolean retainSemiperfect,
			int minSitesToRetain, int maxSitesToRetain){
//		assert(false);
		assert(minSitesToRetain>=1);
		assert(maxSitesToRetain>minSitesToRetain);
		if(ssl==null || ssl.size()<=minSitesToRetain){return;}
		while(ssl.size()>maxSitesToRetain){ssl.remove(ssl.size()-1);}
		
		int removed=0;
		final int maxToRemove=ssl.size()-minSitesToRetain;
		
		assert(minSitesToRetain==1 || inOrder(ssl));
		
		if(retainPaired){
			for(int i=ssl.size()-1; i>=0; i--){
				SiteScore ss=ssl.get(i);
				if(!retainSemiperfect || !ss.semiperfect){
					if(ss.score<cutoff && ss.pairedScore<=0){
						ssl.set(i, null);
						removed++;
						if(removed>=maxToRemove){
							assert(removed==maxToRemove);
							break;
						}
					}
				}
			}
		}else{
			for(int i=ssl.size()-1; i>=0; i--){
				SiteScore ss=ssl.get(i);
				if(!retainSemiperfect || !ss.semiperfect){
					if(ss.score<cutoff){
						ssl.set(i, null);
						removed++;
						if(removed>=maxToRemove){
							assert(removed==maxToRemove);
							break;
						}
					}
				}
			}
		}
		
		if(removed>0){
			condenseStrict(ssl);
		}
		assert(ssl.size()>=minSitesToRetain);
	}
	
	//Messes up order
//	public static final void trimSitesBelowCutoffInplace(ArrayList<SiteScore> ssl, int cutoff, boolean retainPaired){
////		assert(false);
//		if(ssl==null || ssl.size()<2){return;}
//		
//		for(int i=0; i<ssl.size(); i++){
//			SiteScore ss=ssl.get(i);
//			if(ss.score<cutoff && (!retainPaired || ss.pairedScore==0)){
//				SiteScore temp=ssl.remove(ssl.size()-1);
//				if(i<ssl.size()){
//					ssl.set(i, temp);
//					i--;
//				}
//			}
//		}
//	}
	
	
	public static boolean equals(byte[] a, byte[] b){
		if(a==b){return true;}
		if(a==null || b==null){return false;}
		if(a.length!=b.length){return false;}
		for(int i=0; i<a.length; i++){
			if(a[i]!=b[i]){return false;}
		}
		return true;
	}

	public static int compare(byte[] a, byte[] b){
		if(a==b){return 0;}
		if(a==null){return -1;}
		if(b==null){return 1;}
		int lim=min(a.length, b.length);
		for(int i=0; i<lim; i++){
			if(a[i]!=b[i]){return a[i]-b[i];}
		}
		return a.length-b.length;
	}

	public static int sum(byte[] array){
		int x=0;
		for(byte y : array){x+=y;}
		return x;
	}
	
	public static long sum(short[] array){
		long x=0;
		for(short y : array){x+=y;}
		return x;
	}
	
	public static long sum(int[] array){
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}
	
	public static long sum(int[] array, int from, int to){
		long x=0;
		for(int i=from; i<=to; i++){x+=array[i];}
		return x;
	}
	
	public static long sum(long[] array){
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}
	
	public static long sum(AtomicIntegerArray array){
		long x=0;
		for(int i=0; i<array.length(); i++){x+=array.get(i);}
		return x;
	}
	
	public static long sum(AtomicLongArray array){
		long x=0;
		for(int i=0; i<array.length(); i++){x+=array.get(i);}
		return x;
	}
	
	public static int min(int[] array){
		int min=Integer.MAX_VALUE;
		for(int y : array){if(y<min){min=y;}}
		return min;
	}
	
	public static int intSum(int[] array){
		int x=0;
		for(int y : array){x+=y;}
		return x;
	}
	
	public static void reverseInPlace(final byte[] array){
		if(array==null){return;}
		final int max=array.length/2, last=array.length-1;
		for(int i=0; i<max; i++){
			byte temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}
	
	public static void reverseInPlace(final int[] array){
		if(array==null){return;}
		final int max=array.length/2, last=array.length-1;
		for(int i=0; i<max; i++){
			int temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}
	
	public static byte[] reverseAndCopy(final byte[] array){
//		if(array==null){return null;}
//		byte[] copy=Arrays.copyOf(array, array.length);
//		reverseInPlace(copy);
//		return copy;
		return reverseAndCopy(array, null);
	}
	
	public static int[] reverseAndCopy(final int[] array){
//		if(array==null){return null;}
//		int[] copy=Arrays.copyOf(array, array.length);
//		reverseInPlace(copy);
//		return copy;
		return reverseAndCopy(array, null);
	}
	
	public static byte[] reverseAndCopy(final byte[] array, byte[] out){
		if(array==null){
			assert(out==null);
			return null;
		}
		if(out==null){out=new byte[array.length];}
		assert(array.length==out.length && array!=out);
		for(int i=0, last=array.length-1; i<array.length; i++){out[i]=array[last-i];}
		return out;
	}
	
	public static int[] reverseAndCopy(final int[] array, int[] out){
		if(array==null){
			assert(out==null);
			return null;
		}
		if(out==null){out=new int[array.length];}
		assert(array.length==out.length && array!=out);
		for(int i=0, last=array.length-1; i<array.length; i++){out[i]=array[last-i];}
		return out;
	}
	
	public static void cullHighFreqEntries(int[][] data, float fractionToExclude){
		if(fractionToExclude<=0){return;}
		int[] count=new int[data.length];
		
		long numBases=0;
		
		for(int i=0; i<data.length; i++){
			count[i]=(data[i]==null ? 0 : data[i].length);
			numBases+=count[i];
		}
		
		int numIndicesToRemove=((int)(numBases*fractionToExclude));
		
		Arrays.sort(count);
		
		for(int i=1; i<count.length; i++){
			assert(count[i]>=count[i-1]) : "\n\ncount["+i+"]="+count[i]+"\ncount["+(i-1)+"]="+count[i-1]+"\n";
		}
		
		int pos=count.length-1;
		for(int sum=0; pos>1 && sum<numIndicesToRemove; pos--){
			sum+=count[pos];
		}
		int maxLengthToKeep2=count[pos];
		
		for(int i=0; i<data.length; i++){
			if(data[i]!=null && data[i].length>maxLengthToKeep2){data[i]=null;}
		}
	}
	
	public static int findLimitForHighFreqEntries(int[][] data, float fractionToExclude){
		if(fractionToExclude<=0){return Integer.MAX_VALUE;}
		int[] count=new int[data.length];
		
		long numBases=0;
		
		for(int i=0; i<data.length; i++){
			count[i]=(data[i]==null ? 0 : data[i].length);
			numBases+=count[i];
		}
		
		int numIndicesToRemove=((int)(numBases*fractionToExclude));
		
		Arrays.sort(count);
		
		for(int i=1; i<count.length; i++){
			assert(count[i]>=count[i-1]) : "\n\ncount["+i+"]="+count[i]+"\ncount["+(i-1)+"]="+count[i-1]+"\n";
		}
		
		int pos=count.length-1;
		for(int sum=0; pos>1 && sum<numIndicesToRemove; pos--){
			sum+=count[pos];
		}
		int maxLengthToKeep2=count[pos];
		
		return maxLengthToKeep2;
	}
	
	public static void cullClumpyEntries(final int[][] data, final int maxDist, final int minLength, final float fraction){
		
		long total=0;
		long removedSites=0;
		long removedKeys=0;
		
		if(maxDist<=0){return;}
		for(int i=0; i<data.length; i++){
			int[] array=data[i];
			total+=(array==null ? 0 : array.length);
			if(array!=null && array.length>=minLength){
				if(isClumpy(array, maxDist, fraction)){
					removedSites+=array.length;
					removedKeys++;
					data[i]=null;
				}
			}
		}

//		System.err.println("Removed\t"+removedSites+"\t/ "+total+"\tsites," +
//				" or "+String.format("%.4f", (removedSites*100f/total))+"%");
//		System.err.println("Removed\t"+removedKeys+"\t/ "+data.length+"\tkeys," +
//				" or  "+String.format("%.4f", (removedKeys*100f/data.length))+"%");
		
	}
	
	public static HashSet<Integer> banClumpyEntries(final int[][] data, final int maxDist, final int minLength, final float fraction){
		
		HashSet<Integer> set=new HashSet<Integer>(128);
		
		long total=0;
		long removedSites=0;
		long removedKeys=0;
		
		if(maxDist<=0){return set;}
		
		for(int i=0; i<data.length; i++){
			int[] array=data[i];
			total+=(array==null ? 0 : array.length);
			if(array!=null && array.length>=minLength){
				if(isClumpy(array, maxDist, fraction)){
					removedSites+=array.length;
					removedKeys++;
					set.add(i);
				}
			}
		}

//		System.err.println("Banned\t"+removedSites+"\t/ "+total+"\tsites," +
//				" or "+String.format("%.4f", (removedSites*100f/total))+"%");
//		System.err.println("Banned\t"+removedKeys+"\t/ "+data.length+"\tkeys," +
//				" or  "+String.format("%.4f", (removedKeys*100f/data.length))+"%");
		
		return set;
		
	}
	
	public static final boolean isClumpy(final int[] array, final int maxDist, final float fraction){
		if(array==null){return false;}
		int count=0;
		for(int i=1; i<array.length; i++){
			int dif=array[i]-array[i-1];
			if(dif<=maxDist){count++;}
		}
		return count>=(array.length*fraction);
	}

	public static int[] makeLengthHistogram(int[][] x, int buckets) {
		int[] lengths=new int[x.length];
		long total=0;
		for(int i=0; i<x.length; i++){
			int[] list=x[i];
			if(list!=null){
				lengths[i]=list.length;
				total+=list.length;
			}
		}
		Arrays.sort(lengths);
		
		int[] hist=new int[buckets+1];
		
		long sum=0;
		int ptr=0;
		for(int i=0; i<buckets; i++){
			long nextLimit=((total*i)+buckets/2)/buckets;
			while(ptr<lengths.length && sum<nextLimit){
				sum+=lengths[ptr];
				ptr++;
			}
			
			hist[i]=lengths[Tools.max(0, ptr-1)];
		}
		hist[hist.length-1]=lengths[lengths.length-1];
		
//		System.out.println(Arrays.toString(hist));
//		assert(false);
		return hist;
	}
	
	public static String toKMG(long x){
		double div=1;
		String ext="";
		if(x>10000000000000L){
			div=1000000000000L;
			ext="T";
		}else if(x>10000000000L){
			div=1000000000L;
			ext="B";
		}else if(x>10000000){
			div=1000000;
			ext="M";
		}else if(x>100000){
			div=1000;
			ext="K";
		}
		return String.format("%.2f", x/div)+ext;
	}
	
	public static long parseKMG(String b){

		char c=Character.toLowerCase(b.charAt(b.length()-1));
		if(!Character.isLetter(c) && !b.contains(".")){
			return Long.parseLong(b);
		}
		
		long mult=1;
		if(Character.isLetter(c)){
			if(c=='k'){mult=1000;}
			else if(c=='m'){mult=1000000;}
			else if(c=='g' || c=='b'){mult=1000000000;}
			else if(c=='t'){mult=1000000000000L;}
			else{throw new RuntimeException(b);}
			b=b.substring(0, b.length()-1);
		}
		
		return ((long)Double.parseDouble(b))*mult;
		
	}
	
	public static boolean parseBoolean(String s){
		if(s==null || s.length()<1){return true;}
		if(s.length()==1){
			char c=Character.toLowerCase(s.charAt(0));
			return c=='t' || c=='1';
		}
		if(s.equalsIgnoreCase("null") || s.equalsIgnoreCase("none")){return false;}
		return Boolean.parseBoolean(s);
	}
	
	public static int parseInt(byte[] array, int a, int b){
		assert(b>a);
		int r=0;
		final byte z='0';
		boolean negative=false;
		if(array[a]=='-'){negative=true; a++;}
		for(; a<b; a++){
			int x=(array[a]-z);
			assert(x<10 && x>=0) : x+" = "+(char)array[a]+"\narray="+new String(array)+", start="+a+", stop="+b;
			r=(r*10)+x;
		}
		if(negative){r*=-1;}
		return r;
	}
	
	/** TODO:  This (temporarily) uses a lot of memory.  Could be reduced by making an array of length max(x) and counting occurrences. */
	public static int[] makeLengthHistogram2(int[] x, int buckets, boolean verbose) {
		int[] lengths=Arrays.copyOf(x, x.length);
		long total=sum(x);
		Arrays.sort(lengths);
		
		if(verbose){
			System.out.println("Length array size:\t"+x.length);
			System.out.println("Min value:        \t"+lengths[0]);
			System.out.println("Med value:        \t"+lengths[lengths.length/2]);
			System.out.println("Max value:        \t"+lengths[lengths.length-1]);
			System.out.println("Total:            \t"+total);
		}
		
		int[] hist=new int[buckets+1];
		
		long sum=0;
		int ptr=0;
		for(int i=0; i<buckets; i++){
			long nextLimit=((total*i)+buckets/2)/buckets;
			while(ptr<lengths.length && sum<nextLimit){
				sum+=lengths[ptr];
				ptr++;
			}
			
			hist[i]=lengths[Tools.max(0, ptr-1)];
		}
		hist[hist.length-1]=lengths[lengths.length-1];
		
//		System.out.println(Arrays.toString(hist));
//		assert(false);
		return hist;
	}
	
	public static int[] makeLengthHistogram3(int[] x, int buckets, boolean verbose) {
		int max=max(x);
		if(max>x.length){
			Data.sysout.println("Reverted to old histogram mode.");
			return makeLengthHistogram2(x, buckets, verbose);
		}
		
		int[] counts=new int[max+1];
		long total=0;
		for(int i=0; i<x.length; i++){
			int a=x[i];
			if(a>=0){
				counts[a]++;
				total+=a;
			}
		}
		
		return makeLengthHistogram4(counts, buckets, total, verbose);
	}
	
	/** Uses counts of occurrences of lengths rather than raw lengths */
	public static int[] makeLengthHistogram4(int[] counts, int buckets, long total, boolean verbose) {
		if(total<=0){
			total=0;
			for(int i=1; i<counts.length; i++){
				total+=(i*counts[i]);
			}
		}
		
		if(verbose){
//			System.out.println("Length array size:\t"+x.length);
//			System.out.println("Min value:        \t"+lengths[0]);
//			System.out.println("Med value:        \t"+lengths[lengths.length/2]);
//			System.out.println("Max value:        \t"+lengths[lengths.length-1]);
			System.out.println("Total:            \t"+total);
		}
		
		int[] hist=new int[buckets+1];
		
		long sum=0;
		int ptr=0;
		for(int i=0; i<buckets; i++){
			long nextLimit=((total*i)+buckets/2)/buckets;
			while(ptr<counts.length && sum<nextLimit){
				sum+=counts[ptr]*ptr;
				ptr++;
			}
			
			hist[i]=Tools.max(0, ptr-1);
		}
		hist[hist.length-1]=counts.length-1;
		
//		System.out.println(Arrays.toString(hist));
//		assert(false);
		return hist;
	}
	
	/**
	 * @param cov
	 * @return
	 */
	public static int average(short[] array) {
		// TODO Auto-generated method stub
		return (int)(array==null || array.length==0 ? 0 : sum(array)/array.length);
	}
	
	/**
	 * @param cov
	 * @return
	 */
	public static int average(int[] array) {
		// TODO Auto-generated method stub
		return (int)(array==null || array.length==0 ? 0 : sum(array)/array.length);
	}
	
	public static int median(int[] array){return percentile(array, .5);}
	
	public static int percentile(int[] array, double fraction){
		if(array==null || array.length<1){return 0;}
		long target=(long)(sum(array)*fraction);
		long sum=0;
		for(int i=0; i<array.length; i++){
			sum+=array[i];
			if(sum>=target){
				return i;
			}
		}
		return array.length-1;
	}

	public static int absdif(int a, int b) {
		return a>b ? a-b : b-a;
	}

	public static float absdif(float a, float b) {
		return a>b ? a-b : b-a;
	}

	public static double absdif(double a, double b) {
		return a>b ? a-b : b-a;
	}
	
	/** Uses unsigned math */
	public static final int absdifUnsigned(int a, int b){
		return (a<0 == b<0) ? a>b ? a-b : b-a : Integer.MAX_VALUE;
	}
	
	public static final boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1; 
	}
	
	public static final int overlapLength(int a1, int b1, int a2, int b2){
		if(!overlap(a1,b1,a2,b2)){return 0;}
		if(a1<=a2){
			return b1>=b2 ? b2-a2+1 : b1-a2+1;
		}else{
			return b2>=b1 ? b1-a1+1 : b2-a1+1;
		}
	}
	
	/** Is (a1, b1) within (a2, b2) ? */
	public static final boolean isWithin(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a1>=a2 && b1<=b2; 
	}
	
	public static final int constrict(int point, int a, int b){
		assert(a<=b);
		return(point<a ? a : point>b ? b : point);
	}
	
	public static final int indexOf(byte[] array, byte b){
		int i=0;
		while(i<array.length && array[i]!=b){i++;}
		return (i==array.length ? -1 : i);
	}
	
	public static final int indexOf(char[] array, char b){
		int i=0;
		while(i<array.length && array[i]!=b){i++;}
		return (i==array.length ? -1 : i);
	}
	
	public static final int lastIndexOf(byte[] array, byte b){
		int i=array.length-1;
		while(i>=0 && array[i]!=b){i--;}
		return i;
	}
	
	public static final int stringLength(long x){
		if(x<0){
			if(x==Integer.MIN_VALUE){return 11;}
			return lengthOf(-x)+1;
		}
		return lengthOf(x);
	}
	
	public static final int stringLength(int x){
		if(x<0){
			if(x==Long.MIN_VALUE){return 20;}
			return lengthOf(-x)+1;
		}
		return lengthOf(x);
	}
	
	public static final int lengthOf(int x){
		assert(x>=0);
		int i=1;
		while(x>ilens[i]){i++;}
		return i;
	}
	
	public static final int lengthOf(long x){
		assert(x>=0);
		int i=1;
		while(x>llens[i]){i++;}
		return i;
	}
	
	public static final int max(int[] array){return array[maxIndex(array)];}
	
	public static final int maxIndex(int[] array){
		int max=array[0], maxIndex=0;
		for(int i=1; i<array.length; i++){
			if(array[i]>max){max=array[i];maxIndex=i;}
		}
		return maxIndex;
	}
	
	public static final double standardDeviation(long[] numbers){
		if(numbers==null || numbers.length<1){return 0;}
		long sum=sum(numbers);
		double avg=sum/(double)numbers.length;
		double sumdev2=0;
		for(int i=0; i<numbers.length; i++){
			long x=numbers[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/numbers.length);
	}
	
	public static final double standardDeviation(int[] numbers){
		if(numbers==null || numbers.length<1){return 0;}
		long sum=sum(numbers);
		double avg=sum/(double)numbers.length;
		double sumdev2=0;
		for(int i=0; i<numbers.length; i++){
			long x=numbers[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/numbers.length);
	}
	
	public static final double standardDeviation(short[] numbers){
		if(numbers==null || numbers.length<1){return 0;}
		long sum=sum(numbers);
		double avg=sum/(double)numbers.length;
		double sumdev2=0;
		for(int i=0; i<numbers.length; i++){
			long x=numbers[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/numbers.length);
	}
	
	public static final double standardDeviationHistogram(long[] histogram){
		long sum=max(1, sum(histogram));
		long sum2=0;
		for(int i=0; i<histogram.length; i++){
			sum2+=(histogram[i]*i);
		}
		double avg=sum2/(double)sum;
		double sumdev2=0;
		for(int i=0; i<histogram.length; i++){
			double dev=avg-i;
			double dev2=dev*dev;
			sumdev2+=(histogram[i]*dev2);
		}
		return Math.sqrt(sumdev2/sum);
	}
	
	/** Special version that calculates standard deviation based on unique kmers rather than overall events */
	public static final double standardDeviationHistogramKmer(long[] histogram){
		final long sum=sum(histogram);
		double sumU=0;
		for(int i=0; i<histogram.length; i++){
			long x=histogram[i];
			sumU+=(x/(double)max(i, 1));
		}
		double avg=sum/max(sumU, 1);
		double sumdev2=0;
		for(int i=1; i<histogram.length; i++){
			double dev=avg-i;
			double dev2=dev*dev;
			long x=histogram[i];
			sumdev2+=((x/(double)max(i, 1))*dev2);
		}
		return Math.sqrt(sumdev2/sumU);
	}
	
	public static final double standardDeviationHistogram(AtomicLongArray histogram){
		long sum=max(1, sum(histogram));
		long sum2=0;
		for(int i=0; i<histogram.length(); i++){
			sum2+=(histogram.get(i)*i);
		}
		double avg=sum2/(double)sum;
		double sumdev2=0;
		for(int i=0; i<histogram.length(); i++){
			double dev=avg-i;
			double dev2=dev*dev;
			sumdev2+=(histogram.get(i)*dev2);
		}
		return Math.sqrt(sumdev2/sum);
	}
	
	/** Special version that calculates standard deviation based on unique kmers rather than overall events */
	public static final double standardDeviationHistogramKmer(AtomicLongArray histogram){
		final long sum=sum(histogram);
		double sumU=0;
		for(int i=0; i<histogram.length(); i++){
			long x=histogram.get(i);
			sumU+=(x/(double)max(i, 1));
		}
		double avg=sum/max(sumU, 1);
		double sumdev2=0;
		for(int i=1; i<histogram.length(); i++){
			double dev=avg-i;
			double dev2=dev*dev;
			long x=histogram.get(i);
			sumdev2+=((x/(double)max(i, 1))*dev2);
		}
		return Math.sqrt(sumdev2/sumU);
	}
	
	public static final long[] downsample(long[] array, int bins){
		if(array==null || array.length==bins){return array;}
		assert(bins<=array.length);
		assert(bins>=0);
		long[] r=new long[bins];
		if(bins==0){return r;}
		double mult=bins/(double)array.length;
		for(int i=0; i<array.length; i++){
			int j=(int)(mult*i);
			r[j]+=array[i];
		}
		return r;
	}

	
	public static final void pause(int millis){
		try {
			Thread.sleep(millis);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}
	public static final int min(int x, int y, int z){return x<y ? (x<z ? x : z) : (y<z ? y : z);}
	public static final int max(int x, int y, int z){return x>y ? (x>z ? x : z) : (y>z ? y : z);}
	public static final int min(int x, int y, int z, int z2){return min(min(x,y), min(z,z2));}
	public static final int max(int x, int y, int z, int z2){return max(max(x,y), max(z,z2));}
	
	//Median of 3
	public static final int mid(int x, int y, int z){return x<y ? (x<z ? min(y, z) : x) : (y<z ? min(x, z) : y);}

	public static final byte min(byte x, byte y){return x<y ? x : y;}
	public static final byte max(byte x, byte y){return x>y ? x : y;}

	public static final char min(char x, char y){return x<y ? x : y;}
	public static final char max(char x, char y){return x>y ? x : y;}

	public static final byte min(byte x, byte y, byte z){return x<y ? min(x, z) : min(y, z);}
	public static final byte max(byte x, byte y, byte z){return x>y ? max(x, z) : max(y, z);}

	public static final byte min(byte x, byte y, byte z, byte a){return min(min(x, y), min(z, a));}
	public static final byte max(byte x, byte y, byte z, byte a){return max(max(x, y), max(z, a));}
	
	public static final long min(long x, long y){return x<y ? x : y;}
	public static final long max(long x, long y){return x>y ? x : y;}
	
	public static final long min(long x, long y, long z){return x<y ? (x<z ? x : z) : (y<z ? y : z);}
	public static final long max(long x, long y, long z){return x>y ? (x>z ? x : z) : (y>z ? y : z);}
	
	public static final double min(double x, double y){return x<y ? x : y;}
	public static final double max(double x, double y){return x>y ? x : y;}
	
	public static final float min(float x, float y){return x<y ? x : y;}
	public static final float max(float x, float y){return x>y ? x : y;}
	
	public static final int min(int[] array, int fromIndex, int toIndex){
		int min=array[fromIndex];
		for(int i=fromIndex+1; i<=toIndex; i++){
			min=min(min, array[i]);
		}
		return min;
	}
	
	public static final int max(int[] array, int fromIndex, int toIndex){
		int max=array[fromIndex];
		for(int i=fromIndex+1; i<=toIndex; i++){
			max=max(max, array[i]);
		}
		return max;
	}

	public static int minIndex(int[] array) {
		if(array==null || array.length<1){return -1;}
		float min=array[0];
		int index=0;
		for(int i=1; i<array.length; i++){
			if(array[i]<min){
				min=array[i];
				index=i;
			}
		}
		return index;
	}
	
	public static double log2(double d){
		return Math.log(d)*invlog2;
	}
	
	public static double logRoot2(double d){
		return Math.log(d)*invlogRoot2;
	}
	
	public static double log1point2(double d){
		return Math.log(d)*invlog1point2;
	}

	private static final double log2=Math.log(2);
	private static final double invlog2=1/log2;
	private static final double logRoot2=Math.log(Math.sqrt(2));
	private static final double invlogRoot2=1/logRoot2;
	private static final double log1point2=Math.log(1.2);
	private static final double invlog1point2=1/log1point2;
	
	public static final int[] ilens;
	public static final long[] llens;
	
	static{
		ilens=new int[Integer.toString(Integer.MAX_VALUE).length()+1];
		llens=new long[Long.toString(Long.MAX_VALUE).length()+1];
		for(int i=1, x=9; i<ilens.length; i++){
			ilens[i]=x;
			x=(x*10)+9;
		}
		ilens[ilens.length-1]=Integer.MAX_VALUE;
		for(long i=1, x=9; i<llens.length; i++){
			llens[(int)i]=x;
			x=(x*10)+9;
		}
		llens[llens.length-1]=Long.MAX_VALUE;
	}
}
