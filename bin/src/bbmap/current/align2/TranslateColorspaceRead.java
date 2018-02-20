package align2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import stream.Read;
import stream.SiteScore;
import var.Variation;
import var.Varlet;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;
import dna.Timer;

public final class TranslateColorspaceRead {
	
	
	public static void main(String[] args){
		
//		byte[] basesRef="CAGTAGTC".getBytes();
//		byte[] colorsRef=AminoAcid.toColorspace(basesRef);
//		
//		byte[] basesCall="CAGAGTC".getBytes();
//		byte[] colors=AminoAcid.toColorspace(basesCall);
//		
//		byte[] match="mmDSmmm".getBytes();
		
//		byte[] basesRef="CAGAGTC".getBytes();
//		byte[] colorsRef=AminoAcid.toColorspace(basesRef);
//		
//		byte[] basesCall="CAGTAGTC".getBytes();
//		byte[] colors=AminoAcid.toColorspace(basesCall);
//		
//		byte[] match="mmISmmm".getBytes();
//		
//		byte[][] crbmq=new byte[][] {colors, colorsRef, basesRef, match};
//		
//		System.err.println(toString(crbmq));
//		fixIndels(crbmq);
//		System.err.println();
//		System.err.println(toString(crbmq));
		
		int rounds=100000;
		
		int correct=test2(rounds);
		String missed="Missed "+String.format("%.2f", ((rounds-correct)*100f/rounds))+"%";
//		System.err.println(missed);
		System.err.println(missed);
	}
	
	public TranslateColorspaceRead(MSA msa){
		msaBS=msa;
	}
	
	public static int test2(int rounds){
		int length=50;
		
		Timer t=new Timer();
		
		int correct=0;
		
		MSA temp=new MultiStateAligner9ts(500, 2400, false);
		TranslateColorspaceRead tcr=new TranslateColorspaceRead(temp);
		
		RandomReads rr=new RandomReads(false);
		
		int maxSnps=2, maxInss=2, maxDels=2, maxSubs=2, maxErrors=2;
		float snpRate=.3f, insRate=.3f, delRate=.3f, subRate=.3f;
		int maxIndelLen=5;
		int maxSubLen=5;
		
		Read[] reads=rr.makeRandomReadsX(rounds, length, maxSnps, maxInss, maxDels, maxSubs,
				snpRate, insRate, delRate, subRate, 
				maxIndelLen, maxIndelLen, maxSubLen, (byte)21, (byte)21, true, 5, 20, 30);
		
		t.start();
		for(int i=0; i<reads.length; i++){
			Read r=reads[i];

			if(r.strand()==Gene.MINUS){
				Tools.reverseInPlace(r.bases);
			}	
			rr.addColorspaceErrors(r, maxErrors);
			
//			int[] max=msa.fill(r.bases, chacs.array, r.start, r.stop);
//			r.match=msa.traceback(r.bases, chacs.array, r.start, r.stop, max[0], max[1], max[2]);
			
			if(r.strand()==Gene.MINUS){
				Tools.reverseInPlace(r.bases);
			}
			
			Read r2=tcr.translateToBasespace(r);
			boolean ok=r2==null ? false : r.obj.equals(new String(r2.bases));
			if(verbose){System.err.println(r.obj+"\t(original)\t"+(ok ? "" : "FAILED"));}
			if(ok){correct++;}
			
			if(ok){
				ArrayList<Varlet> list=tcr.toVars(r2, true, true, false);
				if(verbose && list!=null){
					for(Varlet v : list){
						System.err.println(v);
					}
				}
			}
			
		}
		t.stop();
		float kbps=(length*(long)rounds)*1000000f/t.elapsed;
		System.err.println("Time: "+t+"\t("+String.format("%.2f", kbps)+" kbps)");
		return correct;
	}
	
	private static CharSequence toString(byte[][] crbmq) {
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<2; i++){
			if(crbmq[i]==null){sb.append("null");}
			else{
				for(byte b : crbmq[i]){
					if(b=='N'){sb.append('N');}
					else{sb.append((char)(b+'0'));}
				}
				sb.append('\n');
			}
		}
		sb.append(new String(crbmq[2]));
		sb.append('\n');
		sb.append(crbmq[3]==null ? "null" : new String(crbmq[3]));
		sb.append('\n');
		return sb;
	}
	
	private static String toStringCS(byte[] colors){
		StringBuilder sb=new StringBuilder(colors.length);
		for(byte b : colors){
			if(b>3){sb.append((char)b);}
			else{sb.append((char)(b+'0'));}
		}
		sb.append('\n');
		return sb.toString();
	}
	
	public void realign_new(final Read r, final int padding, final boolean recur, final int minValidScore, boolean forbidIndels){
		realign_new(r, r.colorspace() ? msaCS : msaBS, padding, recur, minValidScore, forbidIndels);
	}
	
	public void realignByReversingRef(final Read r, final int padding, final boolean recur){
		realignByReversingRef(r, r.colorspace() ? msaCS : msaBS, padding, recur);
	}
	
	/** This aligns a read with the reference, and generates the match string. */
	public static void realignByReversingRef(final Read r, final MSA msa, int padding, final boolean recur){
		if(r.shortmatch()){
			r.match=null;
			r.setShortMatch(false);
		}
//		assert(r.colorspace());
//		assert(msa.colorspace);
		padding=Tools.min(padding, (msa.maxColumns-r.bases.length)/2-20);
		padding=Tools.max(padding, 0);
		assert(r.colorspace()==msa.colorspace);
		final ChromosomeArray chacs=Data.getChromosome(r.chrom);
		if(verbose){
			System.err.println("Realigning.");
			System.err.println("Original: "+r.start+", "+r.stop+", "+Gene.strandCodes[r.strand()]);
		}

		{
			assert(r.stop>=r.start); //Otherwise this is pointless...
			int a=r.bases.length;
			int b=r.stop-r.start+1;
			if(b<a){
				int c=Tools.min(r.bases.length, a-b+10)/2;
				padding=Tools.max(padding, c+1);
			}
		}
		padding=Tools.min(padding, r.bases.length+10);
		padding=Tools.min(padding, (msa.maxColumns-Tools.max(r.bases.length, GapTools.calcGrefLen(r.start, r.stop, r.gaps)))/2-1);
		
//		if(padding==4){System.err.print(".");}
//		else{
//			System.err.println("\npadding="+padding+", \trecur="+recur);
//			if(padding>10){
//				if(r.match!=null){System.err.print(new String(r.match));}
//				if(padding>20){System.err.print("\t"+r.start+", "+r.stop+", "+(r.stop-r.start+1));}
//				System.err.println();
//			}
//		}

		final int maxQ=msa.maxQuality(r.bases.length);
		final int maxI=msa.maxImperfectScore(r.bases.length);

		if(r.strand()==Gene.PLUS){
			assert(maxQ>maxI);

			byte[][] matchR=new byte[1][];
			if(r.match!=null && r.match.length==r.bases.length){
				matchR[0]=r.match;
			}else{
				//				System.err.println(new String(r.match));
				matchR[0]=r.match=new byte[r.bases.length];
			}
			int scoreNoIndel=msa.scoreNoIndelsAndMakeMatchString(r.bases, chacs.array, r.start, matchR);
			r.match=matchR[0];
			
			if(scoreNoIndel>=maxI){
				if(verbose){System.err.println("Quick match.");}
//				assert(r.match[0]!='X') : r.toText(false);
//				assert(r.match[r.match.length-1]!='X') : r.toText(false);
				//				assert(r.stop==r.start+r.bases.length-1);
				r.stop=r.start+r.bases.length-1;
				r.mapScore=scoreNoIndel;
			}else{
				if(verbose){System.err.println("Slow match.");}
				
//				int minLoc=Tools.max(r.start-padding, chacs.minIndex);
				int minLoc=Tools.max(r.start-padding, 0); //It's OK to be off the beginning as long as bases prior to the true start are 'N'
				int maxLoc=Tools.min(r.stop+padding, chacs.maxIndex);

				//These assertions are not too important... they indicate the read mapped off the end of the chromosome.
				assert(minLoc<=r.start) : "\nchr"+r.chrom+": "+minLoc+", "+maxLoc+", "+r.start+", "+r.stop+
					", "+chacs.minIndex+", "+chacs.maxIndex+"\n"+r.toText(false);
				assert(maxLoc>=r.stop) : "\nchr"+r.chrom+": "+minLoc+", "+maxLoc+", "+r.start+", "+r.stop+
					", "+chacs.minIndex+", "+chacs.maxIndex+"\n"+r.toText(false);

				//			System.err.println("Aligning:\n"+new String(r.bases)+"\n"+chacs.getString(minLoc, maxLoc));
				int[] max=msa.fillLimited(r.bases, chacs.array, minLoc, maxLoc, scoreNoIndel, r.gaps);
				//			System.err.println(Arrays.toString(max));
				r.match=msa.traceback(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null);
//				System.err.println(new String(r.match));
				int[] score=msa.score(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null);
//				System.err.println(Arrays.toString(score));
				r.start=score[1];
				r.stop=score[2];
				r.mapScore=score[0];
				//			System.err.println(Arrays.toString(score));
				//			assert(false);
			}
		}else{
			assert(maxQ>maxI);

			byte[][] matchR=new byte[1][];
			if(r.match!=null && r.match.length==r.bases.length){
				matchR[0]=r.match;
			}else{
				//				System.err.println(new String(r.match));
				matchR[0]=r.match=new byte[r.bases.length];
			}
			
			int scoreNoIndel=-9999;
			if(r.bases.length==(r.stop-r.start+1)){
				
				byte[] ref=chacs.getBytes(r.start, r.stop);
				if(r.colorspace()){
					Tools.reverseInPlace(ref);
				}else{
					AminoAcid.reverseComplementBasesInPlace(ref);
				}
				scoreNoIndel=msa.scoreNoIndelsAndMakeMatchString(r.bases, ref, 0, matchR);
				r.match=matchR[0];
			}

			if(scoreNoIndel>=maxI){
				if(verbose){System.err.println("Quick match.");}
				assert(r.match[0]!='X') : r.toText(false);
				assert(r.match[r.match.length-1]!='X') : r.toText(false);
				r.stop=r.start+r.bases.length-1;
				r.mapScore=scoreNoIndel;
			}else{
				if(verbose){System.err.println("Slow match.");}

//				int minLoc=Tools.max(r.start-padding, chacs.minIndex);
				int minLoc=Tools.max(r.start-padding, 0); //It's OK to be off the beginning as long as bases prior to the true start are 'N'
				int maxLoc=Tools.min(r.stop+padding, chacs.maxIndex);

				//These assertions are not too important... they indicate the read mapped off the end of the chromosome.
				assert(minLoc<=r.start) : "\nchr"+r.chrom+": "+minLoc+", "+maxLoc+", "+r.start+", "+r.stop+
					", "+chacs.minIndex+", "+chacs.maxIndex+"\n"+r.toText(false);
				assert(maxLoc>=r.stop) : "\nchr"+r.chrom+": "+minLoc+", "+maxLoc+", "+r.start+", "+r.stop+
					", "+chacs.minIndex+", "+chacs.maxIndex+"\n"+r.toText(false);

				byte[] ref=chacs.getBytes(minLoc, maxLoc);
				//			System.err.println("Aligning:\n"+new String(r.bases)+"\n"+new String(ref));
				if(r.colorspace()){
					Tools.reverseInPlace(ref);
				}else{
					AminoAcid.reverseComplementBasesInPlace(ref);
				}

				//			System.err.println("Aligning:\n"+new String(r.bases)+"\n"+new String(ref));
				int[] max=msa.fillLimited(r.bases, ref, 0, ref.length-1, scoreNoIndel, r.gaps);
				//			System.err.println(Arrays.toString(max));
				r.match=msa.traceback(r.bases, ref, 0, ref.length-1, max[0], max[1], max[2], r.gaps!=null);
//				System.err.println(new String(r.match));
				int[] score=msa.score(r.bases, ref, 0, ref.length-1, max[0], max[1], max[2], r.gaps!=null);
//				System.err.println(Arrays.toString(score));
				//			System.err.println(Arrays.toString(score));
				//			assert(false);

				int start2=minLoc+(ref.length-score[2]-1);
				int stop2=maxLoc-(score[1]);

				r.start=start2;
				r.stop=stop2;
				r.mapScore=score[0];
			}
		}
		if(verbose){System.err.println("Final: "+r.start+", "+r.stop+", "+Gene.strandCodes[r.strand()]);}
		
		if(recur && r.stop<chacs.maxIndex && r.start>0 && (r.match[0]=='X' || r.match[0]=='I' || 
				r.match[r.match.length-1]=='Y' || r.match[r.match.length-1]=='X' || r.match[r.match.length-1]=='I')){
			int xy=0;
			for(int i=0; i<r.match.length; i++){
				byte b=r.match[i];
				if(b=='X' || b=='Y' || b=='I'){xy++;}
			}
//			System.err.println("xy = "+xy);
			realignByReversingRef(r, msa, Tools.min(10+padding+2*xy, msa.maxColumns/2-r.bases.length-20), false);
		}
//		assert(r.mapScore>0) : padding+", "+recur+", "+r.mapScore+", "+r.strand()+", "+r.colorspace()+"\n"+r.toText(false);
//		assert(r.match[0]!='X') : r.toText(false);
//		assert(r.match[r.match.length-1]!='X') : r.toText(false);
	}
	
	/** For some reason realign was making the match string backwards... */
	public static void realign_new(final Read r, final MSA msa, int padding, final boolean recur, int minValidScore, boolean forbidIndels){
		if(r.shortmatch()){
			r.match=null;
			r.setShortMatch(false);
		}
		if(verbose){System.err.println("Padding = "+padding+"; msa.maxColumns = "+msa.maxColumns+"; maplen = "+(r.stop-r.start+1)+"; gaps = "+Arrays.toString(r.gaps));}
		
		assert(padding>=0) : padding+", "+r;
		padding=Tools.min(padding, (msa.maxColumns-r.bases.length)/2-20);
		if(verbose){System.err.println("Padding = "+padding);}
		assert(padding>=0) : padding+", "+r;
		padding=Tools.max(padding, 0);
		if(verbose){System.err.println("Padding = "+padding);}
		assert(padding>=0) : padding+", "+r;

		
		assert(r.colorspace()==msa.colorspace);
		final ChromosomeArray chacs=Data.getChromosome(r.chrom);
		if(verbose){
			System.err.println("Realigning.");
			System.err.println("Original: "+r.start+", "+r.stop+", "+Gene.strandCodes[r.strand()]);
			if(verbose){System.err.println("Estimated greflen: "+GapTools.calcGrefLen(r.start, r.stop, r.gaps));}
		}
		
		{
			int expectedLen=GapTools.calcGrefLen(r.start, r.stop, r.gaps);
			if(expectedLen>msa.maxColumns-20){
				//TODO: Alternately, I could kill the site.
				r.stop=r.start+Tools.min(r.bases.length+40, msa.maxColumns-20);
				if(r.gaps!=null){r.gaps=GapTools.fixGaps(r.start, r.stop, r.gaps, Shared.MINGAP);}
			}
		}
		
		if(r.start<0){r.start=0;} //Prevents assertion errors.  This change should be reset by the realignment so it shouldn't matter.
		if(r.stop>chacs.maxIndex){r.stop=chacs.maxIndex;} //Also to prevent a potential assertion error in unpadded references
		assert(0<=r.start) : "\nchr"+r.chrom+": r.start="+r.start+", r.stop="+r.stop+", padding="+padding+
			", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+r.toText(false);
		assert(chacs.maxIndex>=r.stop) : "\nchr"+r.chrom+": r.start="+r.start+", r.stop="+r.stop+", padding="+padding+
			", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+r.toText(false);

		{
			assert(r.stop>=r.start); //Otherwise this is pointless...
			int a=r.bases.length;
			int b=r.stop-r.start+1;
			if(b<a){
				int c=Tools.min(r.bases.length, a-b+10)/2;
				padding=Tools.max(padding, c+1);
//				if(verbose){System.err.println("Padding = "+padding);}
//				assert(padding>=0) : padding;
			}
		}
		
		if(verbose){System.err.println("Padding = "+padding);}
		
		padding=Tools.min(padding, (msa.maxColumns-Tools.max(r.bases.length, GapTools.calcGrefLen(r.start, r.stop, r.gaps)))/2-1);
		if(forbidIndels){padding=0;}
		if(verbose){System.err.println("Padding = "+padding);}
		assert(padding>=0) : r.numericID+", "+padding;
		
		
		final int maxQ=msa.maxQuality(r.bases.length);
		final int maxI=msa.maxImperfectScore(r.bases.length);

		if(r.strand()==Gene.PLUS){
			assert(maxQ>maxI);

			byte[][] matchR=new byte[1][];
			if(r.match!=null && r.match.length==r.bases.length){
				matchR[0]=r.match;
			}else{
				//				System.err.println(new String(r.match));
				matchR[0]=r.match=new byte[r.bases.length];
			}
			int scoreNoIndel=msa.scoreNoIndelsAndMakeMatchString(r.bases, chacs.array, r.start, matchR);
			r.match=matchR[0];
			
			assert(0<=r.start) : "\nchr"+r.chrom+": r.start="+r.start+", r.stop="+r.stop+", padding="+padding+
				", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+r.toText(false);
			assert(chacs.maxIndex>=r.stop) : "\nchr"+r.chrom+": r.start="+r.start+", r.stop="+r.stop+", padding="+padding+
				", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+r.toText(false);
			
			if(verbose){System.err.println("Estimated greflen: "+GapTools.calcGrefLen(r.start, r.stop, r.gaps));}
			
			if(scoreNoIndel>=maxI || forbidIndels){
				if(verbose){System.err.println("Quick match.");}
//				assert(r.match[0]!='X') : r.toText(false);
//				assert(r.match[r.match.length-1]!='X') : r.toText(false);
				//				assert(r.stop==r.start+r.bases.length-1);
				r.stop=r.start+r.bases.length-1;
				r.mapScore=scoreNoIndel;
			}else{
				if(verbose){System.err.println("Slow match.");}
				
//				int minLoc=Tools.max(r.start-padding, chacs.minIndex);
				int minLoc=Tools.max(r.start-padding, 0); //It's OK to be off the beginning as long as bases prior to the true start are 'N'
				int maxLoc=Tools.min(r.stop+padding, chacs.maxIndex);
				
				if(verbose){System.err.println("minLoc = "+minLoc+", maxLoc = "+maxLoc);}
				if(verbose){System.err.println("A. Estimated greflen: "+GapTools.calcGrefLen(r.start, r.stop, r.gaps));}
				if(verbose){System.err.println("A. Estimated greflen2: "+GapTools.calcGrefLen(minLoc, maxLoc, r.gaps));}
				
				//These assertions are not too important... they indicate the read mapped off the end of the chromosome.
				assert(minLoc<=r.start) : "\nchr"+r.chrom+": minloc="+minLoc+", maxLoc="+maxLoc+", r.start="+r.start+", r.stop="+r.stop+", padding="+padding+
					", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+r.toText(false);
				assert(maxLoc>=r.stop) : "\nchr"+r.chrom+": minloc="+minLoc+", maxLoc="+maxLoc+", r.start="+r.start+", r.stop="+r.stop+", padding="+padding+
					", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+r.toText(false);

				//			System.err.println("Aligning:\n"+new String(r.bases)+"\n"+chacs.getString(minLoc, maxLoc));
				
				int[] max=null;
				int[] score=null;
				try {
					if(verbose){
						System.err.println("Calling fillLimited(bases, chacs, "+minLoc+", "+maxLoc+", "+
								Tools.max(scoreNoIndel, minValidScore)+", "+(r.gaps==null ? "null" : Arrays.toString(r.gaps))+")");
					}
					max=msa.fillLimited(r.bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), r.gaps);
					score=(max==null ? null : msa.score(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null));
					if(verbose){System.err.println("Estimated greflen: "+GapTools.calcGrefLen(r.start, r.stop, r.gaps));}
					
					if(score!=null && score.length>6){
						int[] oldArray=score.clone();
						assert(score.length==8);
						int extraPadLeft=score[6];
						int extraPadRight=score[7];
						
						if(r.gaps==null){
							assert(maxLoc-minLoc+1<=msa.maxColumns);
							int newlen=(maxLoc-minLoc+1+extraPadLeft+extraPadRight);
							if(newlen>=msa.maxColumns){
								while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
								while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
								while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
							}else{
								int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
								extraPadLeft=Tools.max(x, extraPadLeft);
								extraPadRight=Tools.max(x, extraPadRight);
							}
						}else{
							//TODO: In this case the alignment will probably be wrong.
							int greflen=Tools.max(r.bases.length, GapTools.calcGrefLen(minLoc, maxLoc, r.gaps));
							int newlen=(greflen+1+extraPadLeft+extraPadRight);
							if(newlen>=msa.maxColumns){
								while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
								while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
								while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
							}else{
								int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
								extraPadLeft=Tools.max(x, extraPadLeft);
								extraPadRight=Tools.max(x, extraPadRight);
							}
						}
						
						assert(extraPadLeft>=0 && extraPadRight>=0) : extraPadLeft+", "+extraPadRight+"\n"+r;
						minLoc=Tools.max(0, minLoc-extraPadLeft);
						maxLoc=Tools.min(chacs.maxIndex, maxLoc+extraPadRight);

						if(verbose){System.err.println("B. Estimated greflen: "+GapTools.calcGrefLen(r.start, r.stop, r.gaps));}
						if(verbose){System.err.println("B. Estimated greflen2: "+GapTools.calcGrefLen(minLoc, maxLoc, r.gaps));}
						max=msa.fillLimited(r.bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), r.gaps);
						score=(max==null ? null : msa.score(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null));
						
						if(score==null || score[0]<oldArray[0]){
							if(!Shared.anomaly){System.err.println("Read "+r.numericID+": Padded match string alignment result was inferior.  Triple-aligning. :(");}
							
							if(r.gaps==null){
								assert(maxLoc-minLoc+1<=msa.maxColumns);
								int newlen=(maxLoc-minLoc+1+extraPadLeft+extraPadRight);
								if(newlen>=msa.maxColumns){
									while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
									while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
									while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
								}else{
									int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
									extraPadLeft=Tools.max(x, extraPadLeft);
									extraPadRight=Tools.max(x, extraPadRight);
								}
							}else{
								//TODO: In this case the alignment will probably be wrong.
								int greflen=Tools.max(r.bases.length, GapTools.calcGrefLen(minLoc, maxLoc, r.gaps));
								int newlen=(greflen+1+extraPadLeft+extraPadRight);
								if(newlen>=msa.maxColumns){
									while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
									while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
									while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
								}else{
									int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
									extraPadLeft=Tools.max(x, extraPadLeft);
									extraPadRight=Tools.max(x, extraPadRight);
								}
							}
							
							assert(extraPadLeft>=0 && extraPadRight>=0) : extraPadLeft+", "+extraPadRight+"\n"+r;
							minLoc=Tools.max(0, minLoc-extraPadLeft);
							maxLoc=Tools.min(chacs.maxIndex, maxLoc+extraPadRight);
							
							max=msa.fillLimited(r.bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), r.gaps);
							score=(max==null ? null : msa.score(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null));
							
							if(minLoc>0 && maxLoc<chacs.maxIndex && (score==null || score[0]<oldArray[0])){
								if(!Shared.anomaly){System.err.println("Still inferior.");}
								minLoc=Tools.max(r.start-8, 0); //It's OK to be off the beginning as long as bases prior to the true start are 'N'
								maxLoc=Tools.min(r.stop+8, chacs.maxIndex);
								max=msa.fillUnlimited(r.bases, chacs.array, minLoc, maxLoc, r.gaps);
								score=(max==null ? null : msa.score(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null));
							}
						}
					}
				} catch (Exception e) {
					System.err.println("Caught exception:\n");
					e.printStackTrace();
					assert(false) : r.toText(false);
				}
				
				if(max!=null){
					r.match=msa.traceback(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null);
					r.start=score[1];
					r.stop=score[2];
					r.mapScore=score[0];
				}else{
					r.stop=r.start+r.bases.length-1;
					r.mapScore=scoreNoIndel;
				}
			}
		}else{
			assert(maxQ>maxI);

			byte[][] matchR=new byte[1][];
			if(r.match!=null && r.match.length==r.bases.length){
				matchR[0]=r.match;
			}else{
				matchR[0]=r.match=new byte[r.bases.length];
			}
			
			if(verbose){
				System.err.println("Before reversed:");
				System.err.println(toStringCS(r.bases));
			}
			
			if(r.colorspace()){
				Tools.reverseInPlace(r.bases);
			}else{
				AminoAcid.reverseComplementBasesInPlace(r.bases);
			}
			
			if(verbose){
				System.err.println("Reversed.");
				System.err.println(toStringCS(r.bases));
			}
			
			int scoreNoIndel=msa.scoreNoIndelsAndMakeMatchString(r.bases, chacs.array, r.start, matchR);
			r.match=matchR[0];
			
			if(scoreNoIndel>=maxI || forbidIndels){
				if(verbose){System.err.println("Quick match.");}
				assert(r.match[0]!='X') : r.toText(false);
				assert(r.match[r.match.length-1]!='X') : r.toText(false);
				r.stop=r.start+r.bases.length-1;
				r.mapScore=scoreNoIndel;
//				Tools.reverseInPlace(r.match);
			}else{
				if(verbose){System.err.println("Slow match.");}
				
				int minLoc=Tools.max(r.start-padding, 0); //It's OK to be off the beginning as long as bases prior to the true start are 'N'
				int maxLoc=Tools.min(r.stop+padding, chacs.maxIndex);
				if(verbose){System.err.println("Slow match "+minLoc+" ~ "+maxLoc);}

				//These assertions are not too important... they indicate the read mapped off the end of the chromosome.
				assert(minLoc<=r.start) : "\nchr"+r.chrom+": "+minLoc+", "+maxLoc+", "+r.start+", "+r.stop+
					", "+chacs.minIndex+", "+chacs.maxIndex+"\n"+r.toText(false);
				assert(maxLoc>=r.stop) : "\nchr"+r.chrom+": "+minLoc+", "+maxLoc+", "+r.start+", "+r.stop+
					", "+chacs.minIndex+", "+chacs.maxIndex+"\n"+r.toText(false);

				if(verbose){System.err.println("Aligning:\n"+new String(r.bases)+"\n"+chacs.getString(minLoc, maxLoc));}
				int[] max=msa.fillLimited(r.bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), r.gaps);
				if(verbose){System.err.println("Aligned: {rows, maxC, maxS, max} = "+Arrays.toString(max));}
				int[] score=null;
				score=(max==null ? null : msa.score(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null));
				
				if(score!=null && score.length>6){
					if(verbose){System.err.println("Entering condition because score="+Arrays.toString(score));}
					int[] oldArray=score.clone();
					assert(score.length==8);
					int extraPadLeft=score[6];
					int extraPadRight=score[7];
					
					if(r.gaps==null){
						assert(maxLoc-minLoc+1<=msa.maxColumns);
						int newlen=(maxLoc-minLoc+1+extraPadLeft+extraPadRight);
						if(newlen>=msa.maxColumns){
							while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
							while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
							while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
						}
					}else{
						//TODO: In this case the alignment will probably be wrong.
						int greflen=Tools.max(r.bases.length, GapTools.calcGrefLen(minLoc, maxLoc, r.gaps));
						int newlen=(greflen+1+extraPadLeft+extraPadRight);
						if(newlen>=msa.maxColumns){
							while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
							while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
							while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
						}else{
							int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
							extraPadLeft=Tools.max(x, extraPadLeft);
							extraPadRight=Tools.max(x, extraPadRight);
						}
					}
					
					minLoc=Tools.max(0, minLoc-extraPadLeft);
					maxLoc=Tools.min(chacs.maxIndex, maxLoc+extraPadRight);
					if(verbose){System.err.println("Set extraPadLeft="+extraPadLeft+", extraPadRight="+extraPadRight);}
					if(verbose){System.err.println("Set minLoc="+minLoc+", maxLoc="+maxLoc);}
					
					max=msa.fillLimited(r.bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), r.gaps);
					score=(max==null ? null : msa.score(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null));
					
					if(score==null || score[0]<oldArray[0]){
						if(!Shared.anomaly){System.err.println("Read "+r.numericID+": Padded match string alignment result was inferior.  Triple-aligning. :(");}
						
						if(r.gaps==null){
							assert(maxLoc-minLoc+1<=msa.maxColumns);
							int newlen=(maxLoc-minLoc+1+extraPadLeft+extraPadRight);
							if(newlen>=msa.maxColumns){
								while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
								while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
								while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
							}else{
								int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
								extraPadLeft=Tools.max(x, extraPadLeft);
								extraPadRight=Tools.max(x, extraPadRight);
							}
						}else{
							//TODO: In this case the alignment will probably be wrong.
							int greflen=Tools.max(r.bases.length, GapTools.calcGrefLen(minLoc, maxLoc, r.gaps));
							int newlen=(greflen+1+extraPadLeft+extraPadRight);
							if(newlen>=msa.maxColumns){
								while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
								while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
								while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
							}else{
								int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
								extraPadLeft=Tools.max(x, extraPadLeft);
								extraPadRight=Tools.max(x, extraPadRight);
							}
						}
						
						minLoc=Tools.max(0, minLoc-extraPadLeft);
						maxLoc=Tools.min(chacs.maxIndex, maxLoc+extraPadRight);
						max=msa.fillLimited(r.bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), r.gaps);
						score=(max==null ? null : msa.score(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null));
					}
				}
				
				
				if(verbose){System.err.println(Arrays.toString(max));}
				
				if(max!=null){
					r.match=msa.traceback(r.bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], r.gaps!=null);
					r.start=score[1];
					r.stop=score[2];
					r.mapScore=score[0];
					if(verbose){System.err.println("Aligned:\n"+new String(r.bases)+"\n"+chacs.getString(r.start, r.stop)+"\n"+new String(r.match));}
				}else{
					assert(r.match[0]!='X') : r.toText(false);
					assert(r.match[r.match.length-1]!='X') : r.toText(false);
					r.stop=r.start+r.bases.length-1;
					r.mapScore=scoreNoIndel;
				}
			}

			if(r.colorspace()){
				Tools.reverseInPlace(r.bases);
			}else{
				AminoAcid.reverseComplementBasesInPlace(r.bases);
			}
		}
		if(verbose){System.err.println("Final: "+r.start+", "+r.stop+", "+Gene.strandCodes[r.strand()]);}
		
		if(recur && r.stop<chacs.maxIndex && r.start>0 && (r.match[0]=='X' || r.match[0]=='I' || 
				r.match[r.match.length-1]=='Y' || r.match[r.match.length-1]=='X' || r.match[r.match.length-1]=='I')){
			int xy=0;
			for(int i=0; i<r.match.length; i++){
				byte b=r.match[i];
				if(b=='X' || b=='Y' || b=='I'){xy++;}
			}
//			System.err.println("xy = "+xy);
			
			r.gaps=GapTools.fixGaps(r.start, r.stop, r.gaps, Shared.MINGAP);
			
			int p_temp=Tools.min(10+padding+2*xy, (msa.maxColumns-r.bases.length)/2-20);
			
			realign_new(r, msa, p_temp, false, minValidScore, forbidIndels);
		}
	}
	
	/** For some reason realign was making the match string backwards... */
	public static void realign_new(final SiteScore ss, final byte[] bases, final MSA msa, int padding, final boolean recur, int minValidScore, boolean forbidIndels, final long id){
		boolean colorspace=false;
		assert(colorspace==msa.colorspace);
//		final byte[] bases=ss.plus() ? basesP : basesM;
		
		if(verbose){System.err.println("Padding = "+padding+"; msa.maxColumns = "+msa.maxColumns+"; maplen = "+(ss.stop-ss.start+1)+"; gaps = "+Arrays.toString(ss.gaps));}
		
		assert(padding>=0) : padding+", id="+id+", "+ss;
		padding=Tools.min(padding, (msa.maxColumns-bases.length)/2-20);
		if(verbose){System.err.println("Padding = "+padding);}
		assert(padding>=0) : padding+", id="+id+", "+ss;
		padding=Tools.max(padding, 0);
		if(verbose){System.err.println("Padding = "+padding);}
		assert(padding>=0) : padding+", id="+id+", "+ss;

		
		final ChromosomeArray chacs=Data.getChromosome(ss.chrom);
		if(verbose){
			System.err.println("Realigning.");
			System.err.println("Original: "+ss.start+", "+ss.stop+", "+Gene.strandCodes[ss.strand()]);
			if(verbose){System.err.println("Estimated greflen: "+GapTools.calcGrefLen(ss.start, ss.stop, ss.gaps));}
		}
		
		{
			int expectedLen=GapTools.calcGrefLen(ss.start, ss.stop, ss.gaps);
			if(expectedLen>msa.maxColumns-20){
				//TODO: Alternately, I could kill the site.
				ss.stop=ss.start+Tools.min(bases.length+40, msa.maxColumns-20);
				if(ss.gaps!=null){ss.gaps=GapTools.fixGaps(ss.start, ss.stop, ss.gaps, Shared.MINGAP);}
			}
		}
		
		if(ss.start<0){ss.start=0;} //Prevents assertion errors.  This change should be reset by the realignment so it shouldn't mattess.
		if(ss.stop>chacs.maxIndex){ss.stop=chacs.maxIndex;} //Also to prevent a potential assertion error in unpadded references
		assert(0<=ss.start) : "\nchr"+ss.chrom+": ss.start="+ss.start+", ss.stop="+ss.stop+", padding="+padding+
			", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+ss.toText();
		assert(chacs.maxIndex>=ss.stop) : "\nchr"+ss.chrom+": ss.start="+ss.start+", ss.stop="+ss.stop+", padding="+padding+
			", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+ss.toText();

		{
			assert(ss.stop>=ss.start); //Otherwise this is pointless...
			int a=bases.length;
			int b=ss.stop-ss.start+1;
			if(b<a){
				int c=Tools.min(bases.length, a-b+10)/2;
				padding=Tools.max(padding, c+1);
//				if(verbose){System.err.println("Padding = "+padding);}
//				assert(padding>=0) : padding;
			}
		}
		
		if(verbose){System.err.println("Padding = "+padding);}
		
		padding=Tools.min(padding, (msa.maxColumns-Tools.max(bases.length, GapTools.calcGrefLen(ss.start, ss.stop, ss.gaps)))/2-1);
		if(forbidIndels){padding=0;}
		if(verbose){System.err.println("Padding = "+padding);}
		assert(padding>=0) : id+", "+padding;
		
		
		final int maxQ=msa.maxQuality(bases.length);
		final int maxI=msa.maxImperfectScore(bases.length);

		if(ss.strand()==Gene.PLUS){
			assert(maxQ>maxI);

			byte[][] matchR=new byte[1][];
			if(ss.match!=null && ss.match.length==bases.length){
				matchR[0]=ss.match;
			}else{
				//				System.err.println(new String(ss.match));
				matchR[0]=ss.match=new byte[bases.length];
			}
			int scoreNoIndel=msa.scoreNoIndelsAndMakeMatchString(bases, chacs.array, ss.start, matchR);
			ss.match=matchR[0];
			
			assert(0<=ss.start) : "\nchr"+ss.chrom+": ss.start="+ss.start+", ss.stop="+ss.stop+", padding="+padding+
				", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+ss.toText();
			assert(chacs.maxIndex>=ss.stop) : "\nchr"+ss.chrom+": ss.start="+ss.start+", ss.stop="+ss.stop+", padding="+padding+
				", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+ss.toText();
			
			if(verbose){System.err.println("Estimated greflen: "+GapTools.calcGrefLen(ss.start, ss.stop, ss.gaps));}
			
			if(scoreNoIndel>=maxI || forbidIndels){
				if(verbose){System.err.println("Quick match.");}
//				assert(ss.match[0]!='X') : ss.toText();
//				assert(ss.match[ss.match.length-1]!='X') : ss.toText();
				//				assert(ss.stop==ss.start+bases.length-1);
				ss.stop=ss.start+bases.length-1;
				ss.slowScore=scoreNoIndel;
			}else{
				if(verbose){System.err.println("Slow match.");}
				
//				int minLoc=Tools.max(ss.start-padding, chacs.minIndex);
				int minLoc=Tools.max(ss.start-padding, 0); //It's OK to be off the beginning as long as bases prior to the true start are 'N'
				int maxLoc=Tools.min(ss.stop+padding, chacs.maxIndex);
				
				if(verbose){System.err.println("minLoc = "+minLoc+", maxLoc = "+maxLoc);}
				if(verbose){System.err.println("A. Estimated greflen: "+GapTools.calcGrefLen(ss.start, ss.stop, ss.gaps));}
				if(verbose){System.err.println("A. Estimated greflen2: "+GapTools.calcGrefLen(minLoc, maxLoc, ss.gaps));}
				
				//These assertions are not too important... they indicate the read mapped off the end of the chromosome.
				assert(minLoc<=ss.start) : "\nchr"+ss.chrom+": minloc="+minLoc+", maxLoc="+maxLoc+", ss.start="+ss.start+", ss.stop="+ss.stop+", padding="+padding+
					", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+ss.toText();
				assert(maxLoc>=ss.stop) : "\nchr"+ss.chrom+": minloc="+minLoc+", maxLoc="+maxLoc+", ss.start="+ss.start+", ss.stop="+ss.stop+", padding="+padding+
					", chacs.minIndex="+chacs.minIndex+", chacs.maxIndex="+chacs.maxIndex+"\nread:\n"+ss.toText();

				//			System.err.println("Aligning:\n"+new String(bases)+"\n"+chacs.getString(minLoc, maxLoc));
				
				int[] max=null;
				int[] score=null;
				try {
					if(verbose){
						System.err.println("Calling fillLimited(bases, chacs, "+minLoc+", "+maxLoc+", "+
								Tools.max(scoreNoIndel, minValidScore)+", "+(ss.gaps==null ? "null" : Arrays.toString(ss.gaps))+")");
					}
					max=msa.fillLimited(bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), ss.gaps);
					score=(max==null ? null : msa.score(bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], ss.gaps!=null));
					if(verbose){System.err.println("Estimated greflen: "+GapTools.calcGrefLen(ss.start, ss.stop, ss.gaps));}
					
					if(score!=null && score.length>6){
						int[] oldArray=score.clone();
						assert(score.length==8);
						int extraPadLeft=score[6];
						int extraPadRight=score[7];
						
						if(ss.gaps==null){
							assert(maxLoc-minLoc+1<=msa.maxColumns);
							int newlen=(maxLoc-minLoc+1+extraPadLeft+extraPadRight);
							if(newlen>=msa.maxColumns){
								while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
								while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
								while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
							}else{
								int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
								extraPadLeft=Tools.max(x, extraPadLeft);
								extraPadRight=Tools.max(x, extraPadRight);
							}
						}else{
							//TODO: In this case the alignment will probably be wrong.
							int greflen=Tools.max(bases.length, GapTools.calcGrefLen(minLoc, maxLoc, ss.gaps));
							int newlen=(greflen+1+extraPadLeft+extraPadRight);
							if(newlen>=msa.maxColumns){
								while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
								while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
								while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
							}else{
								int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
								extraPadLeft=Tools.max(x, extraPadLeft);
								extraPadRight=Tools.max(x, extraPadRight);
							}
						}
						
						assert(extraPadLeft>=0 && extraPadRight>=0) : extraPadLeft+", "+extraPadRight+"\n"+id+", "+ss;
						minLoc=Tools.max(0, minLoc-extraPadLeft);
						maxLoc=Tools.min(chacs.maxIndex, maxLoc+extraPadRight);

						if(verbose){System.err.println("B. Estimated greflen: "+GapTools.calcGrefLen(ss.start, ss.stop, ss.gaps));}
						if(verbose){System.err.println("B. Estimated greflen2: "+GapTools.calcGrefLen(minLoc, maxLoc, ss.gaps));}
						max=msa.fillLimited(bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), ss.gaps);
						score=(max==null ? null : msa.score(bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], ss.gaps!=null));
						
						if(score==null || score[0]<oldArray[0]){
							if(!Shared.anomaly){System.err.println("Read "+id+": Padded match string alignment result was inferior.  Triple-aligning. :(");}
							
							if(ss.gaps==null){
								assert(maxLoc-minLoc+1<=msa.maxColumns);
								int newlen=(maxLoc-minLoc+1+extraPadLeft+extraPadRight);
								if(newlen>=msa.maxColumns){
									while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
									while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
									while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
								}else{
									int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
									extraPadLeft=Tools.max(x, extraPadLeft);
									extraPadRight=Tools.max(x, extraPadRight);
								}
							}else{
								//TODO: In this case the alignment will probably be wrong.
								int greflen=Tools.max(bases.length, GapTools.calcGrefLen(minLoc, maxLoc, ss.gaps));
								int newlen=(greflen+1+extraPadLeft+extraPadRight);
								if(newlen>=msa.maxColumns){
									while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
									while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
									while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
								}else{
									int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
									extraPadLeft=Tools.max(x, extraPadLeft);
									extraPadRight=Tools.max(x, extraPadRight);
								}
							}
							
							assert(extraPadLeft>=0 && extraPadRight>=0) : extraPadLeft+", "+extraPadRight+"\n"+id+", "+ss;
							minLoc=Tools.max(0, minLoc-extraPadLeft);
							maxLoc=Tools.min(chacs.maxIndex, maxLoc+extraPadRight);
							
							max=msa.fillLimited(bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), ss.gaps);
							score=(max==null ? null : msa.score(bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], ss.gaps!=null));
							
							if(minLoc>0 && maxLoc<chacs.maxIndex && (score==null || score[0]<oldArray[0])){
								if(!Shared.anomaly){System.err.println("Still inferior.");}
								minLoc=Tools.max(ss.start-8, 0); //It's OK to be off the beginning as long as bases prior to the true start are 'N'
								maxLoc=Tools.min(ss.stop+8, chacs.maxIndex);
								max=msa.fillUnlimited(bases, chacs.array, minLoc, maxLoc, ss.gaps);
								score=(max==null ? null : msa.score(bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], ss.gaps!=null));
							}
						}
					}
				} catch (Exception e) {
					System.err.println("Caught exception:\n");
					e.printStackTrace();
					assert(false) : ss.toText();
				}
				
				if(max!=null){
					ss.match=msa.traceback(bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], ss.gaps!=null);
					ss.start=score[1];
					ss.stop=score[2];
					ss.slowScore=score[0];
				}else{
					ss.stop=ss.start+bases.length-1;
					ss.slowScore=scoreNoIndel;
				}
			}
		}else{
			assert(maxQ>maxI);

			byte[][] matchR=new byte[1][];
			if(ss.match!=null && ss.match.length==bases.length){
				matchR[0]=ss.match;
			}else{
				matchR[0]=ss.match=new byte[bases.length];
			}
			
			int scoreNoIndel=msa.scoreNoIndelsAndMakeMatchString(bases, chacs.array, ss.start, matchR);
			ss.match=matchR[0];
			
			if(scoreNoIndel>=maxI || forbidIndels){
				if(verbose){System.err.println("Quick match.");}
				assert(ss.match[0]!='X') : ss.toText();
				assert(ss.match[ss.match.length-1]!='X') : ss.toText();
				ss.stop=ss.start+bases.length-1;
				ss.slowScore=scoreNoIndel;
//				Tools.reverseInPlace(ss.match);
			}else{
				if(verbose){System.err.println("Slow match.");}
				
				int minLoc=Tools.max(ss.start-padding, 0); //It's OK to be off the beginning as long as bases prior to the true start are 'N'
				int maxLoc=Tools.min(ss.stop+padding, chacs.maxIndex);
				if(verbose){System.err.println("Slow match "+minLoc+" ~ "+maxLoc);}

				//These assertions are not too important... they indicate the read mapped off the end of the chromosome.
				assert(minLoc<=ss.start) : "\nchr"+ss.chrom+": "+minLoc+", "+maxLoc+", "+ss.start+", "+ss.stop+
					", "+chacs.minIndex+", "+chacs.maxIndex+"\n"+ss.toText();
				assert(maxLoc>=ss.stop) : "\nchr"+ss.chrom+": "+minLoc+", "+maxLoc+", "+ss.start+", "+ss.stop+
					", "+chacs.minIndex+", "+chacs.maxIndex+"\n"+ss.toText();

				if(verbose){System.err.println("Aligning:\n"+new String(bases)+"\n"+chacs.getString(minLoc, maxLoc));}
				int[] max=msa.fillLimited(bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), ss.gaps);
				if(verbose){System.err.println("Aligned: {rows, maxC, maxS, max} = "+Arrays.toString(max));}
				int[] score=null;
				score=(max==null ? null : msa.score(bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], ss.gaps!=null));
				
				if(score!=null && score.length>6){
					if(verbose){System.err.println("Entering condition because score="+Arrays.toString(score));}
					int[] oldArray=score.clone();
					assert(score.length==8);
					int extraPadLeft=score[6];
					int extraPadRight=score[7];
					
					if(ss.gaps==null){
						assert(maxLoc-minLoc+1<=msa.maxColumns);
						int newlen=(maxLoc-minLoc+1+extraPadLeft+extraPadRight);
						if(newlen>=msa.maxColumns){
							while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
							while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
							while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
						}
					}else{
						//TODO: In this case the alignment will probably be wrong.
						int greflen=Tools.max(bases.length, GapTools.calcGrefLen(minLoc, maxLoc, ss.gaps));
						int newlen=(greflen+1+extraPadLeft+extraPadRight);
						if(newlen>=msa.maxColumns){
							while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
							while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
							while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
						}else{
							int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
							extraPadLeft=Tools.max(x, extraPadLeft);
							extraPadRight=Tools.max(x, extraPadRight);
						}
					}
					
					minLoc=Tools.max(0, minLoc-extraPadLeft);
					maxLoc=Tools.min(chacs.maxIndex, maxLoc+extraPadRight);
					if(verbose){System.err.println("Set extraPadLeft="+extraPadLeft+", extraPadRight="+extraPadRight);}
					if(verbose){System.err.println("Set minLoc="+minLoc+", maxLoc="+maxLoc);}
					
					max=msa.fillLimited(bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), ss.gaps);
					score=(max==null ? null : msa.score(bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], ss.gaps!=null));
					
					if(score==null || score[0]<oldArray[0]){
						if(!Shared.anomaly){System.err.println("Read "+id+": Padded match string alignment result was inferior.  Triple-aligning. :(");}
						
						if(ss.gaps==null){
							assert(maxLoc-minLoc+1<=msa.maxColumns);
							int newlen=(maxLoc-minLoc+1+extraPadLeft+extraPadRight);
							if(newlen>=msa.maxColumns){
								while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
								while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
								while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
							}else{
								int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
								extraPadLeft=Tools.max(x, extraPadLeft);
								extraPadRight=Tools.max(x, extraPadRight);
							}
						}else{
							//TODO: In this case the alignment will probably be wrong.
							int greflen=Tools.max(bases.length, GapTools.calcGrefLen(minLoc, maxLoc, ss.gaps));
							int newlen=(greflen+1+extraPadLeft+extraPadRight);
							if(newlen>=msa.maxColumns){
								while(newlen>=msa.maxColumns && extraPadLeft>extraPadRight){newlen--;extraPadLeft--;}
								while(newlen>=msa.maxColumns && extraPadLeft<extraPadRight){newlen--;extraPadRight--;}
								while(newlen>=msa.maxColumns){newlen-=2;extraPadLeft--;extraPadRight--;}
							}else{
								int x=Tools.min(20, ((msa.maxColumns-newlen)/2)-1);
								extraPadLeft=Tools.max(x, extraPadLeft);
								extraPadRight=Tools.max(x, extraPadRight);
							}
						}
						
						minLoc=Tools.max(0, minLoc-extraPadLeft);
						maxLoc=Tools.min(chacs.maxIndex, maxLoc+extraPadRight);
						max=msa.fillLimited(bases, chacs.array, minLoc, maxLoc, Tools.max(scoreNoIndel, minValidScore), ss.gaps);
						score=(max==null ? null : msa.score(bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], ss.gaps!=null));
					}
				}
				
				
				if(verbose){System.err.println(Arrays.toString(max));}
				
				if(max!=null){
					ss.match=msa.traceback(bases, chacs.array, minLoc, maxLoc, max[0], max[1], max[2], ss.gaps!=null);
					ss.start=score[1];
					ss.stop=score[2];
					ss.slowScore=score[0];
					if(verbose){System.err.println("Aligned:\n"+new String(bases)+"\n"+chacs.getString(ss.start, ss.stop)+"\n"+new String(ss.match));}
				}else{
					assert(ss.match[0]!='X') : ss.toText();
					assert(ss.match[ss.match.length-1]!='X') : ss.toText();
					ss.stop=ss.start+bases.length-1;
					ss.slowScore=scoreNoIndel;
				}
			}
		}
		if(verbose){System.err.println("Final: "+ss.start+", "+ss.stop+", "+Gene.strandCodes[ss.strand()]);}
		
		if(recur && ss.stop<chacs.maxIndex && ss.start>0 && (ss.match[0]=='X' || ss.match[0]=='I' || 
				ss.match[ss.match.length-1]=='Y' || ss.match[ss.match.length-1]=='X' || ss.match[ss.match.length-1]=='I')){
			int xy=0;
			for(int i=0; i<ss.match.length; i++){
				byte b=ss.match[i];
				if(b=='X' || b=='Y' || b=='I'){xy++;}
			}
//			System.err.println("xy = "+xy);
			
			ss.gaps=GapTools.fixGaps(ss.start, ss.stop, ss.gaps, Shared.MINGAP);
			
			int p_temp=Tools.min(10+padding+2*xy, (msa.maxColumns-bases.length)/2-20);
			
			realign_new(ss, bases, msa, p_temp, false, minValidScore, forbidIndels, id);
		}
		ss.setPerfect(bases);
	}
	
	private static final boolean checkArray(byte[] bases){
		for(byte b : bases){
//			assert(b>0) : Arrays.toString(bases);
			if(b<=0){return false;}
		}
		return true;
	}
	

	public Read translateToBasespace(Read r){
//		if(r.strand==Gene.MINUS){
//			r.reverse();
//			r.obj=AminoAcid.reverseComplementBases(r.obj);
//		}
		
		assert(r.colorspace()) : r.toText(false);
		assert(r.start<r.stop) : r.toText(false);
		
		ChromosomeArray chab=Data.getChromosome(r.chrom);
		ChromosomeArray chacs=null;//Data.getChromosome(r.chrom, true);
		
		
		if(verbose){System.err.println("\n"+r.numericID+"\t**************************\n");}

		if(verbose){		
			byte[][] crbmq=new byte[5][];
			crbmq[0]=Arrays.copyOf(r.bases, r.bases.length);
			crbmq[1]=chacs.getBytes(r.start, r.stop);
			crbmq[2]=chab.getBytes(r.start, r.stop+1);
			crbmq[3]=(r.match==null ? null : Arrays.copyOf(r.match, r.match.length));
			crbmq[4]=Arrays.copyOf(r.quality, r.quality.length);
			System.err.println("Original:\n"+toString(crbmq)+"\n");
			assert(checkArray(crbmq[2]));
//			assert(crbmq[3]==null || verifyMatchString(crbmq[0], crbmq[1], crbmq[3], 0));
		}

		final boolean perfect=perfectMatch(r.match);
		
//		assert(r.match!=null);
		
		if(verbose){
			System.err.println(r.toText(false));
		}
		realignByReversingRef(r, msaCS, 4, true); //TODO - can be mostly skipped if perfect==true
		
		
		byte[][] crbmq=new byte[5][];
//		crbmq[0]=Arrays.copyOf(r.bases, r.bases.length);
//		crbmq[1]=chacs.getBytes(r.start, r.stop);
//		crbmq[2]=chab.getBytes(r.start, r.stop+1);
//		crbmq[3]=Arrays.copyOf(r.match, r.match.length);
//		crbmq[4]=Arrays.copyOf(r.quality, r.quality.length);
		
		crbmq[0]=r.bases;
		crbmq[1]=chacs.getBytes(r.start, r.stop);
		crbmq[2]=chab.getBytes(r.start, r.stop+1);
		crbmq[3]=r.match;
		crbmq[4]=r.quality;
		
//		assert(verifyMatchString(crbmq[0], crbmq[1], crbmq[3], 0));
		
		//assert(checkArray(crbmq[2]));
		
//		System.err.println("~\n"+toString(crbmq));
		
		if(verbose){
			System.err.println("Realigned:\n"+toString(crbmq)+"\n");
		}
		
		if(r.strand()==Gene.MINUS){
			Tools.reverseInPlace(crbmq[1]);
			AminoAcid.reverseComplementBasesInPlace(crbmq[2]);
			if(verbose){
				System.err.println("Reversed:\n"+toString(crbmq)+"\n");
			}
			
		}
		//assert(checkArray(crbmq[2]));
//		
		assert(verifyMatchString(crbmq[0], crbmq[1], crbmq[3], 0, true, true)) : "\n"+r.toText(false);
		

		if(verbose){System.err.println("r.start="+r.start+", r.stop="+r.stop+", len="+(r.stop-r.start+1));}
		
//		System.err.println("~\n"+toString(crbmq)+"\n");
		
		int x;
		
//		if(crbmq[3].length<2){return null;}
//		x=trimStart(crbmq, r);
//		assert(checkArray(crbmq[2]));
		
		if(crbmq[3].length<2){return null;}
		x=trimEnd(crbmq, 2, r);
		assert(x>=0);
		//assert(checkArray(crbmq[2]));
		if(verbose && x!=0){
			System.err.println("Trimmed end:");
			System.err.println(toString(crbmq));
		}
		if(crbmq[3].length<2){return null;}
		
//		byte[] original=Arrays.copyOf(crbmq[2], crbmq[2].length);
		
		if(crbmq[3][crbmq[3].length-1]!='m'){
			System.err.println("Failed to trim read "+r.numericID+", x="+x);
			System.err.println(r.toText(false));
			System.err.println(toString(crbmq));
			return null;
		}

		final boolean containsIndels=containsIndels(crbmq[3]);
		final boolean containsNocalls=containsNocalls(crbmq[3]) || 
			containsNocalls(crbmq[0]) || containsNocalls(crbmq[1]) || containsNocalls(crbmq[2]);
		
		if(containsNocalls){
			
//			if(containsIndels){
//				System.err.println("*** Before fixing nocalls:");
//				System.err.println(toString(crbmq));
//			}
			
			x=fixNocallsInline(crbmq, r);
			//assert(checkArray(crbmq[2]));
			if(x<0){
				if(verbose){System.err.println("------------------------ broke decoder ------------------------");}
				return null;
			}
			if(verbose){
				if(x>0){
					System.err.println("Fixed Nocalls:");
					System.err.println(toString(crbmq));
				}
			}

//			if(containsIndels){
//				System.err.println("*** After fixing nocalls:");
//				System.err.println(toString(crbmq)+"\n");
//			}
//			
//			x=fixNocalls(crbmq);
//			assert(checkArray(crbmq[2]));
//			if(x<0){
//				if(verbose){System.err.println("------------------------ broke decoder ------------------------");}
//				return null;
//			}
//			if(verbose){
//				if(x>0){
//					System.err.println("Fixed Nocalls:");
//					System.err.println(toString(crbmq));
//				}
//			}
//			
//
//			if(containsIndels){
//				x=fixNocallsBackward(crbmq);
//				assert(checkArray(crbmq[2]));
//				if(x<0){
//					if(verbose){System.err.println("------------------------ broke decoder ------------------------");}
//					return null;
//				}
//				if(verbose){
//					if(x>0){
//						System.err.println("Fixed Nocalls backwards:");
//						System.err.println(toString(crbmq));
//					}
//				}
//			}
		}
		
		
		
		if(verbose){System.err.println(toString(crbmq));}
		x=fixIndels(crbmq, r);
		if(crbmq[3].length<2){return null;}
		//assert(checkArray(crbmq[2])) : new String(crbmq[2])+"\n"+new String(original)+"\n"+new String(r.match);
		final int indels=x;
		if(x<0){
			if(verbose){System.err.println("------------------------ broke decoder ------------------------");}
			return null;
		}
		if(verbose){
			if(x>0){
				System.err.println("Fixed indels:");
				System.err.println(toString(crbmq));
			}
		}
		assert(crbmq[0].length==crbmq[1].length) : "\n"+toString(crbmq)+"\n"+r.toText(false)+"\n"+chab.maxIndex+"\n";
		assert(crbmq[0].length==crbmq[3].length) : "\n"+toString(crbmq)+"\n"+r.toText(false)+"\n"+chab.maxIndex+"\n";
		assert(crbmq[0].length==(crbmq[2].length-1)) : "\n"+toString(crbmq)+"\n"+r.toText(false)+"\n"+chab.maxIndex+"\n";
		
//		if(containsNocalls && containsIndels){//Indels are gone now, but some nocalls may remain
//			assert(!containsIndels(crbmq[3]));
//			x=fixNocallsInline(crbmq, r);
//			assert(checkArray(crbmq[2]));
//			if(x<0){
//				if(verbose){System.err.println("------------------------ broke decoder ------------------------");}
//				return null;
//			}
//			if(verbose){
//				if(x>0){
//					System.err.println("Fixed Nocalls 2:");
//					System.err.println(toString(crbmq));
//				}
//			}
////			x=fixNocalls(crbmq);
////			assert(checkArray(crbmq[2]));
////			if(x<0){
////				if(verbose){System.err.println("------------------------ broke decoder ------------------------");}
////				return null;
////			}
////			if(verbose){
////				if(x>0){
////					System.err.println("Fixed Nocalls 2:");
////					System.err.println(toString(crbmq));
////				}
////			}
//		}else{assert(!containsNocalls(crbmq[3])) : "\n"+toString(crbmq)+"\n";}
		
		x=fixColorspaceErrors(crbmq, 3, 10);
		r.errors=(x>=0 ? x : 0);
		//assert(checkArray(crbmq[2]));
		if(x<0){
			if(verbose){System.err.println("------------------------ broke decoder ------------------------");}
			return null;
		}
		if(verbose){
			if(x>0){
				System.err.println("Fixed Errors:");
				System.err.println(toString(crbmq));
			}
		}
		x=fixSubs(crbmq);
		//assert(checkArray(crbmq[2]));
		if(x<0){
			if(verbose){System.err.println("------------------------ broke decoder ------------------------");}
			return null;
		}

		if(verbose){
			if(x>0){
				System.err.println("Fixed Substitutions:");
				System.err.println(toString(crbmq));
			}else{
				System.err.println("Final:");
				System.err.println(toString(crbmq));
			}
		}
		
		byte[] qualityBS=translateQuality(crbmq[4]);
		
		int flags=(r.flags&(~Read.COLORMASK));
		
		assert(false) : "TODO: Make sure the next line (read instantiation) is correct.";
		Read r2=new Read(crbmq[2], r.chrom, r.start, r.stop+1, r.id, qualityBS, r.numericID, flags);
		
		r2.originalSite=r.originalSite;
		r2.errors=r.errors;
		r2.mapScore=r.mapScore;
		r2.obj=r.obj;
		r2.sites=r.sites;
		r2.copies=r.copies;
		r2.mapLength=r.mapLength;
		
		if(verbose){
			System.err.println("r2:");
			System.err.println(chab.getString(r2.start, r2.stop)+" (ref) ");
			System.err.println(new String(r2.bases)+" (call) ");
		}
		
		if(indels>0 || (r2.stop-r2.start+1!=r2.bases.length)){
			if(verbose){
				System.err.println("Making slow BS match:");
			}
			realign_new(r2, 4, true, 0, false);
			
//			int padding=4;
//			{
//				int a=r2.bases.length;
//				int b=r2.stop-r2.start+1;
//				if(b<a){
//					int c=Tools.min(r.bases.length, a-b+10)/2;
//					padding=Tools.max(padding, c+1);
//				}
//			}
//			
//			assert(padding>=0);
//			final int minLoc=Tools.max(0, r2.start-padding);
//			final int maxLoc=Tools.min(chab.maxIndex, r2.stop+padding);
//			if(r2.strand()==Gene.PLUS){
//				int[] max=msaBS.fill(r2.bases, chab.array, minLoc, maxLoc);
//				r2.match=msaBS.traceback(r2.bases, chab.array, minLoc, maxLoc, max[0], max[1], max[2]);
//				int[] score=msaBS.score(r2.bases, chab.array, minLoc, maxLoc, max[0], max[1], max[2]);
//				r2.start=score[1];
//				r2.stop=score[2];
//			}else{
//				byte[] bases=AminoAcid.reverseComplementBases(r2.bases);
//				if(verbose){
//					System.err.println("reversed:");
//					System.err.println(chab.getString(minLoc, maxLoc)+" (ref extended) ");
//					System.err.println(chab.getString(r2.start, r2.stop)+" (ref) ");
//					System.err.println(new String(bases)+" (call) ");
//				}
//				int[] max=msaBS.fill(bases, chab.array, minLoc, maxLoc);
//				r2.match=msaBS.traceback(bases, chab.array, minLoc, maxLoc, max[0], max[1], max[2]);
//				int[] score=msaBS.score(bases, chab.array, minLoc, maxLoc, max[0], max[1], max[2]);
//				r2.start=score[1];
//				r2.stop=score[2];
////				Tools.reverseInPlace(r2.match);
//			}
			
		}else{
			if(verbose){
				System.err.println("Making quick BS match:");
			}
			final int maxI=msaBS.maxImperfectScore(r2.bases.length);
			if(r2.strand()==Gene.PLUS){
				byte[][] matchR=new byte[1][];
				if(r2.match!=null && r2.match.length==r2.bases.length){
					matchR[0]=r2.match;
				}else{
//					System.err.println(new String(r2.match));
					matchR[0]=r2.match=new byte[r2.bases.length];
				}
				int scoreNoIndel=msaBS.scoreNoIndelsAndMakeMatchString(r2.bases, chab.array, r2.start, matchR);
				
				//TODO: If scoreNoIndel<maxI, the decode process may have failed.  Should realign slowly.
				
//				assert(r2.stop==r2.start+r2.bases.length-1);
				r2.stop=r2.start+r2.bases.length-1;
			}else{
				byte[] bases=AminoAcid.reverseComplementBases(r2.bases);
				
				byte[][] matchR=new byte[1][];
				if(r2.match!=null && r2.match.length==r2.bases.length){
					matchR[0]=r2.match;
				}else{
//					System.err.println(new String(r2.match));
					matchR[0]=r2.match=new byte[r2.bases.length];
				}
				int scoreNoIndel=msaBS.scoreNoIndelsAndMakeMatchString(bases, chab.array, r2.start, matchR);
				
				//TODO: If scoreNoIndel<maxI, the decode process may have failed.  Should realign slowly.
				
//				assert(r2.stop==r2.start+r2.bases.length-1);
				r2.stop=r2.start+r2.bases.length-1;
				
//				Tools.reverseInPlace(r2.match);
			}
		}
		
		return r2;
	}
	
	
	public static byte[] translateQuality(byte[] qcs){
		byte[] qbs=new byte[qcs.length+1];
		qbs[0]=qcs[0];
		qbs[qbs.length-1]=qcs[qcs.length-1];
		for(int i=1; i<qcs.length; i++){
			int x=Tools.min(qcs[i-1], qcs[i]);
			int y=Tools.max(qcs[i-1], qcs[i]);
			qbs[i]=(byte) ((3*x+y)/4);
		}
		return qbs;
	}
	
	private static int fixIndels(byte[][] crbmq, Read r){
		
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		byte[] quality=crbmq[4];
		
		for(int i=0; i<match.length; i++){
			if(match[i]=='X' || match[i]=='Y'){
//				assert(false) : "\n"+new String(colors)+"\n"+new String(colorRef)+"\n"+new String(baseRef)+"\n"+new String(match)+"\n";
				
				assert(false) : "TODO: Truncate ends.\n"+toString(crbmq)+"\n";
				
				match[i]='I';
//				match[i]='S';
			}
		}
		
		int fixed=0;
		
		for(int loc=0, refloc=0, mloc=0; mloc<match.length; mloc++){
			byte b=match[mloc];
			boolean ok=true;
			if(b=='S' || b== 'm' || b=='N'){
				loc++;
				refloc++;
			}else if(b=='D'){
				fixed++;
				ok=fixDeletion(crbmq, mloc, r);
				match=crbmq[3];
				assert(ok);
			}
//			else if(b=='I' || b=='X' || b=='Y'){
//				ok=fixInsertion(crbmq, mloc);
//				assert(ok);
//			}
			else if(b=='I'){
				fixed++;
				ok=fixInsertion(crbmq, mloc);
				match=crbmq[3];
//				assert(ok);
			}else{
				assert(false) : ""+((char)b);
			}
			if(!ok){
				return -1;
			}
		}
		
		colors=crbmq[0];
		colorRef=crbmq[1];
		baseRef=crbmq[2];
		match=crbmq[3];
		
		assert(baseRef.length==colorRef.length+1);
		if(colorRef.length>colors.length){
			colorRef=Arrays.copyOf(colorRef, colors.length);
			baseRef=Arrays.copyOf(baseRef, colorRef.length+1);
			crbmq[1]=colorRef;
			crbmq[2]=baseRef;
		}
		
		return fixed;
	}
	
	private static boolean fixDeletion(final byte[][] crbmq, int loc, Read r){
		
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		byte[] quality=crbmq[4];
		
		assert(match[loc]=='D') : loc;
		
		int len=1;
		for(int i=loc+1; i<match.length; i++){
			byte b=match[i];
			if(b=='D'){
				len++;
			}else{
				break;
			}
		}
		int b=loc+len-1;
		
		//TODO
		if(loc<=1 || b>match.length-2){return false;} //Indels on very ends need to be processed differently
		
		//Deletion is from a to b, inclusive.  Note that basespace coords are +1 from colorspace coords.

		byte[] colorRef2=new byte[colorRef.length-len];
		byte[] baseRef2=new byte[baseRef.length-len];
		byte[] match2=new byte[match.length-len];
		
		assert(loc<colorRef2.length) : "TODO: Seems odd... "+loc+", "+colorRef2.length+", "+match;
		assert(baseRef2.length==colorRef2.length+1);
		
		for(int i=0; i<=loc; i++){
			colorRef2[i]=colorRef[i];
			baseRef2[i]=baseRef[i];
			match2[i]=match[i];
		}
		
		for(int i=loc+1; i<baseRef2.length; i++){
			baseRef2[i]=baseRef[i+len];
		}
		for(int i=loc+1; i<colorRef2.length; i++){
			colorRef2[i]=colorRef[i+len];
		}
		for(int i=loc+1; i<match2.length; i++){
			match2[i]=match[i+len];
		}
		
		
		colorRef2[loc]=AminoAcid.baseToColor(baseRef2[loc], baseRef2[loc+1]);
		if(colorRef2[loc]==colors[loc]){
			match2[loc]='m';
		}else{
			assert(colorRef2[loc]!='N' && colors[loc]!='N') : "TODO\n"+r.toText(false)+"\n"+toString(crbmq)+"\n";
			match2[loc]='S';
		}
		
		crbmq[1]=colorRef2;
		crbmq[2]=baseRef2;
		crbmq[3]=match2;
		crbmq[4]=quality;
		
		return true;
	}
	
	private static boolean fixInsertion(final byte[][] crbmq, int loc){
		
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		byte[] quality=crbmq[4];
		
		assert(match[loc]=='I');
		
		int len=1;
		for(int i=loc+1; i<match.length; i++){
			byte b=match[i];
			if(b=='I'){
				len++;
			}else{
				break;
			}
		}
		int b=loc+len-1;
		
		byte[] colorRef2=new byte[colorRef.length+len];
		byte[] baseRef2=new byte[baseRef.length+len];
		byte[] match2=new byte[match.length]; //TODO:  Unnecessary duplication\
		
		//TODO
//		if(b>match.length-2){return false;} //Indels on very ends need to be processed differently
		
		
		//Deletion is from a to b, inclusive.  Note that basespace coords are +1 from colorspace coords
		
		assert(loc<colorRef2.length) : "TODO: Seems odd... "+loc+", "+colorRef2.length+", "+match;
		assert(baseRef2.length==colorRef2.length+1);
		
		//Fill first half
		for(int i=0; i<loc; i++){
			colorRef2[i]=colorRef[i];
			baseRef2[i]=baseRef[i];
			match2[i]=match[i];
		}
		baseRef2[loc]=baseRef[loc];
		
		//Fill last half
		for(int i=loc+1; i<colorRef.length; i++){
			colorRef2[i+len]=colorRef[i];
		}
		for(int i=loc; i<baseRef.length; i++){
			baseRef2[i+len]=baseRef[i];
		}
		for(int i=loc+1; i<match.length; i++){
			match2[i]=match[i];
		}
		
		//Now, just fill in the inserted portion
		if(verbose){
			System.err.println("loc="+loc+", colorRef2="+colorRef2.length+", colors="+colors.length+", match2="+match2.length);
			System.err.println("max="+Tools.min(loc+len, Tools.min(colorRef2.length, colors.length)-1));
		}
		for(int i=loc, max=Tools.min(loc+len, Tools.min(colorRef2.length, colors.length)-1); i<=max; i++){
			colorRef2[i]=colors[i];
		}
		for(int i=loc, max=Tools.min(loc+len, match2.length-1); i<=max; i++){
			match2[i]='m';
		}
		


		if(loc==0){
			for(int i=(Tools.min(loc+len, colorRef.length-1)); i>=0; i--){
				if(DISCARD_NOCALLED_INSERTIONS && colorRef2[i]=='N'){return false;} //Fail.
				
//				if(colorRef2[i]=='N'){System.err.println("Keeping no-called insertion:\n"+toString(crbmq)+"\n");}
				
//				assert(colorRef2[i]!='N') : "TODO\n"+toString(crbmq)+"\n";
				
				//			System.err.println(""+(char)AminoAcid.colorToBase(baseRef2[i-1], colorRef2[i-1]));
				//			System.err.println(""+(char)baseRef2[i-1]);
				//			System.err.println(""+(char)colorRef2[i-1]);
				//			System.err.println(new String(baseRef)+"\t"+new String(colorRef)+"\t"+new String(baseRef2)+"\t"+new String(colorRef2));
				//			System.err.println("loc="+loc+", i="+i);
				baseRef2[i]=AminoAcid.colorToBase(baseRef2[i+1], colorRef2[i]);
			}
		}else{

			for(int i=loc+1, max=loc+len; i<=max; i++){
				if(DISCARD_NOCALLED_INSERTIONS && colorRef2[i-1]=='N'){return false;} //Fail.
				
//				if(colorRef2[i-1]=='N'){System.err.println("Keeping no-called insertion:\n"+toString(crbmq)+"\n");}
				
//				assert(colorRef2[i-1]!='N') : "TODO\n"+toString(crbmq)+"\n";
				
				//			System.err.println(""+(char)AminoAcid.colorToBase(baseRef2[i-1], colorRef2[i-1]));
				//			System.err.println(""+(char)baseRef2[i-1]);
				//			System.err.println(""+(char)colorRef2[i-1]);
				//			System.err.println(new String(baseRef)+"\t"+new String(colorRef)+"\t"+new String(baseRef2)+"\t"+new String(colorRef2));
				//			System.err.println("loc="+loc+", i="+i);
				baseRef2[i]=AminoAcid.colorToBase(baseRef2[i-1], colorRef2[i-1]);
			}
		}
		
		
		crbmq[1]=colorRef2;
		crbmq[2]=baseRef2;
		crbmq[3]=match2;
		crbmq[4]=quality;
		
		return true;
	}
	
	
//	private static int fixNocalls(final byte[][] crbmq){
//		byte[] colors=crbmq[0];
//		byte[] colorRef=crbmq[1];
//		byte[] baseRef=crbmq[2];
//		
//		int fixedRef=0;
//		int fixedCall=0;
//		assert(colors.length==colorRef.length) : "\n"+Arrays.toString(colors)+"\n"+Arrays.toString(colorRef)+
//			"\n"+new String(baseRef)+"\n"+new String(crbmq[3])+"\n";
//		for(int i=0; i<colors.length; i++){
//			if(colors[i]=='N' || colors[i]=='.'){
//				colors[i]=colorRef[i];
//				fixedCall++;
//			}
//			if(colorRef[i]=='N' || colorRef[i]=='.'){
//				colorRef[i]=colors[i];
//				fixedRef++;
//			}
//		}
//		
//		if(fixedRef>0){
//			for(int i=1; i<colorRef.length; i++){
//				if(baseRef[i]=='N'){
//					baseRef[i]=AminoAcid.colorToBase(baseRef[i-1], colorRef[i-1]);
//				}
//			}
//			for(int i=colorRef.length-2; i>=0; i--){
//				if(baseRef[i]=='N'){
//					baseRef[i]=AminoAcid.colorToBase(baseRef[i+1], colorRef[i+1]);
//				}
//			}
//		}
//		return fixedRef+fixedCall;
//	}
	
	
	private static int fixNocallsInline(final byte[][] crbmq, Read read){
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		
		int fixedRef=0;
		int fixedCall=0;
		
//		boolean indels=false;
//		
//		int indexOfIndel=colors.length;
//		for(int i=0; i<match.length; i++){
//			if(match[i]=='I' || match[i]=='X' || match[i]=='Y' || match[i]=='D'){
//				indels=true;
//				indexOfIndel=i;
//				break;
//			}
//		}
		
		
		for(int mi=0, ci=0, ri=0; mi<match.length; mi++){

			assert(ci<colors.length) : "\n"+read.toText(false)+"\n"+toString(crbmq);
			
			if(ri>=colorRef.length){
				System.err.println("Failed fixNocallsInline for read "+read.numericID);
				System.err.println(read.toText(false));
				System.err.println(toString(crbmq));
				return -1;
			}
			
			assert(ri<colorRef.length) : "\n"+read.toText(false)+"\n"+toString(crbmq);
			
			final byte m=match[mi];
			final byte c=colors[ci];
			final byte r=colorRef[ri];
			
			
			if(m=='m' || m=='S' || m=='N' || m=='X'){
				
				if(c=='N' || c=='.'){
					if(r!='N' && r!='.'){
						colors[ci]=r;
						fixedCall++;
//						match[mi]='m';
					}
				}
				if(r=='N' || r=='.'){
					if(c!='N' && c!='.'){
						colorRef[ri]=c;
						fixedRef++;
						match[mi]='m';
					}
				}
				if(m=='X'){//Not sure about this
					if(c=='N' || c=='.'){
						match[mi]='N';
					}
					match[mi]='m';
				}
				
				ci++;
				ri++;
			}else if(m=='D'){
				ri++;
			}else if(m=='I'){
				ci++;
			}else{
				assert(false) : "m="+(char)m+"\n"+read.toText(false)+"\n"+toString(crbmq);
			}
				
//			assert(m!='Y') : "m="+(char)m+"\n"+read.toText(false)+"\n"+toString(crbmq);
		}
		
		if(fixedRef>0){
			{//forward

				for(int mi=0, ri=0; mi<match.length; mi++){
					
					assert(ri<colorRef.length) : "\n"+read.toText(false)+"\n"+toString(crbmq);

					byte m=match[mi];
					byte r=colorRef[ri];

					if(m=='m' || m=='S' || m=='N'){
						
						if(baseRef[ri]=='N'){
							baseRef[ri]=AminoAcid.colorToBase(baseRef[ri+1], r);
						}
						if(baseRef[ri+1]=='N'){
							baseRef[ri+1]=AminoAcid.colorToBase(baseRef[ri], r);
						}
						ri++;
					}else if(m=='D'){
						ri++;
					}else if(m=='I'){
					}else{
						assert(false) : "m="+(char)m+"\n"+read.toText(false)+"\n"+toString(crbmq);
					}
				}
			}
			

			{//reverse

				for(int mi=match.length-1, ri=colorRef.length-1; mi>=0; mi--){
					
					assert(ri>=0) : "\n"+read.toText(false)+"\n"+toString(crbmq);

					byte m=match[mi];
					byte r=colorRef[ri];

					if(m=='m' || m=='S' || m=='N'){
						
						if(baseRef[ri]=='N'){
							baseRef[ri]=AminoAcid.colorToBase(baseRef[ri+1], r);
						}
						if(baseRef[ri+1]=='N'){
							baseRef[ri+1]=AminoAcid.colorToBase(baseRef[ri], r);
						}
						ri--;
					}else if(m=='D'){
						ri--;
					}else if(m=='I'){
					}else{
						assert(false) : "m="+(char)m+"\n"+read.toText(false)+"\n"+toString(crbmq);
					}
				}
			}
		}
		return fixedRef+fixedCall;
	}
	
	
	private static int fixNocalls(final byte[][] crbmq){
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		
		int fixedRef=0;
		int fixedCall=0;
		
		boolean indels=false;
		
		int indexOfIndel=colors.length;
		for(int i=0; i<match.length; i++){
			if(match[i]=='I' || match[i]=='X' || match[i]=='Y' || match[i]=='D'){
				indels=true;
				indexOfIndel=i;
				break;
			}
		}
		
//		assert(colors.length==colorRef.length) : "\n"+Arrays.toString(colors)+"\n"+Arrays.toString(colorRef)+
//			"\n"+new String(baseRef)+"\n"+new String(crbmq[3])+"\n";
		for(int i=0; i<indexOfIndel; i++){
//			if(match[i]=='I' || match[i]=='X' || match[i]=='Y' || match[i]=='D'){
//				indels=true;
//				break;
//			}
			if(colors[i]=='N' || colors[i]=='.'){
				if(colorRef[i]!='N' && colorRef[i]!='.'){
					colors[i]=colorRef[i];
					fixedCall++;
					assert(match[i]!='I' && match[i]!='D') : toString(crbmq);
					match[i]='m';
				}
			}
			if(colorRef[i]=='N' || colorRef[i]=='.'){
				if(colors[i]!='N' && colors[i]!='.'){
					colorRef[i]=colors[i];
					fixedRef++;
					assert(match[i]!='I' && match[i]!='D') : toString(crbmq);
					match[i]='m';
				}
			}
		}
		
		assert(indels || colors.length==colorRef.length) : "\n"+toString(crbmq)+"\n";
		
		if(fixedRef>0){

			for(int i=1; i<indexOfIndel; i++){
				if(match[i]=='I' || match[i]=='X' || match[i]=='Y' || match[i]=='D'){
					assert(false);
					break;
				}
				if(baseRef[i]=='N'){
					baseRef[i]=AminoAcid.colorToBase(baseRef[i-1], colorRef[i-1]);
				}
			}
			if(!indels){
				for(int i=colorRef.length-2; i>=0; i--){
					if(baseRef[i]=='N'){
						baseRef[i]=AminoAcid.colorToBase(baseRef[i+1], colorRef[i+1]);
					}
				}
			}
		}
		return fixedRef+fixedCall;
	}
	
	
	private static int fixNocallsBackward(final byte[][] crbmq){
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		
		int fixedRef=0;
		int fixedCall=0;
		
		boolean indels=false;

		int indexOfIndelCall=0;
		int indexOfIndelRef=0;
		int indexOfIndelMatch=0;
		for(int i=match.length-1; i>=0; i--){
			if(match[i]=='I' || match[i]=='X' || match[i]=='Y' || match[i]=='D'){
				indels=true;
				int safe=(match.length-1)-i;
				indexOfIndelMatch=i;
				indexOfIndelCall=colors.length-safe;
				indexOfIndelRef=colorRef.length-safe;
//				System.err.println("indexOfIndelMatch="+indexOfIndelMatch+
//						"\nindexOfIndelCall="+indexOfIndelCall+
//						"\nindexOfIndelRef="+indexOfIndelRef+
//						"\nsafe="+safe);
				break;
			}
		}
		
//		assert(colors.length==colorRef.length) : "\n"+Arrays.toString(colors)+"\n"+Arrays.toString(colorRef)+
//			"\n"+new String(baseRef)+"\n"+new String(crbmq[3])+"\n";
		for(int i=colors.length-1, j=colorRef.length-1, k=match.length-1; i>=indexOfIndelCall; i--, j--, k--){
//			if(match[i]=='I' || match[i]=='X' || match[i]=='Y' || match[i]=='D'){
//				indels=true;
//				break;
//			}
			if(colors[i]=='N' || colors[i]=='.'){
				if(colorRef[j]!='N' && colorRef[j]!='.'){
					colors[i]=colorRef[j];
					fixedCall++;
					assert(match[k]!='I' && match[k]!='D') : "i="+i+", j="+j+", k="+k+"\n"+toString(crbmq);
					match[k]='m';
				}
			}
			if(colorRef[j]=='N' || colorRef[j]=='.'){
				if(colors[i]!='N' && colors[i]!='.'){
					colorRef[j]=colors[i];
					fixedRef++;
					assert(match[k]!='I' && match[k]!='D') : "i="+i+", j="+j+", k="+k+"\n"+toString(crbmq);
					match[k]='m';
				}
			}
		}
		
		assert(indels || colors.length==colorRef.length) : "\n"+toString(crbmq)+"\n";
		
		if(fixedRef>0){
			if(!indels){
				for(int i=1; i<colorRef.length; i++){
					if(baseRef[i]=='N'){baseRef[i]=AminoAcid.colorToBase(baseRef[i-1], colorRef[i-1]);}
				}
			}
			for(int i=colorRef.length-2; i>=indexOfIndelRef; i--){
				if(match[i]=='I' || match[i]=='X' || match[i]=='Y' || match[i]=='D'){
					assert(false);
					break;
				}
				if(baseRef[i]=='N'){
					baseRef[i]=AminoAcid.colorToBase(baseRef[i+1], colorRef[i+1]);
				}
			}
		}
		return fixedRef+fixedCall;
	}
	
	public static boolean perfectMatch(final byte[] match){
		if(match==null){return false;}
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b!='m'){return false;}
		}
		return true;
	}

	private static boolean containsIndels(final byte[] match){
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='I' || b=='D' || b=='X' || b=='Y'){return true;}
		}
		return false;
	}
	
	private static boolean containsNocalls(final byte[] match){
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='N' || b=='X' || b=='Y'){return true;}
		}
		return false;
	}
	
	private static boolean containsXY(final byte[] match){
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='X' || b=='Y'){return true;}
		}
		return false;
	}
	
	
	private static int fixColorspaceErrors(final byte[][] crbmq, int thresh, int minQuality){
		
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		byte[] quality=crbmq[4];

		assert(colors.length==colorRef.length);
		assert(colors.length==match.length);
		
		int total=0;
		int corrected=1;
		
		//First change ALL very-low-quality colors to ref 
		for(int i=0; i<colors.length; i++){
			if(colors[i]!=colorRef[i] && quality[i]<=4){
				if(verbose){System.err.println("Fixed error of quality "+quality[i]);}
				total++;
				colors[i]=colorRef[i];
				match[i]='m';
				quality[i]=1;
			}
		}
		
		if(verbose && total>0){
			System.err.println("Fixed qual-4 Errors:");
			System.err.println(toString(crbmq));
		}

		while(corrected>0){
			corrected=0;
			for(int i=0; i<colors.length; i++){
				if(colors[i]!=colorRef[i]){
					int dist=distToMismatch(colors, colorRef, i, thresh);

					boolean correct=false;
					if(dist==1){
						if(quality[i]<minQuality-2){correct=true;}
					}else if(dist>thresh || quality[i]<minQuality){
						correct=true;
					}

					if(correct){
						if(verbose){System.err.println("Fixed error of quality "+quality[i]+", dist "+dist);}
						corrected++;
						total++;
						colors[i]=colorRef[i];
						match[i]='m';
						quality[i]=1;
					}
				}
			}
		}
		return total;
	}
	
	//TODO: Add support for deletions
	/** thresh: Must see this many consecutive 'm' to stop. */ 
	private static int trimEnd(final byte[][] crbmq, int thresh, Read r){
		
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		byte[] quality=crbmq[4];

//		if(match[0]=='m' && match[match.length-1]=='m'){return 0;}
		{
//			byte a=match[0], b=match[1], c=match[match.length-1];
//			if(a=='m' && b=='m' && c=='m' || c=='S'){return 0;}
			
//			System.err.println(new String(match));
			byte a=match[match.length-1], b=match[match.length-2];
			
//			System.err.println("a="+(char)a+", b="+(char)b);
//			System.err.println("X");
			if((a=='m' || a=='D') && (b=='m')){return 0;}
//			System.err.println("Y");
			if((a=='m' || a=='D') && quality[quality.length-1]>=22 && quality[quality.length-2]>=22){return 0;}
//			System.err.println("Z");
		}

		int last=match.length-1;
		int minBadIndex=last;
		int mcount=0;

		int insertions=0;

		while(last>1 && mcount<thresh){
			byte c=match[last];
			if(c=='m'){mcount++;}
			else if(match[last]=='S' || match[last]=='I' || match[last]=='N'
				|| match[last]=='X'|| match[last]=='Y'){
				minBadIndex=last;
				mcount=0;
				if(match[last]=='I' || match[last]=='X'|| match[last]=='Y'){
					insertions++;
				}
			}else{
				break;
			}
			last--;
		}

		final int trim=match.length-minBadIndex;

		int trim2=insertions-trim;

		colors=Arrays.copyOf(colors, colors.length-trim);
		if(trim2!=0){
			colorRef=Arrays.copyOf(colorRef, colorRef.length-trim+insertions);
			baseRef=Arrays.copyOf(baseRef, baseRef.length-trim+insertions);
		}
		match=Arrays.copyOf(match, match.length-trim);
		quality=Arrays.copyOf(quality, quality.length-trim);

		crbmq[0]=colors;
		crbmq[1]=colorRef;
		crbmq[2]=baseRef;
		crbmq[3]=match;
		crbmq[4]=quality;

		if(r.strand()==Gene.PLUS){
			r.stop-=(trim-insertions);
		}else{
			r.start+=(trim-insertions);
		}
		
//		System.err.println(new String(match));
		return trim;
	}
	
	/** thresh: Must see this many consecutive 'm' to stop. */ 
	private static int trimStart(final byte[][] crbmq, int thresh, Read r){
		assert(false) : "TODO";
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		byte[] quality=crbmq[4];

//		if(match[0]=='m' && match[match.length-1]=='m'){return 0;}
		{
			byte a=match[0];
			if(a=='m'){return 0;}
		}
		if(match[0]=='S' || match[0]=='I' || match[0]=='N'
			|| match[0]=='X'|| match[0]=='Y'){
			int last=match.length-1;
			int minBadIndex=last;
			int mcount=0;
			
			int insertions=0;
			
			while(last>1 && mcount<thresh){
				byte c=match[last];
				if(c=='m'){mcount++;}
				else if(match[last]=='S' || match[last]=='I' || match[last]=='N'
					|| match[last]=='X'|| match[last]=='Y'){
					minBadIndex=last;
					mcount=0;
					if(match[last]=='I' || match[last]=='X'|| match[last]=='Y'){
						insertions++;
					}
				}else{
					break;
				}
				last--;
			}
			
			final int trim=match.length-minBadIndex;
			
			int trim2=insertions-trim;
			
			colors=Arrays.copyOf(colors, colors.length-trim);
			if(trim2!=0){
				colorRef=Arrays.copyOf(colorRef, colorRef.length-trim+insertions);
				baseRef=Arrays.copyOf(baseRef, baseRef.length-trim+insertions);
			}
			match=Arrays.copyOf(match, match.length-trim);
			quality=Arrays.copyOf(quality, quality.length-trim);
			
			crbmq[0]=colors;
			crbmq[1]=colorRef;
			crbmq[2]=baseRef;
			crbmq[3]=match;
			crbmq[4]=quality;
			
			if(r.strand()==Gene.PLUS){
				r.stop-=(trim-insertions);
			}else{
				r.start+=(trim-insertions);
			}
			
			return trim;
		}
		return 0;
	}
	
	
//	private static int trimStart(final byte[][] crbmq, Read r){
//		
//		byte[] colors=crbmq[0];
//		byte[] colorRef=crbmq[1];
//		byte[] baseRef=crbmq[2];
//		byte[] match=crbmq[3];
//		byte[] quality=crbmq[4];
//		
//		if(match[0]=='m' || match[0]=='S'){return 0;}
//		
//		int index=0;
//		int insertions=0;
//		while(index<match.length && (match[index]=='I' || match[index]=='X')){
//			if(match[index]=='I'){insertions++;}
//			index++;
//		}
//		if(index==0){return 0;}
//		System.err.println("*** "+r.toText(false));
//		System.err.println(toString(crbmq));
//		
//		int start2=index-insertions;
//		
//		colors=Arrays.copyOfRange(colors, index, colors.length);
//		if(start2!=0){
//			colorRef=Arrays.copyOfRange(colorRef, index-insertions, colorRef.length);
//			baseRef=Arrays.copyOfRange(baseRef, index-insertions, baseRef.length);
//		}
//		match=Arrays.copyOfRange(match, index, match.length);
//		quality=Arrays.copyOfRange(quality, index, quality.length);
//
//		
//		crbmq[0]=colors;
//		crbmq[1]=colorRef;
//		crbmq[2]=baseRef;
//		crbmq[3]=match;
//		crbmq[4]=quality;
//
//		System.err.println("\n"+toString(crbmq));
//		
//		if(r.strand()==Gene.PLUS){
//			r.start+=(index-insertions);
//		}else{
//			r.stop-=(index-insertions);
//		}
//		
//		return index;
//	}
	
	
	private static int fixSubs(final byte[][] crbmq){
		
		byte[] colors=crbmq[0];
		byte[] colorRef=crbmq[1];
		byte[] baseRef=crbmq[2];
		byte[] match=crbmq[3];
		byte[] quality=crbmq[4];

		assert(colors.length==colorRef.length) : "\n"+toString(crbmq);
		assert(colors.length==match.length) : "\n"+toString(crbmq);
		assert(colors.length==baseRef.length-1) : "\n"+toString(crbmq);
		
		int first=match.length-1, last=0;
		
		for(int i=0; i<match.length; i++){
			if(match[i]=='S'){
				first=Tools.min(first, i);
				last=Tools.max(last, i);
			}
		}
		
		if(verbose){System.err.println("First="+first+", last="+last);}
		
		if(first>last){return 0;} //No subs
		
		
		if(last>=colors.length-1 && first==0){
			return -1; //Cannot decode
		}else if(first>0){ //Go right only
			if(verbose){System.err.println("max="+Tools.min(last+1, colors.length));}
			for(int i=first, max=Tools.min(last+1, colors.length); i<max; i++){
				match[i]='m';
				if(colors[i-1]!='N' && baseRef[i-1]!='N'){
					baseRef[i]=AminoAcid.colorToBase(baseRef[i-1], colors[i-1]);
				}//else do nothing
			}
			if(last==match.length-1 && colors[last]!='N' && baseRef[last]!='N'){
				baseRef[last+1]=AminoAcid.colorToBase(baseRef[last], colors[last]);
			}
		}else if(first==0){ //Left only
			for(int i=last, min=Tools.max(0, first); i>=min; i--){
				match[i]='m';
				if(colors[i]!='N' && baseRef[i+1]!='N'){
					baseRef[i]=AminoAcid.colorToBase(baseRef[i+1], colors[i]);
					assert(baseRef[i]!='N') : i+", "+colors[i]+", "+(char)baseRef[i+1];
				}//else do nothing
			}
		}
		
		return 1;
		
	}
	
	
	private static int distToMismatch(byte[] colors, byte[] colorRef, int loc, int limit) {
		int min=limit+1;
		int left=Tools.max(0, loc-limit);
		int right=Tools.min(colors.length, loc+limit+1);
		for(int i=left; i<loc; i++){
			if(colors[i]!=colorRef[i]){
				min=Tools.min(min, loc-i);
			}
		}
		for(int i=loc+1; i<right; i++){
			if(colors[i]!=colorRef[i]){
				min=Tools.min(min, i-loc);
			}
		}
		return min;
	}
	
	
	public static boolean verifyMatchString2(Read r, boolean loud){
		int maxVars=0;
		
		assert(r.mapped());
		assert(r.valid());
		if(r.match==null){return false;}
		if(r.match.length<r.bases.length){return false;}
		
		byte last='m';
		for(int i=0; i<r.match.length; i++){
			byte b=r.match[i];
			if(b=='X' || b=='Y'){
//				assert(false) : read.toText(false);
//				b=r.match[i]='I';
			} //TODO: Should not be needed, if reads are trimmed...
			
			if(b!='m' && b!=last){
				maxVars++;
			}
			last=b;
		}
		
		if(maxVars==0){
			assert(r.match.length==r.bases.length);
			return true;
		}
		
//		byte[] original=Arrays.copyOf(call, call.length);
		if(r.strand()==Gene.MINUS){
			if(r.colorspace()){
				Tools.reverseInPlace(r.bases);
			}else{
				AminoAcid.reverseComplementBasesInPlace(r.bases);
			}
			Tools.reverseInPlace(r.quality);
		}
		
		
		//assert(checkArray(call)) : 
//			"\n"+new String(original)+"\n"+new String(Tools.reverseAndCopy(call))+"\n"+
//			"\n"+Arrays.toString(original)+"\n"+Arrays.toString(Tools.reverseAndCopy(call))+"\n";
		
//		assert(false) : "TODO: ensure read is aligned with forward strand.";
		
		ChromosomeArray cha=Data.getChromosome(r.chrom);

		
		boolean b=true;
		try{
			b=(verifyMatchString(r.bases, cha.array, r.match, r.start, loud, r.colorspace()));
		}catch(Exception e){
			System.err.println(e);
			System.err.println("This read failed verifyMatchString:\n"+r.toText(false)+"\n");
			b=true;//ignores the problem.
		}
		
		if(r.strand()==Gene.MINUS){
			if(r.colorspace()){
				Tools.reverseInPlace(r.bases);
			}else{
				AminoAcid.reverseComplementBasesInPlace(r.bases);
			}
			Tools.reverseInPlace(r.quality);
		}
		return b;
	}
	
	
	public static boolean verifyMatchString(byte[] call, byte[] ref, byte[] match, int rstart, boolean loud, boolean colorspace){
		
		boolean ok=true;
		for(int ci=0, mi=0, ri=rstart; ok && mi<match.length; mi++){
			byte m=match[mi];
			byte c=(m=='D' ? (byte)'?' : call[ci]);
//			byte r=((m=='I' || m=='X' || m=='Y') ? (byte)'?' : ref[ri]);
			byte r=((m=='I' || m=='X' || m=='Y') ? (byte)'?' : ((ri>=0 && ri<ref.length) ? ref[ri] : (byte)'N'));
			
			if(m=='m' || m=='s'){
				ok=c==r;
				ci++;
				ri++;
			}else if(m=='D'){
				ri++;
			}else if(m=='I' || m=='X' || m=='Y'){
				ci++;
			}else if(m=='S'){
				ok=c!=r;
				ci++;
				ri++;
			}else if(m=='N'){
				ok=(c=='N' || r=='N');
				ci++;
				ri++;
			}else{
				assert(false) : (char)m;
			}
			
		}
		
		if(!ok && loud){
			System.err.println("NOT OK!");
			if(call[0]<4){
				if(ref.length>400){
					System.err.println(toStringCS(Arrays.copyOfRange(ref, rstart, rstart+call.length))+" (ref)");
				}else{
					System.err.println(toStringCS(ref)+" (ref)");
				}
				System.err.println(toStringCS(call)+" (call)");
				System.err.println(new String(match));
			}else{
				if(ref.length>400){
					System.err.println(new String(Arrays.copyOfRange(ref, rstart, rstart+call.length))+" (ref)");
				}else{
					System.err.println(new String(ref)+" (ref)");
				}
				System.err.println(new String(call)+" (call)");
				System.err.println(new String(match));
			}
		}
		
		if(!ok){
			
			ok=true;
			
			if(loud){System.err.println("Attempting to fix and skip error.");}
			for(int ci=0, mi=0, ri=rstart; mi<match.length; mi++){
				byte m=match[mi];
				byte c=(m=='D' ? (byte)'?' : call[ci]);
//				byte r=((m=='I' || m=='X' || m=='Y') ? (byte)'?' : ref[ri]);
				byte r=((m=='I' || m=='X' || m=='Y') ? (byte)'?' : ((ri>=0 && ri<ref.length) ? ref[ri] : (byte)'N'));
				
				if(m=='m' || m=='s'){
					if(!AminoAcid.isFullyDefined(c, colorspace) || !AminoAcid.isFullyDefined(r, colorspace)){
						match[mi]='N';
					}else{
						ok=(ok && c==r);
					}
					ci++;
					ri++;
				}else if(m=='D'){
					ri++;
				}else if(m=='I' || m=='X' || m=='Y'){
					ci++;
				}else if(m=='S'){
					ok=(ok && c!=r);
					ci++;
					ri++;
				}else if(m=='N'){
					ok=(ok && (!AminoAcid.isFullyDefined(c, colorspace) || !AminoAcid.isFullyDefined(r, colorspace)));
					ci++;
					ri++;
				}else{
					assert(false) : (char)m;
				}
			}
			
			if(call[0]<4){
				if(ref.length>400){
					System.err.println(toStringCS(Arrays.copyOfRange(ref, rstart, rstart+call.length))+" (ref)");
				}else{
					System.err.println(toStringCS(ref)+" (ref)");
				}
				System.err.println(toStringCS(call)+" (call)");
				System.err.println(new String(match));
			}else{
				if(ref.length>400){
					System.err.println(new String(Arrays.copyOfRange(ref, rstart, rstart+call.length))+" (ref)");
				}else{
					System.err.println(new String(ref)+" (ref)");
				}
				System.err.println(new String(call)+" (call)");
				System.err.println(new String(match));
			}
			

			
			if(THROW_EXCEPTION_ON_VERIFY_FAILURE){
				System.err.println("Fixed successfully?\t"+ok);
				throw new RuntimeException("Failed VerifyMatchString()");
			}
			
		}
		
		return ok;
	}
	
	//TODO: No-calls and no-ref are currently considered the same.
	/** When this is called, the match string should be plus-oriented */
	public ArrayList<Varlet> toVars(final Read read, final boolean CONDENSE, final boolean CONDENSE_SNPS, final boolean SPLIT_SUBS){
		assert(!read.colorspace());
		assert(read.match!=null);
		byte[] match=read.match;
		byte[] quality=read.quality;
		byte[] call=read.bases;
		
		if(quality==null){quality=Read.getFakeQuality(call.length);}
		
		assert(!read.colorspace()); //or else reverse-complement will mess things up
		
		assert(checkArray(call));
		
		int maxVars=0;
		
		byte last='m';
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='X' || b=='Y'){
//				assert(false) : read.toText(false);
				b=match[i]='I';
			} //TODO: Should not be needed, if reads are trimmed...
			
			if(b!='m' && b!=last){
				maxVars++;
			}
			last=b;
		}
		
		if(maxVars==0){return null;}
		
//		byte[] original=Arrays.copyOf(call, call.length);
		if(read.strand()==Gene.MINUS){
			AminoAcid.reverseComplementBasesInPlace(call);
			Tools.reverseInPlace(quality);
		}
		
		
		//assert(checkArray(call)) : 
//			"\n"+new String(original)+"\n"+new String(Tools.reverseAndCopy(call))+"\n"+
//			"\n"+Arrays.toString(original)+"\n"+Arrays.toString(Tools.reverseAndCopy(call))+"\n";
		
//		assert(false) : "TODO: ensure read is aligned with forward strand.";
		
		ArrayList<Varlet> vars=new ArrayList<Varlet>(maxVars);
		ChromosomeArray cha=Data.getChromosome(read.chrom);


//		assert(verifyMatchString(call, cha.array, match, read.start, true, read.colorspace())) : read.toText(false);
		boolean vms=false;
		try {
			vms=verifyMatchString(call, cha.array, match, read.start, true, read.colorspace());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			vms=false;
			System.err.println("in TranslateColorspace.toVars(), a read failed verification:\n"+read.toText(false)+"\n");
		}
		
		if(verbose){
			System.err.println("Making vars:");
			System.err.println(new String(call));
			System.err.println(cha.getString(read.start, read.stop));
			System.err.println(new String(match));

		}
		
		int readQuality;
		{
			int totalQual=0;
			int minQual=quality[0];
			for(int i=0; i<quality.length; i++){
				totalQual+=quality[i];
				minQual=Tools.min(minQual, quality[i]);
			}
			readQuality=(totalQual+2*minQual)/(read.bases.length+2);
		}
		final float expectedErrors=read.expectedErrors();
		
		last='m';
		int callPos=0;
		int refPos=read.start;
		
		//Make variations, then merge adjacent variations.
		for(int matchPos=0; matchPos<match.length; matchPos++){
			
			if(match[matchPos]=='N'){
				byte a=call[callPos];
				byte b=cha.get(refPos);
				if(a!='N' && b=='N'){match[matchPos]='R';}
			}
			
			final byte type=match[matchPos];
			
			if(type=='m'){
				callPos++;
				refPos++;
			}else{
				byte m;
				int nCount=0; //"no-call": "N" in read
				int rCount=0; //"no-ref":  "N" in ref but read is called
				int iCount=0;
				int dCount=0;
				int sCount=0;
				
				//call string
				StringBuilder cs=new StringBuilder(8);
				
				//ref string
				StringBuilder rs=new StringBuilder(8);
				
				final int mstart=matchPos;
				final int cstart=callPos;
				final int rstart=refPos;
				
				int qualSum=0;
				int qualMin=quality[callPos];
				
				while(matchPos<match.length && (m=match[matchPos])==type){
					
					//TODO: Not very good for deletions...
					qualSum+=quality[callPos];
					qualMin=Tools.min(qualMin, quality[callPos]);
					
					if(m=='I'){
						iCount++;
						cs.append((char)call[callPos]);
						callPos++;
					}else if(m=='D'){
						dCount++;
						rs.append((char)cha.get(refPos));
						refPos++;
					}else if(m=='S'){
						sCount++;
						cs.append((char)call[callPos]);
						rs.append((char)cha.get(refPos));
						assert(call[callPos]!='N');
						assert(cha.get(refPos)!='N');
						callPos++;
						refPos++;
						if(SPLIT_SUBS){
							matchPos++;break;//Forces all subs to be split
						}
					}else if(m=='N'){
						
						assert(call[callPos]=='N') : callPos+"\n"+new String(call)+"\n"+new String(match)+"\n"+cha.getString(read.start, read.stop)+"\n";
						nCount++;
//						cs.append((char)call[callPos]);
						cs.append('N');
						rs.append((char)cha.get(refPos));
						callPos++;
						refPos++;
						
						//This block corrects for a rare situation when both no-calls and no-refs are mixed in a single 'N' block.
						{
							int x=matchPos+1;
							if(x<match.length && match[x]=='N'){
								byte a=call[callPos];
								byte b=cha.get(refPos);
								if(a!='N' && b=='N'){match[x]='R';}
							}
						}
						
					}else if(m=='R'){
						assert(call[callPos]!='N');
						assert(cha.get(refPos)=='N');
						rCount++;
						cs.append((char)call[callPos]);
						rs.append((char)cha.get(refPos));
						callPos++;
						refPos++;
						matchPos++;break; //Output no-ref individually
					}else{
						System.err.println("Detected invalid decode for read "+read.numericID+":");
						System.err.println((char)m+"\n"+new String(rs)+"\n"+new String(cs)+"\n"+new String(match)+"\n"
						+new String(call)+"\n"+cha.getString(read.start, read.stop)+"\n"+read.toText(false)+"\n");
						return null;
//						assert(false) : (char)m+"\n"+new String(rs)+"\n"+new String(cs)+"\n"+new String(match)+"\n"
//						+new String(call)+"\n"+cha.getString(read.start, read.stop)+"\n"+read.toText(false)+"\n";
					}
					matchPos++;
				}
				matchPos--;
				
				int mstop=matchPos;
				
				int mlen=iCount+dCount+sCount+nCount+rCount;
				
				int clen=iCount+sCount+nCount+rCount;
				int rlen=dCount+sCount+nCount+rCount;
				
				Varlet v;
				
				callPos=cstart+clen;
				refPos=rstart+rlen;
				
				final int rstop=Tools.max(rstart, rstart+rlen-1);
				final int cstop=cstart+clen-1;
				
				final byte varType;
				
				if(rlen==0){
					varType=Variation.INS;
					if(verbose){System.err.println("Setting type INS: "+Variation.varTypeMap[varType]);}
				}else if(clen==0){varType=Variation.DEL;}
				else if(rCount>0){varType=Variation.NOREF;}
				else if(cs.charAt(0)=='N'){
					varType=Variation.NOCALL;
					if(verbose){System.err.println("Setting type NOCALL: "+Variation.varTypeMap[varType]);}
				}else if(mlen==1){varType=Variation.SNP;}
				else{varType=Variation.DELINS;}
				
				
				final int headDist, tailDist, endDist;
				{
					int cstart2=cstart, cstop2=cstop;
					if(varType==Variation.DEL){
						cstart2--;
						cstop2++;
					}
					
					assert(cstop2>=cstart2) : Variation.varTypeMap[varType]+", "+cstop2+", "+cstart2+", "+clen+
						"\n'"+cs+"', '"+rs+"'\n"+new String(match);
					assert(cstop2<call.length);
					
					if(read.strand()==Gene.PLUS){
						headDist=cstart2;
						tailDist=call.length-cstop2-1;
					}else{
						tailDist=cstart2;
						headDist=call.length-cstop2-1;
					}
					endDist=Tools.min(headDist, tailDist);
					assert(headDist>=0);
					assert(tailDist>=0);
				}
				
				
				int varQuality;
				if(varType==Variation.DEL){
					varQuality=((qualSum/mlen)+(qualMin))/2;
				}else{
					if(callPos<quality.length-1 && callPos>1){
						qualMin=Tools.min(quality[callPos-2], quality[callPos-2], quality[callPos-2], quality[callPos-2]);
						varQuality=(quality[callPos-2]+quality[callPos-1]+quality[callPos]+quality[callPos+1]+(qualMin))/5;
					}else if(callPos<quality.length && callPos>0){
						qualMin=Tools.min(quality[callPos-1], quality[callPos]);
						varQuality=qualMin;
					}else{
						varQuality=((qualSum/mlen)+(qualMin))/2;
					}
				}
				
				if(verbose){
					System.err.println("mlen="+mlen+", rlen="+rlen+", clen="+clen+", varType="+Variation.varTypeMap[varType]+"\n"+
							", cs="+cs+", nCount="+nCount+", rCount="+rCount+", iCount="+iCount+", dCount="+dCount+", sCount="+sCount);
				}
				
//				assert(read.mapScore>0) : read.toText(false);
				v=new Varlet(read.chrom, read.strand(), rstart, rstop, mstart, mstop, varType, rs.toString(), cs.toString(),
						 varQuality, readQuality, read.mapScore, read.errors, expectedErrors, (read.paired() ? 1 : 0), read.numericID, 
						 read.bases.length, read.mapLength, read.start, read.stop, read.copies, headDist, tailDist, endDist, 
						 read.pairnum());
				
//				if(v.varType==Variation.NOREF){System.err.print("R");}
				
				if(v.varType==Variation.SNP){
					if(v.call.equals(v.ref)){
						System.err.println("\n"+read.toText(false));
						System.err.println("\n"+v.toText());
						System.err.println("\n"+read.strand());
						System.err.println("\n");
						System.err.println(cha.getString(read.start, read.stop));
						System.err.println(new String(call));
						System.err.println(new String(match));
						System.err.println("\n");
						assert(false);
					}
					
				}
				
				vars.add(v);
			}
		}
		//assert(checkArray(call));
		
//		assert(read.numericID!=3448228) : CONDENSE+"\n"+vars;
		
//		boolean fail=false;
//		{
//			int nr=0;
//			for(Variation v : vars){
//				if(v.varType==Variation.NOREF){
//					nr++;
//					fail=nr>0;
//				}
//			}
//			System.err.print(" "+nr);
//		}
//		if(fail){verbose=true;}

//		if(read.numericID==3448228){verbose=true;}
		
		//Optionally, merge nearby variations
		if(CONDENSE && vars.size()>1){
			boolean condense=false;
			
			int mergeDistance=1; //  1 for adjacent, 2 for non-adjacent.
			
			for(int i=1; i<vars.size() && !condense; i++){
				Varlet v1=vars.get(i-1);
				Varlet v2=vars.get(i);
				assert(v1.matchStop<v2.matchStart);
				
				if(!v1.isNR_or_NC() && !v2.isNR_or_NC()){
					if(v1.endLoc>=v2.beginLoc){condense=true;} //To prevent overlapping variations
					else if(CONDENSE_SNPS || (v1.varType!=Variation.SNP && v2.varType!=Variation.SNP)){
						condense|=(v1.matchStop>=v2.matchStart-mergeDistance);
					}
				}

				if(verbose){
					System.err.println("Compared\n"+v1+"\nand\n"+v2+"\ncondense="+condense+"\n"+v1.matchStart+", "+v2.matchStart+", "+mergeDistance);
				}
			}
			
//			condense=false;
			if(condense){
				if(verbose){
					System.err.println("Condensing:");
					for(Varlet v : vars){
						System.err.println(v);
					}
				}
				ArrayList<Varlet> list2=new ArrayList<Varlet>(vars.size()-1);
				for(int i=vars.size()-2; i>=0; i--){
					Varlet prev=vars.get(i);
//					Varlet v=vars.get(i+1);
					Varlet v=vars.remove(i+1);
					
					
					boolean merge=(!v.isNR_or_NC() && !prev.isNR_or_NC() && (prev.matchStop>=v.matchStart-mergeDistance || prev.endLoc>=v.beginLoc));
					if(merge && !CONDENSE_SNPS && prev.endLoc<v.beginLoc){
						if(v.varType==Variation.SNP || prev.varType==Variation.SNP){
							merge=false;
						}
					}

					byte varType;
					
					if(merge){ //then merge.
						
//						if(v.varType==prev.varType){
//							varType=v.varType;
//						}else{
//							varType=Variation.DELINS;
//						}
						varType=Variation.DELINS;
						
						int midstart=prev.endLoc+1;
						int midstop=v.beginLoc-1;
						
						if(prev.varType==Variation.INS){midstart--;}
						
						String middle=(midstart>midstop ? "" : cha.getString(midstart, midstop));

						String cs=(prev.call==null ? "" : prev.call)+middle+(v.call==null ? "" : v.call);
						String rs=(prev.ref==null ? "" : prev.ref)+middle+(v.ref==null ? "" : v.ref);

						final int headDist=Tools.min(v.headDist, prev.headDist);
						final int tailDist=Tools.min(v.tailDist, prev.tailDist);
						final int endDist=Tools.min(v.endDist, prev.endDist);
						
						
						Varlet v2=new Varlet(read.chrom, read.strand(), prev.beginLoc, v.endLoc, prev.matchStart, v.matchStop, varType,
								rs, cs, (prev.avgVarQuality()+v.avgVarQuality())/2, readQuality, read.mapScore, read.errors, expectedErrors,
								(read.paired() ? 1 : 0), read.numericID, read.bases.length, read.mapLength,
								read.start, read.stop, read.copies, headDist, tailDist, endDist, read.pairnum());
						
						vars.remove(i); //prev
						vars.add(v2);
					}else{
						list2.add(v);
					}
				}
				assert(vars.size()==1);
				list2.add(vars.get(0));
				Collections.reverse(list2);
				vars=list2;
				
				if(verbose){
					System.err.println("Condensed:");
					for(Varlet v : vars){
						System.err.println(v);
					}
					System.err.println();
				}
			}
		}
		
//		{
//			int nr=0;
//			for(Variation v : vars){
//				if(v.varType==Variation.NOREF){
//					nr++;
//				}
//			}
//			System.err.println(" "+nr);
//		}
//		
//		assert(!fail);

//		assert(read.numericID!=3448228) : CONDENSE+"\n"+vars;
		
		//assert(checkArray(call));
		//Don't exit early and forget to undo this!
		if(read.strand()==Gene.MINUS){
			AminoAcid.reverseComplementBasesInPlace(call);
			Tools.reverseInPlace(quality);
		}
		//assert(checkArray(call));
		return vars;
	}
	

	protected MSA msaBS;
	protected MSA msaCS;

	public static boolean verbose=false;
	
	public static boolean DISCARD_NOCALLED_INSERTIONS=false;
	public static boolean THROW_EXCEPTION_ON_VERIFY_FAILURE=true; //Throws an exception when "verify match string" fails
	
}
