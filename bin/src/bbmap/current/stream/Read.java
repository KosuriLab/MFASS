package stream;

import java.util.ArrayList;
import java.util.Arrays;

import align2.GapTools;
import align2.QualityTools;
import align2.Shared;
import align2.Tools;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;

public final class Read implements Comparable<Read>, Cloneable{
	
	public static void main(String[] args){
		byte[] a=args[0].getBytes();
		System.out.println(new String(a));
		byte[] b=toShortMatchString(a);
		System.out.println(new String(b));
		byte[] c=toLongMatchString(b);
		System.out.println(new String(c));
		byte[] d=toLongMatchString(c);
		System.out.println(new String(d));
//		byte[] e=toShortMatchString(b);
//		System.out.println(new String(e));
		
	}
	
	public Read(byte[] bases_, byte[] quals_, long id_){
		this(bases_, -1, (byte)0, 0, 0, Long.toString(id_), quals_, false, id_);
	}
	
	public Read(byte[] s_, int chrom_, byte strand_, int start_, int stop_, long id_, byte[] quality_, boolean cs_){
		this(s_, chrom_, strand_, start_, stop_, Long.toString(id_), quality_, cs_, id_);
	}
	
	public Read(byte[][] fasta_, byte[][] qual_, boolean cs_, long numericID_){
		this(fasta_[1], 0, (byte)0, 0, 0, new String(fasta_[0]), qual_[1], cs_, numericID_);
	}
	
	public Read(byte[] s_, int chrom_, byte strand_, int start_, int stop_, String id_, byte[] quality_, boolean cs_, long numericID_){
		this(s_, chrom_, start_, stop_, id_, quality_, numericID_, (strand_|(cs_ ? COLORMASK : 0)));
		assert(strand_==0 || strand_==1);
		assert(start_<=stop_) : chrom_+", "+start_+", "+stop_+", "+numericID_;
	}
	
	public Read(byte[] s_, int chrom_, int start_, int stop_, String id_, byte[] quality_, long numericID_, int flags_){
		assert(quality_==null || quality_[0]<=80 || !FASTQ.DETECT_QUALITY) : "\n"+Arrays.toString(quality_)+
			"\n"+Arrays.toString(s_)+"\n"+numericID_+"\n"+id_+"\n"+FASTQ.ASCII_OFFSET;
		
		flags=flags_;
		byte[] basesOriginal=s_;
		byte[] qualityOriginal=quality_;
		if(qualityOriginal!=null && qualityOriginal.length==4){
			if(qualityOriginal[0]=='n' && 
					qualityOriginal[0]=='u' && 
					qualityOriginal[0]=='l' && 
					qualityOriginal[0]=='l'){
				qualityOriginal=null;
			}
		}
		
		if(NULLIFY_BROKEN_QUALITY && qualityOriginal!=null && qualityOriginal.length!=s_.length){
			qualityOriginal=null;
			setDiscarded(true);
		}
		assert(qualityOriginal==null || absdif(qualityOriginal.length, s_.length)<=2) : 
			"\nMismatch between length of bases and qualities for read "+numericID_+" (id="+id_+").\n"+
					"# qualities="+qualityOriginal.length+", # bases="+s_.length+"\n\n"+
					FASTQ.qualToString(qualityOriginal)+"\n"+new String(s_)+"\n";
		
		assert(basesOriginal.length<2 || basesOriginal[1]=='N' || basesOriginal[1]=='.' || basesOriginal[1]=='-' || colorspace()!=Character.isLetter(basesOriginal[1])) : 
			"\nAn input file appears to be misformatted.  The character with ASCII code "+basesOriginal[1]+" appeared where a base was expected.\n" +
					colorspace()+", "+Arrays.toString(basesOriginal);
		
		if(colorspace()){ //Trim tips of reads that have a primer base attached
			int x=0, y=basesOriginal.length;
			int xq=0, yq=qualityOriginal.length;
			if(basesOriginal[0]>3 && basesOriginal[0]!='N' && basesOriginal[0]!='.'){
				assert(basesOriginal[0]=='T' || basesOriginal[0]=='G') : "Just an assumption based on SOLiD - safe to disable.";
				assert(basesOriginal[1]<=3 || basesOriginal[1]=='N' || basesOriginal[1]=='.') : "Just an assumption based on SOLiD - safe to disable.";
				x+=2;
				xq+=1;
			}
			byte last=basesOriginal[basesOriginal.length-1];
			if(last>3 && last!='N' && last!='.'){
				//Might be 'G' too
				assert(last=='T' || last=='G') : "Just an assumption based on SOLiD - safe to disable.\n"+new String(basesOriginal);
//				assert(basesOriginal[basesOriginal.length-2]<=3) : "Just an assumption based on SOLiD - safe to disable.";
				y-=2;
				yq-=1;
			}
			
			if(x!=0 || y!=basesOriginal.length){
				bases=Arrays.copyOfRange(basesOriginal, x, y);
				assert(bases.length==y-x);
				assert(bases[bases.length-1]==basesOriginal[bases.length-1+x]);
				
				if(qualityOriginal!=null){
					quality=Arrays.copyOfRange(qualityOriginal, xq, yq);
					assert(bases.length==quality.length);
//					quality=new byte[yq-xq];
//					for(int i=0; i<quality.length; i++){
//						quality[i]=qualityOriginal[i+xq];
//					}
				}else{
					quality=null;
				}
//				bases=new byte[y-x];
//				for(int i=0; i<bases.length; i++){
//					bases[i]=basesOriginal[i+x];
//				}
			}else{
				bases=basesOriginal;
				quality=qualityOriginal;
			}
			
			if(quality!=null){
				
				for(int i=0; i<quality.length; i++){
					byte b=bases[i];
					byte q=quality[i];
					assert(b!='.' && b!='T');
					if(b>=0 && b<=3){
						if(q<MIN_CALLED_QUALITY){
//							assert(false) : b+", "+q;
							quality[i]=MIN_CALLED_QUALITY;
						}
					}else{
//						assert(false) : b+", "+q;
						quality[i]=0;
					}
				}
			}
		}else{
			
			bases=basesOriginal;
			quality=qualityOriginal;
			
			if(quality!=null){
				for(int i=0; i<quality.length; i++){
					byte b=bases[i];
					byte q=quality[i];
					if(AminoAcid.isFullyDefined(b)){
						if(q<MIN_CALLED_QUALITY){
//							assert(false) : (char)b+", "+q;
							quality[i]=MIN_CALLED_QUALITY;
						}
					}else{
//						assert(false) : (char)b+", "+q;
						quality[i]=0;
						if(b=='-' || b=='.'){bases[i]='N';}
					}
					if(TO_UPPER_CASE && b>90){bases[i]-=32;}
				}
			}else if(TO_UPPER_CASE){
				for(int i=0; i<bases.length; i++){if(bases[i]>90){bases[i]-=32;}}
			}
		}
		
		chrom=chrom_;
		start=start_;
		stop=stop_;
		
		id=id_;
		assert(quality==null || quality.length==bases.length) : "\n"+new String(bases)+"\n"+bases.length+"!="+quality.length+"\n"+
				Arrays.toString(quality)+"\n"+Arrays.toString(bases)+"\n";
		numericID=numericID_;
		
		mapLength=bases.length;
	}
	
	private static final int absdif(int a, int b){
		return a>b ? a-b : b-a;
	}
	
	
	
	/** Returns true if these reads are identical, allowing at most n no-calls and m mismatches of max quality q*/
	public boolean isDuplicateByBases(Read r, int nmax, int mmax, byte qmax, boolean banSameQualityMismatch){
		return isDuplicateByBases(r, nmax, mmax, qmax, false, false);
	}
	
	
	
	/** Returns true if these reads are identical, allowing at most n no-calls and m mismatches of max quality q*/
	public boolean isDuplicateByBases(Read r, int nmax, int mmax, byte qmax, boolean banSameQualityMismatch, boolean allowDifferentLength){
		int n=0, m=0;
		assert(r.bases.length==bases.length) : "Merging different-length reads is supported but seems to be not useful.";
		if(!allowDifferentLength && r.bases.length!=bases.length){return false;}
		int minLen=Tools.min(bases.length, r.bases.length);
		for(int i=0; i<minLen; i++){
			byte b1=bases[i];
			byte b2=r.bases[i];
			if(b1=='N' || b2=='N'){
				n++;
				if(n>nmax){return false;}
			}else if(b1!=b2){
				m++;
				if(m>mmax){return false;}
				if(quality[i]>qmax && r.quality[i]>qmax){return false;}
				if(banSameQualityMismatch && quality[i]==r.quality[i]){return false;}
			}
		}
		return true;
	}
	
	public boolean isDuplicateByMapping(Read r, boolean bothEnds, boolean checkAlignment){
		if(bases.length!=r.bases.length){
			return isDuplicateByMappingDifferentLength(r, bothEnds, checkAlignment);
		}
		assert(this!=r && mate!=r);
		assert(!bothEnds || bases.length==r.bases.length);
		if(!mapped() || !r.mapped()){return false;}
//		if(chrom==-1 && start==-1){return false;}
		if(chrom<1 && start<1){return false;}
		
//		if(chrom!=r.chrom || strand()!=r.strand() || start!=r.start){return false;}
////		if(mate==null && stop!=r.stop){return false;} //For unpaired reads, require both ends match
//		if(stop!=r.stop){return false;} //For unpaired reads, require both ends match
//		return true;
		
		if(chrom!=r.chrom || strand()!=r.strand()){return false;}
		if(bothEnds){
			if(start!=r.start || stop!=r.stop){return false;}
		}else{
			if(strand()==Gene.PLUS){
				if(start!=r.start){return false;}
			}else{
				if(stop!=r.stop){return false;}
			}
		}
		if(checkAlignment){
			if(perfect() && r.perfect()){return true;}
			if(match!=null && r.match!=null){
				if(match.length!=r.match.length){return false;}
				for(int i=0; i<match.length; i++){
					byte a=match[i];
					byte b=r.match[i];
					if(a!=b){
						if((a=='D') != (b=='D')){return false;}
						if((a=='I' || a=='X' || a=='Y') != (b=='I' || b=='X' || b=='Y')){return false;}
					}
				}
			}
		}
		return true;
	}
	
	public boolean isDuplicateByMappingDifferentLength(Read r, boolean bothEnds, boolean checkAlignment){
		assert(this!=r && mate!=r);
		assert(bases.length!=r.bases.length);
		if(bothEnds){return false;}
//		assert(!bothEnds || bases.length==r.bases.length);
		if(!mapped() || !r.mapped()){return false;}
//		if(chrom==-1 && start==-1){return false;}
		if(chrom<1 && start<1){return false;}
		
//		if(chrom!=r.chrom || strand()!=r.strand() || start!=r.start){return false;}
////		if(mate==null && stop!=r.stop){return false;} //For unpaired reads, require both ends match
//		if(stop!=r.stop){return false;} //For unpaired reads, require both ends match
//		return true;
		
		if(chrom!=r.chrom || strand()!=r.strand()){return false;}

		if(strand()==Gene.PLUS){
			if(start!=r.start){return false;}
		}else{
			if(stop!=r.stop){return false;}
		}
		
		if(checkAlignment){
			if(perfect() && r.perfect()){return true;}
			if(match!=null && r.match!=null){
				int minLen=Tools.min(match.length, r.match.length);
				for(int i=0; i<minLen; i++){
					byte a=match[i];
					byte b=r.match[i];
					if(a!=b){
						if((a=='D') != (b=='D')){return false;}
						if((a=='I' || a=='X' || a=='Y') != (b=='I' || b=='X' || b=='Y')){return false;}
					}
				}
			}
		}
		return true;
	}
	
	public void merge(Read r, boolean mergeVectors, boolean mergeN){mergePrivate(r, mergeVectors, mergeN, true);}
	
	private void mergePrivate(Read r, boolean mergeVectors, boolean mergeN, boolean mergeMate){
		assert(r!=this);
		assert(r!=this.mate);
		assert(r!=r.mate);
		assert(this!=this.mate);
		assert(r.mate==null || r.mate.mate==r);
		assert(this.mate==null || this.mate.mate==this);
		assert(r.mate==null || r.numericID==r.mate.numericID);
		assert(mate==null || numericID==mate.numericID);
		mergeN=(mergeN||mergeVectors);
		
		assert(r.bases.length==bases.length) : "Merging different-length reads is supported but seems to be not useful.";
		
		if((mergeN || mergeVectors) && bases.length<r.bases.length){
			int oldLenB=bases.length;
			start=Tools.min(start, r.start);
			stop=Tools.max(stop, r.stop);
			mapScore=Tools.max(mapScore, r.mapScore);
			mapLength=Tools.max(mapLength, r.mapLength);
			
			bases=Arrays.copyOfRange(bases, 0, r.bases.length);
			quality=Arrays.copyOfRange(quality, 0, r.quality.length);
			for(int i=oldLenB; i<bases.length; i++){
				bases[i]='N';
				quality[i]=0;
			}
			match=null;
			r.match=null;
		}
		
		assert(r.colorspace()==colorspace());
		copies+=r.copies;
		
		
//		if(numericID==11063941 || r.numericID==11063941 || numericID==8715632){
//			System.err.println("***************");
//			System.err.println(this.toText()+"\n");
//			System.err.println(r.toText()+"\n");
//			System.err.println(mergeVectors+", "+mergeN+", "+mergeMate+"\n");
//		}
		
		boolean pflag1=perfect();
		boolean pflag2=r.perfect();

		final int minLenB=Tools.min(bases.length, r.bases.length);
		
		if(mergeN){
			if(quality==null){
				for(int i=0; i<minLenB; i++){
					byte b=r.bases[i];
					if(bases[i]=='N' && b!='N'){bases[i]=b;}
				}
			}else{
				for(int i=0; i<minLenB; i++){
					final byte b1=bases[i];
					final byte b2=r.bases[i];
					final byte q1=Tools.max((byte)0, quality[i]);
					final byte q2=Tools.max((byte)0, r.quality[i]);
					if(b1==b2){
						if(b1=='N'){
							//do nothing
						}else if(mergeVectors){
							//merge qualities
							//						quality[i]=(byte) Tools.min(40, q1+q2);
							if(q1>=q2){
								quality[i]=(byte) Tools.min(48, q1+1+q2/4);
							}else{
								quality[i]=(byte) Tools.min(48, q2+1+q1/4);
							}
						}
					}else if(b1=='N'){
						bases[i]=b2;
						quality[i]=q2;
					}else if(b2=='N'){
						//do nothing
					}else if(mergeVectors){
						if(q1<1 && q2<1){
							//Special case - e.g. Illumina calls bases at 0 quality.
							//Possibly best to keep the matching allele if one matches the ref.
							//But for now, do nothing.
							//This was causing problems changing perfect match strings into imperfect matches.
						}else if(q1==q2){
							assert(b1!=b2);
							bases[i]='N';
							quality[i]=0;
						}else if(q1>q2){
							bases[i]=b1;
							quality[i]=(byte)(q1-q2/2);
						}else{
							bases[i]=b2;
							quality[i]=(byte)(q2-q1/2);
						}
						assert(quality[i]>=0 && quality[i]<=48);
					}
				}
			}
		}
		
		//TODO:
		//Note that the read may need to be realigned after merging, so the match string may be rendered incorrect.
		
		if(mergeN && match!=null){
			if(r.match==null){match=null;}
			else{
				if(match.length!=r.match.length){match=null;}
				else{
					boolean ok=true;
					for(int i=0; i<match.length && ok; i++){
						byte a=match[i], b=r.match[i];
						if(a!=b){
							if((a=='m' || a=='S') && b=='N'){
								//do nothing;
							}else if(a=='N' && (b=='m' || b=='S')){
								match[i]=b;
							}else{
								ok=false;
							}
						}
					}
					if(!ok){match=null;}
				}
			}
		}
		
		if(mergeMate && mate!=null){
			mate.mergePrivate(r.mate, mergeVectors, mergeN, false);
			assert(copies==mate.copies);
		}
		assert(copies>1);
		
		assert(r!=this);
		assert(r!=this.mate);
		assert(r!=r.mate);
		assert(this!=this.mate);
		assert(r.mate==null || r.mate.mate==r);
		assert(this.mate==null || this.mate.mate==this);
		assert(r.mate==null || r.numericID==r.mate.numericID);
		assert(mate==null || numericID==mate.numericID);
		
		
//		if(numericID==11063941 || r.numericID==11063941 || numericID==8715632){
//			System.err.println("\nAfter:\n");
//			System.err.println(this.toText()+"\n");
//			System.err.println(r.toText()+"\n");
//			System.err.println("***************");
//		}
	}
	
	public Read translateToColorspace(boolean appendT){
		assert(!colorspace());
		byte[] temp=(appendT ? AminoAcid.toColorspaceSimulated(bases) : AminoAcid.toColorspace(bases));
//		assert(false) : "\n"+new String(bases)+"\n->\n"+new String(temp);
		Read r=new Read(temp, chrom, start, stop-1, id, quality, numericID, (flags|COLORMASK));
		r.sites=sites;
		r.originalSite=originalSite;
		r.errors=errors;
		r.errorsCorrected=errorsCorrected;
		r.mapScore=mapScore;
		r.obj=obj;
		assert(false) : "TODO: Be sure to copy ALL fields, like flags, etc.";
		return r;
	}
	
	public String toString(){return toText(false).toString();}
	
	public StringBuilder toSites(){
		StringBuilder sb;
		if(numSites()==0){
			sb=new StringBuilder(2);
			sb.append('.');
		}else{
			sb=new StringBuilder(sites.size()*20);
			int appended=0;
			for(SiteScore ss : sites){
				if(appended>0){sb.append('\t');}
				if(ss!=null){
					sb.append(ss.toText());
					appended++;
				}
			}
			if(appended==0){sb.append('.');}
		}
		return sb;
	}
	
	public ByteBuilder toSitesB(ByteBuilder sb){
		if(numSites()==0){
			if(sb==null){sb=new ByteBuilder(2);}
			sb.append('.');
		}else{
			if(sb==null){sb=new ByteBuilder(sites.size()*20);}
			int appended=0;
			for(SiteScore ss : sites){
				if(appended>0){sb.append('\t');}
				if(ss!=null){
					ss.toBytes(sb);
					appended++;
				}
			}
			if(appended==0){sb.append('.');}
		}
		return sb;
	}
	
	public StringBuilder toInfo(){
		if(obj==null){return new StringBuilder();}
		if(obj.getClass()==StringBuilder.class){return (StringBuilder)obj;}
		return new StringBuilder(obj.toString());
	}
	
	public ByteBuilder toInfoB(){
		if(obj==null){return new ByteBuilder();}
		if(obj.getClass()==ByteBuilder.class){return (ByteBuilder)obj;}
		return new ByteBuilder(obj.toString());
	}
	
	public StringBuilder toFastq(){
		return FASTQ.toFASTQ(this, (StringBuilder)null);
	}
	
	public ByteBuilder toFastq(ByteBuilder bb){
		return FASTQ.toFASTQ(this, bb);
	}

	public StringBuilder toFasta(){return toFasta(FastaReadInputStream.DEFAULT_WRAP);}
	public ByteBuilder toFasta(ByteBuilder bb){return toFasta(FastaReadInputStream.DEFAULT_WRAP, bb);}
	
	public StringBuilder toFasta(int wrap){
		if(wrap<1){wrap=Integer.MAX_VALUE;}
		int len=(id==null ? Tools.stringLength(numericID) : id.length())+(bases==null ? 0 : bases.length+bases.length/wrap)+5;
		StringBuilder sb=new StringBuilder(len);
		sb.append('>');
		if(id==null){sb.append(numericID);}
		else{sb.append(id);}
		sb.append('\n');
		if(bases!=null){
//			for(int i=0, j=0; i<bases.length; i++, j++){
//				if(j==wrap){
//					sb.append('\n');
//					j=0;
//				}
//				sb.append((char)bases[i]);
//			}
			final char[] buffer=Shared.getTLCB(Tools.min(wrap, bases.length));
			int j=0;
			for(int i=0; i<bases.length; i++, j++){
				if(j==wrap){
					sb.append(buffer, 0, j);
					sb.append('\n');
					j=0;
				}
				buffer[j]=(char)bases[i];
			}
			if(j>0){sb.append(buffer, 0, j);}
		}
		return sb;
	}
	
	public ByteBuilder toFasta(int wrap, ByteBuilder sb){
		if(wrap<1){wrap=Integer.MAX_VALUE;}
		int len=(id==null ? Tools.stringLength(numericID) : id.length())+(bases==null ? 0 : bases.length+bases.length/wrap)+5;
		if(sb==null){sb=new ByteBuilder(len+1);}
		sb.append('>');
		if(id==null){sb.append(numericID);}
		else{sb.append(id);}
		if(bases!=null){
			int pos=0;
			while(pos<bases.length-wrap){
				sb.append('\n');
				sb.append(bases, pos, wrap);
				pos+=wrap;
			}
			if(pos<bases.length){
				sb.append('\n');
				sb.append(bases, pos, bases.length-pos);
			}
		}
		return sb;
	}
	
	public StringBuilder toSam(){
		return new SamLine(this, pairnum()).toText();
	}
	
	public static CharSequence header(){

		StringBuilder sb=new StringBuilder();
		sb.append("id");
		sb.append('\t');
		sb.append("numericID");
		sb.append('\t');
		sb.append("chrom");
		sb.append('\t');
		sb.append("strand");
		sb.append('\t');
		sb.append("start");
		sb.append('\t');
		sb.append("stop");
		sb.append('\t');
		
//		sb.append("colorspace");
//		sb.append('\t');
//		sb.append("paired");
//		sb.append('\t');

		sb.append("flags");
		sb.append('\t');
		
		sb.append("copies");
		sb.append('\t');
		
		sb.append("errors,fixed");
		sb.append('\t');
		sb.append("mapScore");
		sb.append('\t');
		sb.append("length");
		sb.append('\t');
		
		sb.append("bases");
		sb.append('\t');
		sb.append("quality");
		sb.append('\t');
		
		sb.append("insert");
		sb.append('\t');
		{
			//These are not really necessary...
			sb.append("avgQual");
			sb.append('\t');
		}
		
		sb.append("match");
		sb.append('\t');
		sb.append("SiteScores: "+SiteScore.header());
		return sb;
	}
	
	public StringBuilder toText(boolean okToCompressMatch){
		
		final byte[] oldmatch=match;
		final boolean oldshortmatch=this.shortmatch();
		if(COMPRESS_MATCH_BEFORE_WRITING && !shortmatch() && okToCompressMatch){
			match=toShortMatchString(match);
			setShortMatch(true);
		}
		
		StringBuilder sb=new StringBuilder();
		sb.append(id);
		sb.append('\t');
		sb.append(numericID);
		sb.append('\t');
		sb.append(chrom);
		sb.append('\t');
		sb.append(Gene.strandCodes[strand()]);
		sb.append('\t');
		sb.append(start);
		sb.append('\t');
		sb.append(stop);
		sb.append('\t');
		
//		sb.append(colorspace() ? 1 : 0);
//		sb.append('\t');
//		
//		sb.append(paired() ? 1 : 0);
//		sb.append('\t');
		
		for(int i=maskArray.length-1; i>=0; i--){
			sb.append(flagToNumber(maskArray[i]));
		}
		sb.append('\t');
		
		sb.append(copies);
		sb.append('\t');

		sb.append(errors);
		if(errorsCorrected>0){
			sb.append(',');
			sb.append(errorsCorrected);
		}
		sb.append('\t');
		sb.append(mapScore);
		sb.append('\t');
		sb.append(mapLength);
		sb.append('\t');
		
		if(bases==null){sb.append('.');}
		else{
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				if(b<4){
					assert(b>=0);
					b=(byte) (b+'0');
				}
				sb.append((char)b);
			}
		}
		sb.append('\t');
		
		int qualSum=0;
		int qualMin=99999;
		
		if(quality==null){
			sb.append('.');
		}else{
			for(int i=0; i<quality.length; i++){
				byte q=quality[i];
				qualSum+=q;
				qualMin=Tools.min(q, qualMin);
				q=(byte) (q+ASCII_OFFSET);
				sb.append((char)q);
			}
		}
		sb.append('\t');
		
		if(insert<1){sb.append('.');}else{sb.append(insert);};
		sb.append('\t');
		
		if(quality==null){
			sb.append('.');
			sb.append('\t');
		}else{
			//These are not really necessary...
			sb.append(qualSum/quality.length);
			sb.append('\t');
		}
		
		if(match==null){sb.append('.');}
		else{for(byte b : match){sb.append((char)b);}}
		sb.append('\t');
		
		if(gaps==null){
			sb.append('.');
		}else{
			for(int i=0; i<gaps.length; i++){
				if(i>0){sb.append('~');}
				sb.append(gaps[i]);
			}
		}
		
		if(sites!=null && sites.size()>0){
			
			assert(absdif(start, stop)<3000 || (gaps==null) == (sites.get(0).gaps==null)) : 
				"\n"+this.numericID+"\n"+Arrays.toString(gaps)+"\n"+sites.toString()+"\n";
			
			for(SiteScore ss : sites){
				sb.append('\t');
				sb.append(ss==null ? "null" : ss.toText());
			}
		}
		
		if(originalSite!=null){
			sb.append('\t');
			sb.append('*');
			sb.append(originalSite.toText());
		}
		
		match=oldmatch;
		setShortMatch(oldshortmatch);
		
		return sb;
	}
	
	public ByteBuilder toText(boolean okToCompressMatch, ByteBuilder bb){
		
		final byte[] oldmatch=match;
		final boolean oldshortmatch=this.shortmatch();
		if(COMPRESS_MATCH_BEFORE_WRITING && !shortmatch() && okToCompressMatch){
			match=toShortMatchString(match);
			setShortMatch(true);
		}
		
		if(bb==null){bb=new ByteBuilder();}
		bb.append(id);
		bb.append('\t');
		bb.append(numericID);
		bb.append('\t');
		bb.append(chrom);
		bb.append('\t');
		bb.append(Gene.strandCodes2[strand()]);
		bb.append('\t');
		bb.append(start);
		bb.append('\t');
		bb.append(stop);
		bb.append('\t');
		
//		sb.append(colorspace() ? 1 : 0);
//		sb.append('\t');
//		
//		sb.append(paired() ? 1 : 0);
//		sb.append('\t');
		
		for(int i=maskArray.length-1; i>=0; i--){
			bb.append(flagToNumber(maskArray[i]));
		}
		bb.append('\t');
		
		bb.append(copies);
		bb.append('\t');

		bb.append(errors);
		if(errorsCorrected>0){
			bb.append(',');
			bb.append(errorsCorrected);
		}
		bb.append('\t');
		bb.append(mapScore);
		bb.append('\t');
		bb.append(mapLength);
		bb.append('\t');
		
		if(bases==null){bb.append('.');}
		else{
			if(colorspace()){
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					if(b<4){
						assert(b>=0);
						b=(byte) (b+'0');
					}
					bb.append((char)b);
				}
			}else{
				bb.append(bases);
			}
		}
		bb.append('\t');
		
//		int qualSum=0;
//		int qualMin=99999;
		
		if(quality==null){
			bb.append('.');
		}else{
			bb.ensureExtra(quality.length);
			for(int i=0, j=bb.length; i<quality.length; i++, j++){
				byte q=quality[i];
				bb.array[j]=(byte)(q+ASCII_OFFSET);
//				qualSum+=q;
//				qualMin=Tools.min(q, qualMin);
			}
			bb.length+=quality.length;
		}
		bb.append('\t');
		
		if(insert<1){bb.append('.');}else{bb.append(insert);};
		bb.append('\t');
		
		if(true || quality==null){
			bb.append('.');
			bb.append('\t');
		}else{
//			//These are not really necessary...
//			sb.append(qualSum/quality.length);
//			sb.append('\t');
		}
		
		if(match==null){bb.append('.');}
		else{bb.append(match);}
		bb.append('\t');
		
		if(gaps==null){
			bb.append('.');
		}else{
			for(int i=0; i<gaps.length; i++){
				if(i>0){bb.append('~');}
				bb.append(gaps[i]);
			}
		}
		
		if(sites!=null && sites.size()>0){
			
			assert(absdif(start, stop)<3000 || (gaps==null) == (sites.get(0).gaps==null)) : 
				"\n"+this.numericID+"\n"+Arrays.toString(gaps)+"\n"+sites.toString()+"\n";
			
			for(SiteScore ss : sites){
				bb.append('\t');
				if(ss==null){
					bb.append((byte[])null);
				}else{
					ss.toBytes(bb);
				}
				bb.append(ss==null ? "null" : ss.toText());
			}
		}
		
		if(originalSite!=null){
			bb.append('\t');
			bb.append('*');
			originalSite.toBytes(bb);
		}
		
		match=oldmatch;
		setShortMatch(oldshortmatch);
		
		return bb;
	}
	
	public static Read fromText(String line){
		if(line.length()==1 && line.charAt(0)=='.'){return null;}
		
		String[] split=line.split("\t");
		
		if(split.length<17){
			throw new RuntimeException("Error parsing read from text.\n" +
					"This may be caused be attempting to parse the wrong format.\n" +
					"Please ensure that the file extension is correct:\n" +
					"\tFASTQ should end in .fastq or .fq\n" +
					"\tFASTA should end in .fasta or .fa, .fas, .fna, .ffn, .frn, .seq, .fsa\n" +
					"\tSAM should end in .sam\n" +
					"\tNative format should end in .txt or .bread\n" +
					"If a file is compressed, there must be a compression extension after the format extension:\n" +
					"\tgzipped files should end in .gz or .gzip\n" +
					"\tzipped files should end in .zip and have only 1 file per archive\n" +
					"\tbz2 files should end in .bz2");
		}
		
		final String id=new String(split[0]);
		long numericID=Long.parseLong(split[1]);
		int chrom=Byte.parseByte(split[2]);
//		byte strand=Byte.parseByte(split[3]);
		int start=Integer.parseInt(split[4]);
		int stop=Integer.parseInt(split[5]);
		
//		boolean cs=(Integer.parseInt(split[6])==1);
//		boolean paired=(Integer.parseInt(split[7])==1);
		
		int flags=Integer.parseInt(split[6], 2);
		boolean cs=((flags&COLORMASK)!=0);
		
		int copies=Integer.parseInt(split[7]);

		int errors;
		int errorsCorrected;
		if(split[8].indexOf(',')>=0){
			String[] estring=split[8].split(",");
			errors=Integer.parseInt(estring[0]);
			errorsCorrected=Integer.parseInt(estring[1]);
		}else{
			errors=Integer.parseInt(split[8]);
			errorsCorrected=0;
		}
		
		int mapScore=Integer.parseInt(split[9]);
		int mapLen=Integer.parseInt(split[10]);
		
		byte[] basesOriginal=split[11].getBytes();
		byte[] qualityOriginal=(split[12].equals(".") ? null : split[12].getBytes());
		
		if(cs){
			for(int i=0; i<basesOriginal.length; i++){
				byte b=basesOriginal[i];
				if(b>='0' && b<='3'){
					b=(byte) (b-'0');
				}
				basesOriginal[i]=b;
			}
		}
		
		if(qualityOriginal!=null){
			for(int i=0; i<qualityOriginal.length; i++){
				byte b=qualityOriginal[i];
				b=(byte) (b-ASCII_OFFSET);
				assert(b>=-1) : b;
				qualityOriginal[i]=b;
			}
		}
		
		int insert=-1;
		if(!split[13].equals(".")){insert=Integer.parseInt(split[13]);}
		
		byte[] match=null;
		if(!split[15].equals(".")){match=split[15].getBytes();}
		int[] gaps=null;
		if(!split[16].equals(".")){
			
			String[] gstring=split[16].split("~");
			gaps=new int[gstring.length];
			for(int i=0; i<gstring.length; i++){
				gaps[i]=Integer.parseInt(gstring[i]);
			}
		}
		
//		assert(false) : split[16];
		
		Read r=new Read(basesOriginal, chrom, start, stop, id, qualityOriginal, numericID, flags);
		r.match=match;
		r.errors=errors;
		r.errorsCorrected=errorsCorrected;
		r.mapScore=mapScore;
		r.copies=copies;
		r.mapLength=mapLen;
		r.gaps=gaps;
		r.insert=insert;
		
		int firstScore=(ADD_BEST_SITE_TO_LIST_FROM_TEXT) ? 17 : 18;
		
		int scores=split.length-firstScore;
		
		int mSites=0;
		for(int i=firstScore; i<split.length; i++){
			if(split[i].charAt(0)!='*'){mSites++;}
		}
		if(mSites>0){r.sites=new ArrayList<SiteScore>(mSites);}
		for(int i=firstScore; i<split.length; i++){
			SiteScore ss=SiteScore.fromText(split[i]);
			if(split[i].charAt(0)=='*'){r.originalSite=ss;}
			else{r.sites.add(ss);}
		}
		
		if(DECOMPRESS_MATCH_ON_LOAD && r.shortmatch()){
			r.match=toLongMatchString(match);
			r.setShortMatch(false);
		}

		assert(r.numSites()==0 || absdif(r.start, r.stop)<3000 || (r.gaps==null) == (r.topSite().gaps==null)) : 
			"\n"+r.numericID+", "+r.chrom+", "+r.strand()+", "+r.start+", "+r.stop+", "+Arrays.toString(r.gaps)+"\n"+r.sites+"\n"+line+"\n";
		
		return r;
	}

//	/** Reverses the read.  Mainly for testing. 
//	 * Seems to be unused. */
//	@Deprecated
//	protected void reverse() {
//		Tools.reverseInPlace(bases);
//		Tools.reverseInPlace(quality);
//		Tools.reverseInPlace(match);
//	}

	/** Reverse-complements the read. */
	public void reverseComplement() {
		AminoAcid.reverseComplementBasesInPlace(bases);
		Tools.reverseInPlace(quality);
		setStrand(strand()^1);
	}
	
	@Override
	public int compareTo(Read o) {
		if(chrom!=o.chrom){return chrom-o.chrom;}
		if(start!=o.start){return start-o.start;}
		if(stop!=o.stop){return stop-o.stop;}
		if(strand()!=o.strand()){return strand()-o.strand();}
		return 0;
	}
	
	public SiteScore toSite(){
		assert(start<=stop) : this.toText(false);
		SiteScore ss=new SiteScore(chrom, strand(), start, stop, 0, 0, rescued(), perfect());
		ss.slowScore=mapScore;
		ss.gaps=gaps;
		ss.match=match;
		originalSite=ss;
		return ss;
	}
	
	public SiteScore topSite(){
		final SiteScore ss=(sites==null || sites.isEmpty()) ? null : sites.get(0);
		assert(sites==null || sites.isEmpty() || ss!=null) : "Top site is null for read "+this;
		return ss;
	}
	
	public int numSites(){
		return (sites==null ? 0 : sites.size());
	}
	
	public SiteScore makeOriginalSite(){
		originalSite=toSite();
		return originalSite;
	}
	
	public void setFromSite(SiteScore ss){
		assert(ss!=null);
		chrom=ss.chrom;
		setStrand(ss.strand);
		start=ss.start;
		stop=ss.stop;
		mapScore=ss.slowScore;
		setRescued(ss.rescued);
		gaps=ss.gaps;
		setPerfect(ss.perfect);
		
		match=ss.match;
		
		if(gaps!=null){
			gaps=ss.gaps=GapTools.fixGaps(start, stop, gaps, Shared.MINGAP);
//			gaps[0]=Tools.min(gaps[0], start);
//			gaps[gaps.length-1]=Tools.max(gaps[gaps.length-1], stop);
		}
	}
	
//	public static int[] fixGaps(int a, int b, int[] gaps, int minGap){
////		System.err.println("fixGaps input: "+a+", "+b+", "+Arrays.toString(gaps)+", "+minGap);
//		int[] r=GapTools.fixGaps(a, b, gaps, minGap);
////		System.err.println("fixGaps output: "+Arrays.toString(r));
//		return r;
//	}

	public void setFromOriginalSite(){
		setFromSite(originalSite);
	}
	public void setFromTopSite(){
		final SiteScore ss=topSite();
		if(ss==null){
			clearSite();
			setMapped(false);
			return;
		}
		setMapped(true);
		setFromSite(ss);
	}
	
	public void setFromTopSite(boolean randomIfAmbiguous, boolean primary, int maxPairDist){
		final SiteScore ss0=topSite();
		if(ss0==null){
			clearSite();
			setMapped(false);
			return;
		}
		setMapped(true);
		
		if(sites.size()==1 || !randomIfAmbiguous || !ambiguous()){
			setFromSite(ss0);
			return;
		}
		
		if(primary || mate==null || !mate.mapped() || !mate.paired()){
			int count=1;
			for(int i=1; i<sites.size(); i++){
				SiteScore ss=sites.get(i);
				if(ss.score<ss0.score || (ss0.perfect && !ss.perfect) || (ss0.semiperfect && !ss.semiperfect)){break;}
				count++;
			}

			int x=(int)(numericID%count);
			if(x>0){
				SiteScore ss=sites.get(x);
				sites.set(0, ss);
				sites.set(x, ss0);
			}
			setFromSite(sites.get(0));
			return;
		}
		
//		assert(false) : "TODO: Proper strand orientation, and more.";
		//TODO: Also, this code appears to sometimes duplicate sitescores(?)
//		for(int i=0; i<list.size(); i++){
//			SiteScore ss=list.get(i);
//			if(ss.chrom==mate.chrom && Tools.min(Tools.absdifUnsigned(ss.start, mate.stop), Tools.absdifUnsigned(ss.stop, mate.start))<=maxPairDist){
//				list.set(0, ss);
//				list.set(i, ss0);
//				setFromSite(ss);
//				return;
//			}
//		}
		
		//If unsuccessful, recur unpaired.
		
		this.setPaired(false);
		mate.setPaired(false);
		setFromTopSite(randomIfAmbiguous, true, maxPairDist);
	}
	
	public void clearPairMapping(){
		clearMapping();
		if(mate!=null){mate.clearMapping();}
	}
	
	public void clearMapping(){
		clearSite();
		match=null;
		sites=null;
		setMapped(false);
		setPaired(false);
	}
	
	public void clearSite(){
		chrom=-1;
		setStrand(0);
		start=-1;
		stop=-1;
//		errors=0;
		mapScore=0;
		gaps=null;
	}


	public void clearAnswers(boolean clearMate) {
//		assert(mate==null || (pairnum()==0 && mate.pairnum()==1)) : pairnum()+", "+mate.pairnum();
		clearSite();
		match=null;
		sites=null;
		flags=(flags&(SYNTHMASK|COLORMASK|PAIRNUMMASK|SWAPMASK));
		if(clearMate && mate!=null){
			mate.clearSite();
			mate.match=null;
			mate.sites=null;
			mate.flags=(mate.flags&(SYNTHMASK|COLORMASK|PAIRNUMMASK|SWAPMASK));
		}
//		assert(mate==null || (pairnum()==0 && mate.pairnum()==1)) : pairnum()+", "+mate.pairnum();
	}
	
	
	public boolean isBadPair(boolean requireCorrectStrands, boolean sameStrandPairs, int maxdist){
		if(mate==null || paired()){return false;}
		if(!mapped() || !mate.mapped()){return false;}
		if(chrom!=mate.chrom){return true;}
		
		{
			int inner;
			if(start<=mate.start){inner=mate.start-stop;}
			else{inner=start-mate.stop;}
			if(inner>maxdist){return true;}
		}
//		if(absdif(start, mate.start)>maxdist){return true;}
		if(requireCorrectStrands){
			if((strand()==mate.strand())!=sameStrandPairs){return true;}
		}
		if(!sameStrandPairs){
			if(strand()==Gene.PLUS && mate.strand()==Gene.MINUS){
				if(start>=mate.stop){return true;}
			}else if(strand()==Gene.MINUS && mate.strand()==Gene.PLUS){
				if(mate.start>=stop){return true;}
			}
		}
		return false;
	}
	
	public int countMismatches(){
		assert(match!=null);
		int x=0;
		for(byte b : match){
			if(b=='S'){x++;}
		}
		return x;
	}
	

	
	/**
	 * @param match string
	 * @return Total number of match, sub, del, ins, or clip symbols
	 */
	public static final int[] matchToMsdicn(byte[] match) {
		if(match==null || match.length<1){return null;}
		int[] msdicn=new int[6];
		
		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						msdicn[0]+=current;
					}else if(mode=='S'){
						msdicn[1]+=current;
					}else if(mode=='D'){
						msdicn[2]+=current;
					}else if(mode=='I'){
						msdicn[3]+=current;
					}else if(mode=='C' || mode=='X' || mode=='Y'){
						msdicn[4]+=current;
					}else if(mode=='N' || mode=='R'){
						msdicn[5]+=current;
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Character.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				msdicn[0]+=current;
			}else if(mode=='S'){
				msdicn[1]+=current;
			}else if(mode=='D'){
				msdicn[2]+=current;
			}else if(mode=='I'){
				msdicn[3]+=current;
			}else if(mode=='C' || mode=='X' || mode=='Y'){
				msdicn[4]+=current;
			}else if(mode=='N' || mode=='R'){
				msdicn[5]+=current;
			}
		}
		return msdicn;
	}
	

	
	/**
	 * Handles short or long mode.
	 * @param match string
	 * @return Total number of match, sub, del, ins, or clip symbols
	 */
	public static final float identity(byte[] match) {
//		assert(false) : new String(match);
		if(match==null || match.length<1){return 0;}
		
		int good=0, bad=0, n=0;
		
		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						good+=current;
//						System.out.println("G: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}else if(mode=='D'){
						int x;
						x=(int)Math.ceil(Math.sqrt(current));
						//x=(int)Math.ceil(Tools.log2(current));
						//May need special handling....
						bad+=(Tools.min(x, current));
						
//						System.out.println("D: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad+", x="+x);
					}else if(mode=='R' || mode=='N'){
						n+=current;
					}else if(mode=='C'){
						//Do nothing
						//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
					}else if(mode!='0'){
						assert(mode=='S' || mode=='I' || mode=='X' || mode=='Y') : (char)mode;
						bad+=current;
//						System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Character.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				good+=current;
			}else if(mode=='R' || mode=='N'){
				n+=current;
			}else if(mode=='C'){
				//Do nothing
				//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
			}else if(mode!='0'){
				assert(mode=='S' || mode=='I' || mode=='X' || mode=='Y') : (char)mode;
				bad+=current;
//				System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
			}
		}
		
		
		n=(n+3)/4;
		good+=n;
		bad+=n;
		float r=good/(float)Tools.max(good+bad, 1);
//		assert(false) : new String(match)+"\nmode='"+(char)mode+"', current="+current+", good="+good+", bad="+bad;

//		System.out.println("match="+new String(match)+"\ngood="+good+", bad="+bad+", r="+r);
//		System.out.println(Arrays.toString(matchToMsdicn(match)));
		
		return r;
	}
	
	
	public int avgQuality(){
		if(quality==null){return 40;}
		assert(quality!=null);
		int x=0;
		for(byte b : quality){
			x+=(b<0 ? 0 : b);
		}
		return x/quality.length;
	}
	
	public int avgQualityFirstNBases(int n){
		if(quality==null || n<1){return 30;}
		assert(quality!=null);
		int x=0;
		if(n>quality.length){return 0;}
		for(int i=0; i<n; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/n;
	}
	
	public int avgQualityLastNBases(int n){
		if(quality==null || n<1){return 30;}
		assert(quality!=null);
		int x=0;
		if(n>quality.length){return 0;}
		for(int i=bases.length-n; i<bases.length; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/n;
	}
	
	public byte minQualityFirstNBases(int n){
		if(quality==null || n<1){return 30;}
		assert(quality!=null && n>0);
		if(n>quality.length){return 0;}
		byte x=quality[0];
		for(int i=1; i<n; i++){
			byte b=quality[i];
			if(b<x){x=b;}
		}
		return x;
	}
	
	public byte minQualityLastNBases(int n){
		if(quality==null || n<1){return 30;}
		assert(quality!=null && n>0);
		if(n>quality.length){return 0;}
		byte x=quality[bases.length-n];
		for(int i=bases.length-n; i<bases.length; i++){
			byte b=quality[i];
			if(b<x){x=b;}
		}
		return x;
	}
	
	public boolean containsNonM(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m'){return true;}
		}
		return false;
	}
	
	public boolean containsNonNM(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='N'){return true;}
		}
		return false;
	}
	
	public boolean containsNonNMXY(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='N' && b!='X' && b!='Y'){return true;}
		}
		return false;
	}
	
	public boolean containsSDI(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b=='S' || b=='s' || b=='D' || b=='I'){return true;}
		}
		return false;
	}
	
	public boolean containsNonNMS(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='s' && b!='N' && b!='S'){return true;}
		}
		return false;
	}
	
	public boolean containsConsecutiveS(int num){
		assert(match!=null && valid() && !shortmatch());
		int cnt=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b=='S'){
				cnt++;
				if(cnt>=num){return true;}
			}else{
				cnt=0;
			}
		}
		return false;
	}
	
	public boolean containsIndels(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='I' || b=='D' || b=='X' || b=='Y'){return true;}
		}
		return false;
	}
	
	public boolean containsInMatch(char c){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b==c){return true;}
		}
		return false;
	}
	
	public boolean containsNocalls(){
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(b=='N'){return true;}
		}
		return false;
	}
	
	public int countNocalls(){
		int n=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(b=='N'){n++;}
		}
		return n;
	}
	
	public boolean containsUndefined(){
		for(byte b : bases){
			if(AminoAcid.baseToNumber[b]<0){return true;}
		}
		return false;
	}
	
	public int countUndefined(){
		int n=0;
		for(byte b : bases){
			if(AminoAcid.baseToNumber[b]<0){n++;}
		}
		return n;
	}
	
	public boolean containsXY(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='X' || b=='Y'){return true;}
		}
		return false;
	}
	
	public boolean containsXY2(){
		if(match==null || match.length<1){return false;}
		boolean b=(match[0]=='X' || match[match.length-1]=='Y');
		assert(!valid() || b==containsXY());
		return b;
	}
	
	/** Replaces 'B' in match string with 'S', 'm', or 'N' */
	public boolean fixMatchB(){
		assert(match!=null);
		final ChromosomeArray ca;
		if(Data.GENOME_BUILD>=0){
			ca=Data.getChromosome(chrom);
		}else{
			ca=null;
		}
		boolean originallyShort=shortmatch();
		if(originallyShort){match=toLongMatchString(match);}
		int mloc=0, cloc=0, rloc=start;
		for(; mloc<match.length; mloc++){
			byte m=match[mloc];
			
			if(m=='B'){
				byte r=(ca==null ? (byte)'?' : ca.get(rloc));
				byte c=bases[cloc];
				if(r=='N' || c=='N'){
					match[mloc]='N';
				}else if(r==c || Character.toUpperCase(r)==Character.toUpperCase(c)){
					match[mloc]='m';
				}else{
					if(ca==null){
						if(originallyShort){match=toShortMatchString(match);}
						return false;
					}
					match[mloc]='S';
				}
				cloc++;
				rloc++;
			}else if(m=='m' || m=='S' || m=='N' || m=='s' || m=='C'){
				cloc++;
				rloc++;
			}else if(m=='D'){
				rloc++;
			}else if(m=='I' || m=='X' || m=='Y'){
				cloc++;
			}
		}
		if(originallyShort){match=toShortMatchString(match);}
		return true;
	}
	
	public float expectedErrors(){
		if(quality==null){return 0;}
		final float[] array=QualityTools.PROB_ERROR;
		assert(array[0]>0 && array[0]<1);
		float sum=0;
		for(int i=0; i<quality.length; i++){
			byte b=bases[i];
			byte q=quality[i];
			if(b=='N'){
				assert(q==0);
//				sum+=1;
			}else{
				sum+=array[q];
			}
		}
		return sum;
	}

	public int estimateErrors() {
		if(quality==null){return 0;}
		assert(match!=null) : this.toText(false);
		
		int count=0;
		for(int ci=0, mi=0; ci<bases.length && mi<match.length; mi++){
			
//			byte b=bases[ci];
			byte q=quality[ci];
			byte m=match[mi];
			if(m=='m' || m=='s' || m=='N'){
				ci++;
			}else if(m=='X' || m=='Y'){
				ci++;
				count++;
			}else if(m=='I'){
				ci++;
			}else if(m=='D'){
				
			}else if(m=='S'){
				ci++;
				if(q<19){
					count++;
				}
			}
			
		}
		return count;
	}
	
	/** {M, S, D, I, N} */
	public int[] countErrors() {
		assert(match!=null) : this.toText(false);
		int m=0;
		int s=0;
		int d=0;
		int i=0;
		int n=0;
		
		for(byte b : match){
			if(b=='m'){
				m++;
			}else if(b=='N' || b=='C'){
				n++;
			}else if(b=='X' || m=='Y'){
				i++;
			}else if(b=='I'){
				i++;
			}else if(b=='D'){
				d++;
			}else if(b=='S'){
				s++;
			}else{
				if(Character.isDigit(b) && shortmatch()){
					System.err.println("Warning! Found read in shortmatch form during countErrors():\n"+this); //Usually caused by verbose output.
					if(mate!=null){System.err.println("mate:\n"+mate.id+"\t"+new String(mate.bases));}
					System.err.println("Stack trace: ");
					new Exception().printStackTrace();
					match=toLongMatchString(match);
					setShortMatch(false);
					return countErrors();
				}else{
					throw new RuntimeException("\nUnknown symbol "+(char)b+":\n"+new String(match)+"\n"+this+"\nshortmatch="+this.shortmatch());
				}
			}
			
		}
		return new int[] {m, s, d, i, n};
	}
	
	public static byte[] toShortMatchString(byte[] match){
		if(match==null){return null;}
		assert(match.length>0);
		StringBuilder sb=new StringBuilder(10);
		
		byte prev=match[0];
		int count=1;
		for(int i=1; i<match.length; i++){
			byte m=match[i];
			assert(Character.isLetter(m) || m==0) : new String(match);
			if(m==0){System.err.println("Warning! Converting empty match string to short form.");}
			if(m==prev){count++;}
			else{
				sb.append((char)prev);
				if(count>2){sb.append(count);}
				else if(count==2){sb.append((char)prev);}
				prev=m;
				count=1;
			}
		}
		sb.append((char)prev);
		if(count>2){sb.append(count);}
		else if(count==2){sb.append((char)prev);}
		
		byte[] r=new byte[sb.length()];
		for(int i=0; i<sb.length(); i++){r[i]=(byte)sb.charAt(i);}
		return r;
	}
	
	public static byte[] toLongMatchString(byte[] shortmatch){
		if(shortmatch==null){return null;}
		assert(shortmatch.length>0);
		
		int count=0;
		int current=0;
		for(int i=0; i<shortmatch.length; i++){
			byte m=shortmatch[i];
			if(Character.isLetter(m)){
				count++;
				count+=(current>0 ? current-1 : 0);
				current=0;
			}else{
				assert(Character.isDigit(m));
				current=(current*10)+(m-48); //48 == '0'
			}
		}
		count+=(current>0 ? current-1 : 0);
		
		
		byte[] r=new byte[count];
		current=0;
		byte lastLetter='?';
		int j=0;
		for(int i=0; i<shortmatch.length; i++){
			byte m=shortmatch[i];
			if(Character.isLetter(m)){
				while(current>1){
					r[j]=lastLetter;
					current--;
					j++;
				}
				current=0;
				
				r[j]=m;
				j++;
				lastLetter=m;
			}else{
				assert(Character.isDigit(m));
				current=(current*10)+(m-48); //48 == '0'
			}
		}
		while(current>1){
			r[j]=lastLetter;
			current--;
			j++;
		}
		
		assert(r[r.length-1]>0);
		return r;
	}
	
	
//	/** Original bases of the read.  Do not modify! */
//	public byte[] basesOriginal;
//	
//	/** Quality code for read bases of the read.  Do not modify! */
//	public byte[] qualityOriginal;

	/** Bases of the read, after trimming (for colorspace). */
	public byte[] bases;
	
	/** Quality of the read, after trimming (for colorspace). */
	public byte[] quality;
	
	/** Alignment string.  E.G. mmmmDDDmmm would have 4 matching bases, then a 3-base deletion, then 3 matching bases. */
	public byte[] match;
	
	public int[] gaps;
	
	public String id;
	public long numericID;
	public int chrom;
	public int start;
	public int stop;

	public int mapLength; //Length used for mapping, before trimming
	
	public int copies=1;

	/** Errors detected (remaining) */
	public int errors=0;
	
	/** Errors corrected.  Total initial errors should be errors+errorsCorrected. */
	public int errorsCorrected=0;
	
	/** Alignment score from BBMap.  Assumed to max at approx 100*bases.length */
	public int mapScore=0;
	
	public ArrayList<SiteScore> sites;
	public SiteScore originalSite; //Origin site for synthetic reads
	public Object obj=null; //For testing only
	public Read mate;
	
	public int flags;
	
	/** -1 if invalid.  TODO: Currently not retained through most processes. */
	private int insert=-1;
	
	/** A random number for deterministic usage.
	 * May decrease speed in multithreaded applications.
	 */
	public double rand=-1;
	
	public byte strand(){return (byte)(flags&1);}
	public boolean mapped(){return (flags&MAPPEDMASK)==MAPPEDMASK;}
	public boolean paired(){return (flags&PAIREDMASK)==PAIREDMASK;}
	public boolean synthetic(){return (flags&SYNTHMASK)==SYNTHMASK;}
	public boolean colorspace(){return (flags&COLORMASK)==COLORMASK;}
	public boolean ambiguous(){return (flags&AMBIMASK)==AMBIMASK;}
	public boolean perfect(){return (flags&PERFECTMASK)==PERFECTMASK;}
//	public boolean semiperfect(){return perfect() ? true : list!=null && list.size()>0 ? list.get(0).semiperfect : false;} //TODO: This is a hack.  Add a semiperfect flag.
	public boolean rescued(){return (flags&RESCUEDMASK)==RESCUEDMASK;}
	public boolean discarded(){return (flags&DISCARDMASK)==DISCARDMASK;}
	public boolean invalid(){return (flags&INVALIDMASK)==INVALIDMASK;}
	public boolean swapped(){return (flags&SWAPMASK)==SWAPMASK;}
	public boolean shortmatch(){return (flags&SHORTMATCHMASK)==SHORTMATCHMASK;}
	public boolean insertvalid(){return (flags&INSERTMASK)==INSERTMASK;}
	public boolean hasadapter(){return (flags&ADAPTERMASK)==ADAPTERMASK;}
	public boolean secondary(){return (flags&SECONDARYMASK)==SECONDARYMASK;}
	/** For paired ends: 0 for read1, 1 for read2 */
	public int pairnum(){return (flags&PAIRNUMMASK)>>PAIRNUMSHIFT;}
	public boolean valid(){return !invalid();}

	public boolean getFlag(int mask){return (flags&mask)==mask;}
	public int flagToNumber(int mask){return (flags&mask)==mask ? 1 : 0;}
	
	public void setFlag(int mask, boolean b){
		flags=(flags&~mask);
		if(b){flags|=mask;}
	}
	
	public void setStrand(int b){
		assert(b==1 || b==0);
		flags=(flags&(~1))|b;
	}
	
	/** For paired ends: 0 for read1, 1 for read2 */
	public void setPairnum(int b){
		assert(b==1 || b==0);
		flags=(flags&(~PAIRNUMMASK))|(b<<PAIRNUMSHIFT);
//		assert(pairnum()==b);
	}
	
	public void setPaired(boolean b){
		flags=(flags&~PAIREDMASK);
		if(b){flags|=PAIREDMASK;}
	}
	
	public void setSynthetic(boolean b){
		flags=(flags&~SYNTHMASK);
		if(b){flags|=SYNTHMASK;}
	}
	
	public void setColorspace(boolean b){
		flags=(flags&~COLORMASK);
		if(b){flags|=COLORMASK;}
	}
	
	public void setAmbiguous(boolean b){
		flags=(flags&~AMBIMASK);
		if(b){flags|=AMBIMASK;}
	}
	
	public boolean setPerfectFlag(int maxScore){
		final SiteScore ss=topSite();
		if(ss==null){
			setPerfect(false);
		}else{
			assert(ss.slowScore<=maxScore) : maxScore+", "+ss.slowScore+", "+ss.toText();
			
			if(ss.slowScore==maxScore || ss.perfect){
				assert(testMatchPerfection(true)) : "\n"+ss+"\n"+maxScore+"\n"+this+"\n"+mate+"\n";
				setPerfect(true);
			}else{
				boolean flag=testMatchPerfection(false);
				setPerfect(flag);
				assert(flag || !ss.perfect) : "flag="+flag+", ss.perfect="+ss.perfect+"\nmatch="+new String(match)+"\n"+this.toText(false);
				assert(!flag || ss.slowScore>=maxScore) : "\n"+ss+"\n"+maxScore+"\n"+this+"\n"+mate+"\n";
			}
		}
		return perfect();
	}
	
	private boolean testMatchPerfection(boolean returnIfNoMatch){
		if(match==null){return returnIfNoMatch;}
		boolean flag=(match.length==bases.length);
		if(shortmatch()){
			flag=(match.length==0 || match[0]=='m');
			for(int i=0; i<match.length && flag; i++){flag=(match[i]=='m' || Character.isDigit(match[i]));}
		}else{
			for(int i=0; i<match.length && flag; i++){flag=(match[i]=='m');}
		}
		for(int i=0; i<bases.length && flag; i++){flag=(bases[i]!='N');}
		return flag;
	}
	
	public int insert(){
		return insertvalid() ? insert : -1;
	}
	
	public int insertSizeMapped(boolean ignoreStrand){
		return insertSizeMapped(this, mate, ignoreStrand);
	}
	
	public static int insertSizeMapped(Read r1, Read r2, boolean ignoreStrand){
		if(ignoreStrand || r2==null || !r1.mapped() || !r2.mapped() || r1.strand()==r2.strand()){return insertSizeMapped_Unstranded(r1, r2);}
		return insertSizeMapped_PlusLeft(r1, r2);
	}
	
	/** TODO: This is not correct when the insert is shorter than a read's bases with same-strand reads */
	public static int insertSizeMapped_PlusLeft(Read r1, Read r2){
		if(r1.strand()>r2.strand()){return insertSizeMapped_PlusLeft(r2, r1);}
		if(r1.strand()==r2.strand() || r1.start>r2.stop){return insertSizeMapped_Unstranded(r2, r1);} //So r1 is always on the left.
		
//		if(!mapped() || !mate.mapped()){return 0;}
		if(r1.chrom!=r2.chrom){return 0;}
		if(r1.start==r1.stop || r2.start==r2.stop){return 0;} //???
		
		int a=(r1.bases==null ? 0 : r1.bases.length);
		int b=(r2.bases==null ? 0 : r2.bases.length);
		int mid=r2.start-r1.stop-1;
		if(-mid>=a+b){return insertSizeMapped_Unstranded(r1, r2);} //Not properly oriented; plus read is to the right of minus read
		return mid+a+b;
	}
	
	public static int insertSizeMapped_Unstranded(Read r1, Read r2){
		if(r2==null){return r1.start==r1.stop ? 0 : r1.stop-r1.start+1;}
		
		if(r1.start>r2.start){return insertSizeMapped_Unstranded(r2, r1);} //So r1 is always on the left side.
		
//		if(!mapped() || !mate.mapped()){return 0;}
		if(r1.start==r1.stop || r2.start==r2.stop){return 0;} //???
		
		if(r1.chrom!=r2.chrom){return 0;}
		int a=(r1.bases==null ? 0 : r1.bases.length);
		int b=(r2.bases==null ? 0 : r2.bases.length);
		if(false && Tools.overlap(r1.start, r1.stop, r2.start, r2.stop)){
			//This does not handle very short inserts
			return Tools.max(r1.stop, r2.stop)-Tools.min(r1.start, r2.start)+1;
			
		}else{
			if(r1.start<r2.start){
				int mid=r2.start-r1.stop-1;
//				assert(false) : mid+", "+a+", "+b;
//				if(-mid>a && -mid>b){return Tools.min(a, b);} //Strange situation, no way to guess insert size
				if(-mid>=a+b){return 0;} //Strange situation, no way to guess insert size
				return mid+a+b;
			}else{
				assert(r1.start==r2.start);
				return Tools.min(a, b);
			}
		}
	}
	
	public int insertSizeOriginalSite(){
		if(mate==null){
//			System.err.println("A: "+(originalSite==null ? "null" : (originalSite.stop-originalSite.start+1)));
			return (originalSite==null ? 0 : originalSite.stop-originalSite.start+1);
		}
		
		final SiteScore ssa=originalSite, ssb=mate.originalSite;
		final int x;
		if(ssa==null || ssb==null){
//			System.err.println("B: 0");
			x=0;
		}else{
			x=insertSize(ssa, ssb, bases.length, mate.bases.length);
		}
		
		assert(pairnum()>=mate.pairnum() || x==mate.insertSizeOriginalSite());
		return x;
	}
	
	public static int insertSize(SiteScore ssa, SiteScore ssb, int lena, int lenb){
		return insertSize(ssa.chrom, ssb.chrom, ssa.start, ssb.start, ssa.stop, ssb.stop, lena, lenb);
	}
	
	public static int insertSize(int chroma, int chromb, int starta, int startb, int stopa, int stopb, int lena, int lenb){
		
		final int x;

		//		if(mate==null || ){return bases==null ? 0 : bases.length;}
		if(chroma!=chromb){x=0;}
		else{

			if(Tools.overlap(starta, stopa, startb, stopb)){
				x=Tools.max(stopa, stopb)-Tools.min(starta, startb)+1;
//				System.err.println("C: "+x);
			}else{
				if(starta<=startb){
					int mid=startb-stopa-1;
					//				assert(false) : mid+", "+a+", "+b;
					x=mid+lena+lenb;
//					System.err.println("D: "+x);
				}else{
					int mid=starta-stopb-1;
					//				assert(false) : mid+", "+a+", "+b;
					x=mid+lena+lenb;
//					System.err.println("E: "+x);
				}
			}
		}
		return x;
	}
	
	public Read joinRead(){
		if(insert<1 || mate==null || !insertvalid()){return this;}
		assert(insert>9 || bases.length<20) : "Perhaps old read format is being used?  This appears to be a quality value, not an insert.\n"+this+"\n\n"+mate+"\n";
		return joinRead(this, mate, insert);
	}
	
	public Read joinRead(int x){
		if(x<1 || mate==null){return this;}
		assert(x>9 || bases.length<20) : "Perhaps old read format is being used?  This appears to be a quality value, not an insert.\n"+this+"\n\n"+mate+"\n";
		return joinRead(this, mate, x);
	}
	
	public static Read joinRead(Read a, Read b, int insert){
		assert(a!=null && b!=null && insert>0);
		final int lengthSum=a.bases.length+b.bases.length;
		final int overlap=Tools.min(insert, lengthSum-insert);
		
		final byte[] bases=new byte[insert];
		final byte[] quals=new byte[insert];
		
		int mismatches=0;
		
		int start, stop;
		
		if(overlap<=0){//Simple join
			for(int i=0; i<a.bases.length; i++){
				bases[i]=a.bases[i];
				quals[i]=a.quality[i];
			}
			int lim=insert-b.bases.length;
			for(int i=a.bases.length; i<lim; i++){
				bases[i]='N';
				quals[i]=0;
			}
			for(int i=0; i<b.bases.length; i++){
				bases[i+lim]=b.bases[i];
				quals[i+lim]=b.quality[i];
			}
			
			start=Tools.min(a.start, b.start);
//			stop=start+insert-1;
			stop=Tools.max(a.stop, b.stop);
			
//		}else if(insert>=a.bases.length && insert>=b.bases.length){ //Overlapped join, proper orientation
//			final int lim1=a.bases.length-overlap;
//			final int lim2=a.bases.length;
//			for(int i=0; i<lim1; i++){
//				bases[i]=a.bases[i];
//				quals[i]=a.quality[i];
//			}
//			for(int i=lim1, j=0; i<lim2; i++, j++){
//				assert(false) : "TODO";
//				bases[i]='N';
//				quals[i]=0;
//			}
//			for(int i=lim2, j=overlap; i<bases.length; i++, j++){
//				bases[i]=b.bases[j];
//				quals[i]=b.quality[j];
//			}
		}else{ //reads go off ends of molecule.
			for(int i=0; i<a.bases.length && i<bases.length; i++){
				bases[i]=a.bases[i];
				quals[i]=a.quality[i];
			}
			for(int i=bases.length-1, j=b.bases.length-1; i>=0 && j>=0; i--, j--){
				byte ca=bases[i], cb=b.bases[j];
				byte qa=quals[i], qb=b.quality[j];
				if(ca==0 || ca=='N'){
					bases[i]=cb;
					quals[i]=qb;
				}else if(ca==cb){
					quals[i]=(byte)Tools.min((Tools.max(qa, qb)+Tools.min(qa, qb)/4), 50);
				}else{
					bases[i]=(ca>=cb ? ca : cb);
					quals[i]=(byte)(Tools.max(ca, cb)-Tools.min(ca, cb));
					if(ca!='N' && cb!='N'){mismatches++;}
				}
			}
			
			if(a.strand()==0){
				start=a.start;
//				stop=start+insert-1;
				stop=b.stop;
			}else{
				stop=a.stop;
//				start=stop-insert+1;
				start=b.start;
			}
			if(start>stop){
				start=Tools.min(a.start, b.start);
				stop=Tools.max(a.stop, b.stop);
			}
		}
//		assert(mismatches>=countMismatches(a, b, insert, 999));
//		System.err.println(mismatches);
		if(a.chrom==0 || start==stop || (!a.mapped() && !a.synthetic())){start=stop=a.chrom=0;}
		
		Read r=new Read(bases, a.chrom, start, stop, a.id, quals, a.numericID, a.flags);
		if(a.chrom==0 || start==stop || (!a.mapped() && !a.synthetic())){r.setMapped(true);}
		r.setInsert(insert);
		r.setPaired(false);
		r.copies=a.copies;
		r.mapScore=a.mapScore+b.mapScore;
		r.mapLength=a.mapLength+b.mapLength;
		if(overlap<=0){
			r.mapScore=a.mapScore+b.mapScore;
			r.mapLength=a.mapLength+b.mapLength;
			r.errors=a.errors+b.errors;
			r.errorsCorrected=a.errorsCorrected+b.errorsCorrected;
			//TODO r.gaps=?
		}else{//Hard to calculate
			r.mapScore=(int)((a.mapScore*(long)a.mapLength+b.mapScore*(long)b.mapLength)/insert);
			r.mapLength=insert;
			r.errors=a.errors;
			r.errorsCorrected=a.errorsCorrected;
		}
		
		
		assert(r.insertvalid()) : "\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";
		assert(r.insert()==r.bases.length) : r.insert()+"\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";
//		assert(false) : "\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";
		
		assert(Shared.anomaly || (a.insertSizeMapped(false)>0 == r.insertSizeMapped(false)>0)) : 
			"\n"+r.insert()+"\n"+r.insertSizeMapped(false)+"\n"+a.insert()+"\n"+a.insertSizeMapped(false)+
			"\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n"; 
		
		return r;
	}

	/**
	 * @param minlen
	 * @param maxlen
	 * @return
	 */
	public ArrayList<Read> split(int minlen, int maxlen) {
		int len=bases==null ? 0 : bases.length;
		if(len<minlen){return null;}
		int parts=(len+maxlen-1)/maxlen;
		ArrayList<Read> subreads=new ArrayList<Read>(parts);
		if(len<=maxlen){
			subreads.add(this);
		}else{
			float ideal=Tools.max(minlen, len/(float)parts);
			int actual=(int)ideal;
			assert(false) : "TODO"; //Some assertion goes here, I forget what
			for(int i=0; i<parts; i++){
				int a=i*actual;
				int b=(i+1)*actual;
				if(b>bases.length){b=bases.length;}
//				if(b-a<)
				byte[] subbases=Arrays.copyOfRange(bases, a, b);
				byte[] subquals=(quality==null ? null : Arrays.copyOfRange(quality, a, b+1));
				Read r=new Read(subbases, -1, -1, -1, id+"_"+i, subquals, numericID, flags);
				subreads.add(r);
			}
		}
		return subreads;
	}
	
	/** Generate and return an array of canonical kmers for this read */
	public long[] toKmers(final int k, final int gap, long[] kmers, boolean makeCanonical) {
		if(gap>0){throw new RuntimeException("Gapped reads: TODO");}
		if(k>31){return toLongKmers(k, gap, kmers, makeCanonical);}
		if(bases==null || bases.length<k+gap){return null;}
		
		final int kbits=2*k;
		final long mask=~((-1L)<<(kbits));
		
		int len=0;
		long kmer=0;
		final int arraylen=bases.length-k+1;
		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
		Arrays.fill(kmers, -1);
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;

				if(len>=k){
					kmers[i-k+1]=kmer;
				}
			}
		}
		
//		System.out.println(new String(bases));
//		System.out.println(Arrays.toString(kmers));
		
		if(makeCanonical){
			this.reverseComplement();
			len=0;
			kmer=0;
			for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
				byte b=bases[i];
				int x=AminoAcid.baseToNumber[b];
				if(x<0){
					len=0;
					kmer=0;
				}else{
					kmer=((kmer<<2)|x)&mask;
					len++;

					if(len>=k){
						assert(kmer==AminoAcid.reverseComplementBinaryFast(kmers[j], k));
						kmers[j]=Tools.max(kmers[j], kmer);
					}
				}
			}
			this.reverseComplement();
			
//			System.out.println(Arrays.toString(kmers));
		}
		
		
		return kmers;
	}
	
	/** Generate and return an array of canonical kmers for this read */
	public long[] toLongKmers(final int k, final int gap, long[] kmers, boolean makeCanonical) {
		if(gap>0){throw new RuntimeException("Gapped reads: TODO");}
		assert(k>31) : k;
		if(bases==null || bases.length<k+gap){return null;}
		
		final int kbits=2*k;
		final long mask=Long.MAX_VALUE;
		
		int len=0;
		long kmer=0;
		final int arraylen=bases.length-k+1;
		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
		Arrays.fill(kmers, -1);
		
		
		final int tailshift=k%32;
		final int tailshiftbits=tailshift*2;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=Long.rotateLeft(kmer, 2);
				kmer=kmer^x;
				len++;

				if(len>=k){
					long x2=AminoAcid.baseToNumber[bases[i-k]];
					kmer=kmer^(x2<<tailshiftbits);
					kmers[i-k+1]=kmer;
				}
			}
		}
		if(makeCanonical){
			this.reverseComplement();
			len=0;
			kmer=0;
			for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
				byte b=bases[i];
				int x=AminoAcid.baseToNumber[b];
				if(x<0){
					len=0;
					kmer=0;
				}else{
					kmer=Long.rotateLeft(kmer, 2);
					kmer=kmer^x;
					len++;

					if(len>=k){
						long x2=AminoAcid.baseToNumber[bases[i-k]];
						kmer=kmer^(x2<<tailshiftbits);
						kmers[j]=mask&(Tools.max(kmers[j], kmer));
					}
				}
			}
			this.reverseComplement();
		}else{
			assert(false) : "Long kmers should be made canonical here because they cannot be canonicized later.";
		}
		
		return kmers;
	}
	
	
	public static final boolean CHECKSITES(Read r, byte[] basesM){
		return CHECKSITES(r.sites, r.bases, basesM, r.numericID);
	}
	
	public static final boolean CHECKSITES(ArrayList<SiteScore> list, byte[] basesP, byte[] basesM, long id){
		if(list==null || list.isEmpty()){return true;}
		for(int i=0; i<list.size(); i++){
			SiteScore ss=list.get(i);
			if(ss.strand==Gene.MINUS && basesM==null && basesP!=null){basesM=AminoAcid.reverseComplementBases(basesP, false);}
			byte[] bases=(ss.strand==Gene.PLUS ? basesP : basesM);
//			System.err.println("Checking site "+i);
			boolean b=CHECKSITE(ss, bases, id);
			assert(b);
//			System.err.println("Checked site "+i+" = "+ss+"\nss.p="+ss.perfect+", ss.sp="+ss.semiperfect);
			if(!b){
//				System.err.println("Error at SiteScore "+i+": ss.p="+ss.perfect+", ss.sp="+ss.semiperfect);
				return false;
			}
		}
		return true;
	}
	
	/** Make sure 'bases' is for correct strand! */
	public static final boolean CHECKSITE(SiteScore ss, byte[] basesP, byte[] basesM, long id){
		return CHECKSITE(ss, ss.plus() ? basesP : basesM, id);
	}
	
	/** Make sure 'bases' is for correct strand! */
	public static final boolean CHECKSITE(SiteScore ss, byte[] bases, long id){
		if(ss==null){return true;}
//		System.err.println("Checking site "+ss+"\nss.p="+ss.perfect+", ss.sp="+ss.semiperfect+", bases="+new String(bases));
		if(ss.perfect){assert(ss.semiperfect) : ss+"\n"+new String(bases);}
		if(bases!=null){

			final boolean p=ss.perfect;
			final boolean sp=ss.semiperfect;
			final boolean p1=ss.isPerfect(bases);
			final boolean sp1=(p1 ? true : ss.isSemiPerfect(bases));

			assert(p==p1) : p+"->"+p1+", "+sp+"->"+sp1+", "+ss.isSemiPerfect(bases)+
				"\nnumericID="+id+"\n"+new String(bases)+"\n\n"+Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n\n";
			assert(sp==sp1) : p+"->"+p1+", "+sp+"->"+sp1+", "+ss.isSemiPerfect(bases)+
				"\nnumericID="+id+"\n"+new String(bases)+"\n\n"+Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n\n";
			
//			ss.setPerfect(bases, false);
			
			assert(p==ss.perfect) : 
				p+"->"+ss.perfect+", "+sp+"->"+ss.semiperfect+", "+ss.isSemiPerfect(bases)+"\nnumericID="+id+"\n"+new String(bases)+"\n\n"+
				Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n\n";
			assert(sp==ss.semiperfect) : 
				p+"->"+ss.perfect+", "+sp+"->"+ss.semiperfect+", "+ss.isSemiPerfect(bases)+"\nnumericID="+id+"\n"+new String(bases)+"\n\n"+
				Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n\n";
			if(ss.perfect){assert(ss.semiperfect);}
		}
		return true;
	}
	
	
	public void setPerfect(boolean b){
		flags=(flags&~PERFECTMASK);
		if(b){flags|=PERFECTMASK;}
	}
	
	public void setRescued(boolean b){
		flags=(flags&~RESCUEDMASK);
		if(b){flags|=RESCUEDMASK;}
	}
	
	public void setMapped(boolean b){
//		assert(false) : mapped()+"->"+b;
		flags=(flags&~MAPPEDMASK);
		if(b){flags|=MAPPEDMASK;}
	}
	
	public void setDiscarded(boolean b){
		flags=(flags&~DISCARDMASK);
		if(b){flags|=DISCARDMASK;}
	}
	
	public void setInvalid(boolean b){
		flags=(flags&~INVALIDMASK);
		if(b){flags|=INVALIDMASK;}
	}
	
	public void setSwapped(boolean b){
		flags=(flags&~SWAPMASK);
		if(b){flags|=SWAPMASK;}
	}
	
	public void setShortMatch(boolean b){
		flags=(flags&~SHORTMATCHMASK);
		if(b){flags|=SHORTMATCHMASK;}
	}
	
	public void setInsertValid(boolean b){
		flags=(flags&~INSERTMASK);
		if(b){flags|=INSERTMASK;}
	}
	
	public void setHasAdapter(boolean b){
		flags=(flags&~ADAPTERMASK);
		if(b){flags|=ADAPTERMASK;}
	}
	
	public void setSecondary(boolean b){
		flags=(flags&~SECONDARYMASK);
		if(b){flags|=SECONDARYMASK;}
	}
	
	public void setInsert(int x){
		if(x<1){x=-1;}
		assert(x==-1 || x>9 || bases.length<20);
		insert=x;
		setInsertValid(x>0);
		if(mate!=null){
			mate.insert=x;
			mate.setInsertValid(x>0);
		}
	}

	private static int[] makeMaskArray(int max) {
		int[] r=new int[max+1];
		for(int i=0; i<r.length; i++){r[i]=(1<<i);}
		return r;
	}
	


	public static byte[] getFakeQuality(int len){
		if(len>=QUALCACHE.length){
			byte[] r=new byte[len];
			Arrays.fill(r, (byte)30);
			return r;
		}
		if(QUALCACHE[len]==null){
			synchronized(QUALCACHE){
				if(QUALCACHE[len]==null){
					QUALCACHE[len]=new byte[len];
					Arrays.fill(QUALCACHE[len], (byte)30);
				}
			}
		}
		return QUALCACHE[len];
	}
	
	
	public byte[] getScaffoldName(boolean requireSingleScaffold){
		byte[] name=null;
		if(mapped()){
			if(!requireSingleScaffold || Data.isSingleScaffold(chrom, start, stop)){
				int idx=Data.scaffoldIndex(chrom, (start+stop)/2);
				name=Data.scaffoldNames[chrom][idx];
//				int scaflen=Data.scaffoldLengths[chrom][idx];
//				a1=Data.scaffoldRelativeLoc(chrom, start, idx);
//				b1=a1-start1+stop1;
			}
		}
		return name;
	}
	
	public Read clone(){
		try {
			return (Read) super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException();
	}
	
	private static final byte[][] QUALCACHE=new byte[1000][];
	

	public static final int STRANDMASK=1;
	public static final int MAPPEDMASK=(1<<1);
	public static final int PAIREDMASK=(1<<2);
	public static final int PERFECTMASK=(1<<3);
	public static final int AMBIMASK=(1<<4);
	public static final int RESCUEDMASK=(1<<5);
	public static final int COLORMASK=(1<<6);
	public static final int SYNTHMASK=(1<<7);
	public static final int DISCARDMASK=(1<<8);
	public static final int INVALIDMASK=(1<<9);
	public static final int SWAPMASK=(1<<10);
	public static final int SHORTMATCHMASK=(1<<11);
	
	public static final int PAIRNUMSHIFT=12;
	public static final int PAIRNUMMASK=(1<<PAIRNUMSHIFT);

	public static final int INSERTMASK=(1<<13);
	public static final int ADAPTERMASK=(1<<14);
	public static final int SECONDARYMASK=(1<<15);
	
	private static final int[] maskArray=makeMaskArray(15); //Be sure this is big enough for all flags!

//	public static byte ASCII_OFFSET=33;
	private static final byte ASCII_OFFSET=33;
	public static byte MIN_CALLED_QUALITY=2;
	public static boolean TO_UPPER_CASE=false;
	
	public static boolean COMPRESS_MATCH_BEFORE_WRITING=true;
	public static boolean DECOMPRESS_MATCH_ON_LOAD=true; //Set to false for some applications, like sorting, perhaps
	
	public static boolean ADD_BEST_SITE_TO_LIST_FROM_TEXT=true;
	public static boolean NULLIFY_BROKEN_QUALITY=false;
}
