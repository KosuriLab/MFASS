package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import align2.Shared;
import align2.Tools;


import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;
import dna.ScafLoc;
import fileIO.TextStreamWriter;


public class SamLine {
	
//	426_647_582	161	chr1	10159	0	26M9H	chr3	170711991	0	TCCCTAACCCTAACCCTAACCTAACC	IIFIIIIIIIIIIIIIIIIIICH2<>	RG:Z:20110708003021394	NH:i:3	CM:i:2	SM:i:1	CQ:Z:A9?(BB?:<A?>=>B67=:7A);.%8'%))/%*%'	CS:Z:G12002301002301002301023010200000003	XS:A:+
	
//	1 QNAME String [!-?A-~]f1,255g Query template NAME
//	2 FLAG Int [0,216-1] bitwise FLAG
//	3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
//	4 POS Int [0,229-1] 1-based leftmost mapping POSition
//	5 MAPQ Int [0,28-1] MAPping Quality
//	6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
//	7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next fragment
//	8 PNEXT Int [0,229-1] Position of the mate/next fragment
//	9 TLEN Int [-229+1,229-1] observed Template LENgth
//	10 SEQ String \*|[A-Za-z=.]+ fragment SEQuence
//	11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
	
	
//	2. FLAG: bitwise FLAG. Each bit is explained in the following table:
//		Bit Description
//		0x1 template having multiple fragments in sequencing
//		0x2 each fragment properly aligned according to the aligner
//		0x4 fragment unmapped
//		0x8 next fragment in the template unmapped
//		0x10 SEQ being reverse complemented
//		0x20 SEQ of the next fragment in the template being reversed
//		0x40 the first fragment in the template
//		0x80 the last fragment in the template
//		0x100 secondary alignment
//		0x200 not passing quality controls
//		0x400 PCR or optical duplicate
	
	
//	FCB062MABXX:1:1101:1177:2115#GGCTACAA	147	chr11	47765857	29	90M	=	47765579	-368	CCTCTGTGGCCCGGGTTGGAGTGCAGTGTCATGATCATGGCTCGCTGTAGCTACACCCTTCTGAGCTCAAGCAATCCTCCCACCTCTCCC	############################################################A@@><D<AAAB<=A2BD/BC<7:<4<%679	XT:A:M	NM:i:5	SM:i:29	AM:i:29	XM:i:5	XO:i:0	XG:i:0	MD:Z:7T4A15G26A30A3
//	FCB062MABXX:1:1101:1193:2122#GGCTACAA	77	*	    0	         0	*	*	0	           0	TATATATGTGCTATGTACAGCATTGGAATTCACACCCTACACTTTCAAAAGNGAGCCCTAAATAAATGTTAGATCGGAAGAGCACACGTC	FCFCFDDDADDEDEBDAEDFEDEFFGGFGGHEEFHHHHHHEDDDEDFFEFB#CBBA@B8BGGFGEEEC>DGGGDFBGGGGHHHHH9<@##

	
//	public static StringBuilder header(int minChrom, int maxChrom, boolean firstLine, boolean lastLine){
//		StringBuilder sb=new StringBuilder(1000);
//		if(firstLine){
//			sb.append("@HD\tVN:1.0\tSO:unsorted\n");
//		}
////		assert(false) : minChrom+", "+maxChrom+", "+Data.numChroms;
//		for(int i=minChrom; i<=maxChrom && i<=Data.numChroms; i++){
//			//			sb.append("\n"+Data.getChromosome(i).getString(0, 2100)+"\n");
//
//			//			sb.append("@SQ\tSN:chr"+i);
//
//			for(int j=0; j<Data.chromScaffolds[i]; j++){
//				sb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
//				if(Data.scaffoldNames[i][j]==null){
//					assert(false) : "scaffoldName["+i+"]["+j+"] = null";
//					sb.append("null");
//				}else{
//					byte[] scn=Data.scaffoldNames[i][j];
//					for(byte b : scn){sb.append((char)b);}
//				}
//
//				//			if(Data.chromosomePlusMatrix[i]!=null){sb.append("\tLN:"+Data.getChromosome(i).maxIndex);}
//				//			else if(Data.chromosomePlusMatrixColorspace[i]!=null){sb.append("\tLN:"+Data.getChromosome(i, true).maxIndex);}
//
//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
//
//				//			sb.append("\tLN:"+(Data.getChromosome(i).maxIndex+100));
//				//			if(UNLOAD_CHROMS_WHEN_MAKING_HEADER){Data.unload(i, true);}
//
//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"build "+Data.GENOME_BUILD).replace('\t', ' '));
//
//				sb.append("\n");
//			}
//		}
//		if(lastLine){
//			sb.append("@RG\tID:unknownRG\tSM:unknownSM\tPL:ILLUMINA\n");
//			sb.append("@PG\tID:BBMap\tPN:BBMap\tVN:"+Shared.BBMAP_VERSION_STRING);
//		}
////		sb.append("@CO\tno comment\n");
//		
////		assert(false) : sb;
//		
//		return sb;
//	}
	
	public static ByteBuilder header0B(ByteBuilder bb){
//		if(MAKE_TOPHAT_TAGS){
//			return new ByteBuilder("@HD\tVN:"+(VERSION<1.4f ? "1.0" : "1.4")+"\tSO:unsorted");
//		}
		bb.append("@HD\tVN:");
		bb.append((VERSION<1.4f ? "1.3" : "1.4"));
		bb.append("\tSO:unsorted");
		return bb;
	}
	
	public static StringBuilder header0(){
//		if(MAKE_TOPHAT_TAGS){
//			return new StringBuilder("@HD\tVN:"+(VERSION<1.4f ? "1.0" : "1.4")+"\tSO:unsorted");
//		}
		StringBuilder sb=new StringBuilder("@HD\tVN:"+(VERSION<1.4f ? "1.3" : "1.4")+"\tSO:unsorted");
		return sb;
	}
	
	private static ArrayList<String> scaffolds(int minChrom, int maxChrom, boolean sort){
		final ArrayList<String> list=new ArrayList<String>(4000);
		final StringBuilder sb=new StringBuilder(1000);
		for(int i=minChrom; i<=maxChrom && i<=Data.numChroms; i++){
			final byte[][] inames=Data.scaffoldNames[i];
			for(int j=0; j<Data.chromScaffolds[i]; j++){
				final byte[] scn=inames[j];
				sb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scn==null){
					assert(false) : "scaffoldName["+i+"]["+j+"] = null";
					sb.append("null");
				}else{
					appendScafName(sb, scn);
				}
				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j])));
//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"b"+Data.GENOME_BUILD).replace('\t', ' '));

				sb.append('\n');
				list.add(sb.toString());
				sb.setLength(0);
			}
		}
		if(sort){Collections.sort(list);}
		return list;
	}
	
	public static StringBuilder header1(int minChrom, int maxChrom){
		StringBuilder sb=new StringBuilder(20000);
		if(SORT_SCAFFOLDS){
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true);
			for(int i=0; i<scaffolds.size(); i++){
				sb.append(scaffolds.get(i));
				scaffolds.set(i, null);
			}
			return sb;
		}
		
		for(int i=minChrom; i<=maxChrom && i<=Data.numChroms; i++){
			final byte[][] inames=Data.scaffoldNames[i];
			for(int j=0; j<Data.chromScaffolds[i]; j++){
				byte[] scn=inames[j];
				sb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scn==null){
					assert(false) : "scaffoldName["+i+"]["+j+"] = null";
					sb.append("null");
				}else{
					appendScafName(sb, scn);
				}

				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j])));
//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"build "+Data.GENOME_BUILD).replace('\t', ' '));

				sb.append('\n');
			}
		}
		
		return sb;
	}
	
	public static void printHeader1(int minChrom, int maxChrom, PrintWriter pw){
		if(SORT_SCAFFOLDS){
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true);
			for(int i=0; i<scaffolds.size(); i++){
				pw.print(scaffolds.set(i, null));
			}
			return;
		}
		
		for(int i=minChrom; i<=maxChrom && i<=Data.numChroms; i++){
			final byte[][] inames=Data.scaffoldNames[i];
			StringBuilder sb=new StringBuilder(256);
			for(int j=0; j<Data.chromScaffolds[i]; j++){
				final byte[] scn=inames[j];
//				StringBuilder sb=new StringBuilder(7+(scn==null ? 4 : scn.length)+4+10+4+/*(Data.name==null ? 0 : Data.name.length()+1)+11*/+4);//last one could be 1
				sb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scn==null){
					assert(false) : "scaffoldName["+i+"]["+j+"] = null";
					sb.append("null");
				}else{
					appendScafName(sb, scn);
				}
				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j])));
//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"b"+Data.GENOME_BUILD).replace('\t', ' '));

				sb.append('\n');
				
				pw.print(sb);
				sb.setLength(0);
			}
		}
	}
	
	public static void printHeader1B(int minChrom, int maxChrom, ByteBuilder bb, OutputStream os){
		if(verbose){System.err.println("printHeader1B("+minChrom+", "+maxChrom+")");}
		
		if(SORT_SCAFFOLDS){
			if(verbose){System.err.println("Sorting scaffolds");}
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true);
			for(int i=0; i<scaffolds.size(); i++){
				String s=scaffolds.set(i, null);
				bb.append(s);
				if(bb.length>=32768){
					try {
						os.write(bb.array, 0, bb.length);
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
					bb.setLength(0);
				}
			}
			return;
		}
		
		if(verbose){System.err.println("Iterating over chroms");}
		for(int chrom=minChrom; chrom<=maxChrom && chrom<=Data.numChroms; chrom++){
//			if(verbose){System.err.println("chrom "+chrom);}
			final byte[][] inames=Data.scaffoldNames[chrom];
//			if(verbose){System.err.println("inames"+(inames==null ? " = null" : ".length = "+inames.length));}
			final int numScafs=Data.chromScaffolds[chrom];
//			if(verbose){System.err.println("scaffolds: "+numScafs);}
			assert(inames.length==numScafs) : "Mismatch between number of scaffolds and names for chrom "+chrom+": "+inames.length+" != "+numScafs;
			for(int scaf=0; scaf<numScafs; scaf++){
//				if(verbose){System.err.println("chromScaffolds["+scaf+"] = "+(inames==null ? "=null" : ".length="+inames.length));}
				final byte[] scafName=inames[scaf];
//				if(verbose){System.err.println("scafName = "+(scafName==null ? "null" : new String(scafName)));}
				bb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scafName==null){
					assert(false) : "scaffoldName["+chrom+"]["+scaf+"] = null";
					bb.append(scafName);
				}else{
					appendScafName(bb, scafName);
				}
				bb.append("\tLN:");
				bb.append(Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[chrom][scaf])));
//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"b"+Data.GENOME_BUILD).replace('\t', ' '));

				bb.append('\n');
				
				if(bb.length>=32768){
					try {
						os.write(bb.array, 0, bb.length);
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
					bb.setLength(0);
				}
			}
		}
	}
	
	public static void printHeader1(int minChrom, int maxChrom, TextStreamWriter tsw){
		if(SORT_SCAFFOLDS){
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true);
			for(int i=0; i<scaffolds.size(); i++){
				tsw.print(scaffolds.set(i, null));
			}
			return;
		}
		
		for(int i=minChrom; i<=maxChrom && i<=Data.numChroms; i++){
			final byte[][] inames=Data.scaffoldNames[i];
			final StringBuilder sb=new StringBuilder(256);
			for(int j=0; j<Data.chromScaffolds[i]; j++){
				final byte[] scn=inames[j];
//				StringBuilder sb=new StringBuilder(7+(scn==null ? 4 : scn.length)+4+10+4+/*(Data.name==null ? 0 : Data.name.length()+1)+11*/+4);//last one could be 1
				sb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scn==null){
					assert(false) : "scaffoldName["+i+"]["+j+"] = null";
					sb.append("null");
				}else{
					appendScafName(sb, scn);
				}
				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j])));
//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"b"+Data.GENOME_BUILD).replace('\t', ' '));

				sb.append('\n');
				
				tsw.print(sb);
				sb.setLength(0);
			}
		}
	}
	
	private static void appendScafName(StringBuilder sb, byte[] scn){
		if(Data.scaffoldPrefixes){
			int k=0;
			while(k<scn.length && scn[k]!='$'){k++;}
			k++;
			while(k<scn.length){
				sb.append((char)scn[k]);
				k++;
			}
		}else{
			final char[] buffer=Shared.getTLCB(scn.length);
			for(int i=0; i<scn.length; i++){buffer[i]=(char)scn[i];}
			sb.append(buffer, 0, scn.length);
		}
	}
	
	private static void appendScafName(ByteBuilder sb, byte[] scn){
		if(Data.scaffoldPrefixes){
			int k=0;
			while(k<scn.length && scn[k]!='$'){k++;}
			k++;
			while(k<scn.length){
				sb.append(scn[k]);
				k++;
			}
		}else{
			sb.append(scn);
		}
	}
	
	public static StringBuilder header2(){
		StringBuilder sb=new StringBuilder(1000);
//		sb.append("@RG\tID:unknownRG\tSM:unknownSM\tPL:ILLUMINA\n"); //Can cause problems.  If RG is in the header, reads may need extra fields.
		
//		if(MAKE_TOPHAT_TAGS){
////			sb.append("@PG\tID:TopHat\tVN:2.0.6\tCL:/usr/common/jgi/aligners/tophat/2.0.6/bin/tophat -p 16 -r 0 --max-multihits 1 Creinhardtii_236 reads_1.fa reads_2.fa");
//			sb.append("@PG\tID:TopHat\tVN:2.0.6");
//		}else{
//			sb.append("@PG\tID:BBMap\tPN:BBMap\tVN:"+Shared.BBMAP_VERSION_STRING);
//		}
		sb.append("@PG\tID:BBMap\tPN:BBMap\tVN:"+Shared.BBMAP_VERSION_STRING);
		
		if(Shared.BBMAP_CLASS!=null){
			sb.append("\tCL:java");
			{
				List<String> list=null;
				list=Shared.JVM_ARGS();
				if(list!=null){
					for(String s : list){
						sb.append(' ');
						sb.append(s);
					}
				}
			}
			sb.append(' ');
			sb.append(Shared.BBMAP_CLASS);
			if(Shared.COMMAND_LINE!=null){
				for(String s : Shared.COMMAND_LINE){
					sb.append(' ');
					sb.append(s);
				}
			}
		}
			
		return sb;
	}
	
	public static ByteBuilder header2B(ByteBuilder sb){
		sb.append("@PG\tID:BBMap\tPN:BBMap\tVN:"+Shared.BBMAP_VERSION_STRING);
		
		if(Shared.BBMAP_CLASS!=null){
			sb.append("\tCL:java");
			{
				List<String> list=null;
				list=Shared.JVM_ARGS();
				if(list!=null){
					for(String s : list){
						sb.append(' ');
						sb.append(s);
					}
				}
			}
			sb.append(' ');
			sb.append(Shared.BBMAP_CLASS);
			if(Shared.COMMAND_LINE!=null){
				for(String s : Shared.COMMAND_LINE){
					sb.append(' ');
					sb.append(s);
				}
			}
		}
			
		return sb;
	}

	public static final boolean KILL_BAD_PAIRS=false;
	public static final boolean REQUIRE_CORRECT_STRANDS_PAIRS=true;
	public static final boolean SAME_STRAND_PAIRS=false;
	
	public SamLine(String s){
		this(s.split("\t"));
	}
	
	/** Prevents references to original string, in case of e.g. very long MD tags. */
	public SamLine toSamLine(String s){
		String[] split=s.split("\t");
		split[0]=new String(split[0]);
		split[5]=new String(split[5]);
		split[9]=new String(split[9]);
		split[10]=new String(split[10]);
		for(int i=11; i<split.length; i++){
			split[i]=new String(split[i]);
		}
		return new SamLine(split);
	}
	
//	@Deprecated
//	private static String toCigarOld(byte[] match){
//		if(match==null){return stringstar;}
//		StringBuilder sb=new StringBuilder(8);
//		int count=0;
//		char mode='M';
//		char lastMode='M';
//		for(int mpos=0; mpos<match.length; mpos++){
//
//			byte m=match[mpos];
//
//			if(m=='m' || m=='s' || m=='S' || m=='N'){
//				mode='M';
//			}else if(m=='I' || m=='X' || m=='Y'){
//				mode='I';
//			}else if(m=='D'){
//				mode='D';
//			}
//
//			if(mode!=lastMode){
//				if(count>0){//Prevents an initial length-0 match
//					sb.append(count);
//					sb.append(lastMode);
//				}
//				count=0;
//				lastMode=mode;
//			}
//
//			count++;
//		}
//		sb.append(count);
//		sb.append(mode);
//		return sb.toString();
//	}
	
	public static String toCigar13(byte[] match, int readStart, int readStop, int reflen, byte[] bases){
		if(match==null || readStart==readStop){return null;}
		StringBuilder sb=new StringBuilder(8);
		int count=0;
		char mode='=';
		char lastMode='=';
		
		int refloc=readStart;

		int cigarlen=0; //for debugging
		int opcount=0; //for debugging
		
		for(int mpos=0; mpos<match.length; mpos++){

			byte m=match[mpos];
			
			boolean sfdflag=false;
			if(SOFT_CLIP && (refloc<0 || refloc>=reflen)){
				mode='S'; //soft-clip out-of-bounds
				if(m!='I'){refloc++;}
				if(m=='D'){sfdflag=true;} //Don't add soft-clip count for deletions!
			}else if(m=='m' || m=='s' || m=='S' || m=='N' || m=='B'){//Little 's' is for a match classified as a sub to improve the affine score.  
				mode='M';
				refloc++;
			}else if(m=='I' || m=='X' || m=='Y'){
				mode='I';
			}else if(m=='D'){
				mode='D';
				refloc++;
			}else if(m=='C'){
				mode='S';
				refloc++;
			}else{
				throw new RuntimeException("Invalid match string character '"+(char)m+"' = "+m+" (ascii).  " +
						"Match string should be in long format here.");
			}

			if(mode!=lastMode){
				if(count>0){//Prevents an initial length-0 match
					sb.append(count);
//					sb.append(lastMode);
					if(lastMode=='D' && count>INTRON_LIMIT){sb.append('N');}
					else{sb.append(lastMode);}
					if(lastMode!='D'){cigarlen+=count;}
					opcount+=count;
				}
				count=0;
				lastMode=mode;
			}

			count++;
			if(sfdflag){count--;}
		}
		sb.append(count);
		if(mode=='D' && count>INTRON_LIMIT){sb.append('N');}
		else{sb.append(mode);}
		if(mode!='D'){cigarlen+=count;}
		opcount+=count;
		
		assert(bases==null || cigarlen==bases.length) : "\n(cigarlen = "+cigarlen+") != (bases.length = "+(bases==null ? -1 : bases.length)+")\n" +
				"cigar = "+sb+"\nmatch = "+new String(match)+"\nbases = "+new String(bases)+"\n";
		
		return sb.toString();
	}
	
	public static String toCigar14(byte[] match, int readStart, int readStop, int reflen, byte[] bases){
		if(match==null || readStart==readStop){return null;}
		StringBuilder sb=new StringBuilder(8);
		int count=0;
		char mode='=';
		char lastMode='=';
		
		int refloc=readStart;

		int cigarlen=0; //for debugging
		int opcount=0; //for debugging
		
		for(int mpos=0; mpos<match.length; mpos++){

			byte m=match[mpos];
			
			boolean sfdflag=false;
			if(SOFT_CLIP && (refloc<0 || refloc>=reflen)){
				mode='S'; //soft-clip out-of-bounds
				if(m!='I'){refloc++;}
				if(m=='D'){sfdflag=true;} //Don't add soft-clip count for deletions!
			}else if(m=='m' || m=='s'){//Little 's' is for a match classified as a sub to improve the affine score.  
				mode='=';
				refloc++;
			}else if(m=='S'){
				mode='X';
				refloc++;
			}else if(m=='I' || m=='X' || m=='Y'){
				mode='I';
			}else if(m=='D'){
				mode='D';
				refloc++;
			}else if(m=='C'){
				mode='S';
				refloc++;
			}else if(m=='N' || m=='B'){
				mode='M';
				refloc++;
			}else{
				throw new RuntimeException("Invalid match string character '"+(char)m+"' = "+m+" (ascii).  " +
						"Match string should be in long format here.");
			}

			if(mode!=lastMode){
				if(count>0){//Prevents an initial length-0 match
					sb.append(count);
					if(lastMode=='D' && count>INTRON_LIMIT){sb.append('N');}
					else{sb.append(lastMode);}
					if(lastMode!='D'){cigarlen+=count;}
					opcount+=count;
				}
				count=0;
				lastMode=mode;
			}

			count++;
			if(sfdflag){count--;}
		}
		sb.append(count);
		if(mode=='D' && count>INTRON_LIMIT){
			sb.append('N');
		}else{
			sb.append(mode);
		}
		if(mode!='D'){cigarlen+=count;}
		opcount+=count;
		
		assert(bases==null || cigarlen==bases.length) : "\n(cigarlen = "+cigarlen+") != (bases.length = "+(bases==null ? -1 : bases.length)+")\n" +
				"cigar = "+sb+"\nmatch = "+new String(match)+"\nbases = "+new String(bases)+"\n";
		
		return sb.toString();
	}
	

	public static String makeStopTag(int pos, int seqLength, String cigar, boolean perfect){
//		return String.format("YS:i:%d", pos+(cigar==null ? seqLength : -countLeadingClip(cigar)+calcCigarLength(cigar))-1);
		return "YS:i:"+(pos+((cigar==null || perfect) ? seqLength : -countLeadingClip(cigar)+calcCigarLength(cigar))-1);
	}
	
	public static String makeIdentityTag(byte[] match, boolean perfect){
		if(perfect){return "YI:f:100";}
		float f=Read.identity(match);
		return String.format("YI:f:%.2f", (100*f));
	}
	
	
	public static String makeMdTag(int chrom, int refstart, byte[] match, byte[] call, boolean colorspace){
		if(match==null || chrom<0 || colorspace){return null;}
		StringBuilder md=new StringBuilder(8);
		md.append("MD:Z:");
		
		if(colorspace){throw new RuntimeException("Colorspace is incompatible with MD tags.");}
		
		ChromosomeArray cha=Data.getChromosome(chrom);
		
		byte prevM='?';
		int count=0;
		int dels=0;
//		inr cpos=0;
		for(int mpos=0, rpos=refstart; mpos<match.length; mpos++){
//			byte c=call[cpos];
			byte m=match[mpos];
			
			if(prevM=='D' && m!='D'){
				if(dels<=INTRON_LIMIT){//Otherwise, ignore it
					md.append(count);
					count=0;
					md.append('^');
					for(int i=rpos-dels; i<rpos; i++){
						md.append((char)cha.get(i));
					}
					dels=0;
				}
			}
			
			if(m=='m' || m=='s'){
				count++;
				rpos++;
//				cpos++;
			}else if(m=='S'){
				md.append(count);
				md.append((char)cha.get(rpos));

				count=0;
				rpos++;
//				cpos++;
			}else if(m=='N'){
				
				byte r=cha.get(rpos);
				if(r=='N'){
					//do nothing
				}else{
					md.append(count);
					md.append((char)cha.get(rpos));
					count=0;
				}
				
				rpos++;
//				cpos++;	
			
			}else if(m=='I' || m=='X' || m=='Y'){
//				cpos++;
//				count++;
			}else if(m=='D'){
//				if(prevM!='D'){
//					md.append(count);
//					count=0;
//					md.append('^');
//				}
//				md.append((char)cha.get(rpos));
				
				rpos++;
				dels++;
			}
			prevM=m;
		}
//		if(count>0){
			md.append(count);
//		}
		
		return md.toString();
	}
	
	
//	public static String makeMdTag(int chrom, int refstart, byte[] match, byte[] call, boolean colorspace){
//		if(match==null || chrom<0 || colorspace){return null;}
//		StringBuilder md=new StringBuilder(8);
//		md.append("MD:Z:");
//		
//		if(colorspace){throw new RuntimeException("Colorspace is incompatible with MD tags.");}
//		
//		ChromosomeArray cha=Data.getChromosome(chrom);
//		
//		byte prevM='?';
//		int count=0;
//		for(int mpos=0, cpos=0, rpos=refstart; mpos<match.length; mpos++){
//			byte c=call[cpos];
//			byte m=match[mpos];
//			
//			if(m=='m' || m=='s'){
//				count++;
//				rpos++;
//				cpos++;
//			}else if(m=='S'){
//				md.append(count);
//				md.append((char)cha.get(rpos));
//
//				count=0;
//				rpos++;
//				cpos++;
//			}else if(m=='N'){
//				
//				byte r=cha.get(rpos);
//				if(r=='N'){
//					//do nothing
//				}else{
//					md.append(count);
//					md.append((char)cha.get(rpos));
//					count=0;
//				}
//				
//				rpos++;
//				cpos++;	
//			
//			}else if(m=='I' || m=='X' || m=='Y'){
//				cpos++;
////				count++;
//			}else if(m=='D'){
//				if(prevM!='D'){
//					md.append(count);
//					count=0;
//					md.append('^');
//				}
//				md.append((char)cha.get(rpos));
//				rpos++;
//			}
//			prevM=m;
//		}
//		md.append(count);
//		
//		return md.toString();
//	}
	
	
	public SamLine(Read r, int fragNum){
		
		if(verbose){
			System.err.println("new SamLine for read with match "+(r.match==null ? "null" : new String(r.match)));
		}
		
		assert(!r.colorspace());
		Read r2=r.mate;
		final boolean perfect=r.perfect();
		
		if(Data.scaffoldLocs==null && r.obj!=null){
			assert(SET_FROM_OK) : "Sam format cannot be used as input to this program when no genome build is loaded.\n" +
					"Please index the reference first and rerun with e.g. 'build=1', or use a different input format.";
			setFrom((SamLine)r.obj);
			return;
		}
		
//		qname=r.id.replace(' ', '_').replace('\t', '_');
//		qname=r.id.split("\\s+")[0];
		qname=r.id.replace('\t', '_');
//		if(!KEEP_NAMES && qname.length()>2 && r2!=null){
//			if(qname.endsWith("/1") || qname.endsWith("/2") || qname.endsWith(" 1") || qname.endsWith(" 2")){}
//		}
		
		if(!KEEP_NAMES && qname.length()>2 && r2!=null){
			char c=qname.charAt(qname.length()-2);
			int num=(qname.charAt(qname.length()-1))-'1';
			if((num==0 || num==1) && (c==' ' || c=='/')){qname=qname.substring(0, qname.length()-2);}
//			if(r.pairnum()==num && (c==' ' || c=='/')){qname=qname.substring(0, qname.length()-2);}
		}
//		flag=Integer.parseInt(s[1]);
		
		int idx1=-1, idx2=-1;
		int chrom1=-1, chrom2=-1;
		int start1=-1, start2=-1, a1=0, a2=0;
		int stop1=-1, stop2=-1, b1=0, b2=0;
		int scaflen=0;
		byte[] name1=bytestar, name2=bytestar;
		if(r.mapped()){
			assert(r.chrom>=0);
			chrom1=r.chrom;
			start1=r.start;
			stop1=r.stop;
			if(Data.isSingleScaffold(chrom1, start1, stop1)){
				assert(Data.scaffoldLocs!=null) : "\n\n"+r+"\n\n"+r.obj+"\n\n";
				idx1=Data.scaffoldIndex(chrom1, (start1+stop1)/2);
				name1=Data.scaffoldNames[chrom1][idx1];
				scaflen=Data.scaffoldLengths[chrom1][idx1];
				a1=Data.scaffoldRelativeLoc(chrom1, start1, idx1);
				b1=a1-start1+stop1;
			}else{
				if(verbose){System.err.println("------------- Found multi-scaffold alignment! -------------");}
				r.setMapped(false);
				r.setPaired(false);
				r.match=null;
				if(r2!=null){r2.setPaired(false);}
			}
		}
		if(r2!=null && r2.mapped()){
			chrom2=r2.chrom;
			start2=r2.start;
			stop2=r2.stop;
			if(Data.isSingleScaffold(chrom2, start2, stop2)){
				idx2=Data.scaffoldIndex(chrom2, (start2+stop2)/2);
				name2=Data.scaffoldNames[chrom2][idx2];
				a2=Data.scaffoldRelativeLoc(chrom2, start2, idx2);
				b2=a2-start2+stop2;
			}else{
				if(verbose){System.err.println("------------- Found multi-scaffold alignment for r2! -------------");}
				r2.setMapped(false);
				r2.setPaired(false);
				r2.match=null;
				if(r!=null){r.setPaired(false);}
			}
		}
		
		flag=0;
		if(r2!=null){
			flag|=0x1;
			
			if(r.mapped() && r.valid() && r.match!=null &&
				(r2==null || (idx1==idx2 && r.paired() && r2.mapped() && r2.valid() && r2.match!=null))){flag|=0x2;}
			if(fragNum==0){flag|=0x40;}
			if(fragNum>0){flag|=0x80;}
		}
		if(!r.mapped()){flag|=0x4;}
		if(r2!=null && !r2.mapped()){flag|=0x8;}
		if(r.strand()==Gene.MINUS){flag|=0x10;}
		if(r2!=null && r2.strand()==Gene.MINUS){flag|=0x20;}
		if(r.secondary()){flag|=0x100;}
		if(r.discarded()){flag|=0x200;}
//		if(){flag|=0x400;}
		
//		assert(!r.secondary()) : r.mapScore;
		
//		2. FLAG: bitwise FLAG. Each bit is explained in the following table:
//		Bit Description
//		0x1 template having multiple fragments in sequencing
//		0x2 each fragment properly aligned according to the aligner
//		0x4 fragment unmapped
//		0x8 next fragment in the template unmapped
//		0x10 SEQ being reverse complemented
//		0x20 SEQ of the next fragment in the template being reversed
//		0x40 the first fragment in the template
//		0x80 the last fragment in the template
//		0x100 secondary alignment
//		0x200 not passing quality controls
//		0x400 PCR or optical duplicate
		
		rname=r.mapped() ? name1 : ((r2!=null && r2.mapped()) ? name2 : null);
		pos=r.mapped() ? a1+1 : ((r2!=null && r2.mapped()) ? Tools.max(a2+1, 1) : 0);
		
//		assert(false) : pos+", "+a1+", "+a2;
		
//		pos=Tools.max(pos, 1);
		
//		mapq=r.mapped() ? Data.max(1, -20+(int)(r.mapScore*60f/(100*r.mapLength))) : 0;//Scale of 0-40
		mapq=r.mapped() ? Data.max(1, r.mapScore/r.mapLength) : 0;//Scale of 0-100

		if(verbose){
			System.err.println("Making cigar for "+(r.match==null ? "null" : new String(r.match)));
		}
		
		if(r.bases!=null && r.mapped() && r.match!=null){
			final boolean inbounds=(a1>=0 && b1<scaflen);
			if(VERSION>1.3f){
				if(inbounds && perfect && !r.containsNonM()){//r.containsNonM() should be unnecessary...  it's there in case of clipping...
					cigar=(r.bases.length+"=");
//					System.err.println("SETTING cigar14="+cigar);
//					
//					byte[] match=r.match;
//					if(r.shortmatch()){match=Read.toLongMatchString(match);}
//					cigar=toCigar13(match, a1, b1, scaflen, r.bases);
//					System.err.println("RESETTING cigar14="+cigar+" from toCigar14("+new String(Read.toShortMatchString(match))+", "+a1+", "+b1+", "+scaflen+", "+r.bases+")");
				}else{
					byte[] match=r.match;
					if(r.shortmatch()){match=Read.toLongMatchString(match);}
					cigar=toCigar14(match, a1, b1, scaflen, r.bases);
//					System.err.println("CALLING toCigar14("+Read.toShortMatchString(match)+", "+a1+", "+b1+", "+scaflen+", "+r.bases+")");
				}
			}else{
				if(inbounds && (perfect || !r.containsNonNMS())){
					cigar=(r.bases.length+"M");
//					System.err.println("SETTING cigar13="+cigar);
//					
//					byte[] match=r.match;
//					if(r.shortmatch()){match=Read.toLongMatchString(match);}
//					cigar=toCigar13(match, a1, b1, scaflen, r.bases);
//					System.err.println("RESETTING cigar13="+cigar+" from toCigar13("+new String(Read.toShortMatchString(match))+", "+a1+", "+b1+", "+scaflen+", "+r.bases+")");
				}else{
					byte[] match=r.match;
					if(r.shortmatch()){match=Read.toLongMatchString(match);}
					cigar=toCigar13(match, a1, b1, scaflen, r.bases);
//					System.err.println("CALLING toCigar13("+Read.toShortMatchString(match)+", "+a1+", "+b1+", "+scaflen+", "+r.bases+")");
				}
			}
		}
		
		if(verbose){
			System.err.println("cigar="+cigar);
		}
		
//		assert(false);
		
//		assert(primary() || cigar.equals(stringstar)) : cigar;
		
		if(r.mapped()){
			int leadingClip=countLeadingClip(cigar);
			int clippedIndels=(r.match==null ? 0 : countLeadingIndels(a1, r.match));
			if(verbose){
				System.err.println("leadingClip="+leadingClip);
				System.err.println("clippedDels="+clippedIndels);
			}
			pos=pos+leadingClip+clippedIndels;
		}
		
		if((!r.mapped() || cigar==null || stringstar.equals(cigar)) && pos<1){
			//This is necessary to prevent unmapped reads from having negative POS,
			//and mapped reads without cigar strings from having POS less than 1.
			pos=Tools.max(pos, r.mapped() ? 1 : 0);
		}
		assert(pos>=0) : "Negative coordinate "+pos+" for read:\n\n"+r+"\n\n"+r2+"\n\n"+this+"\n\na1="+a1+", a2="+a2+", clip="+countLeadingClip(cigar);
//		if(pos<0){pos=0;cigar=null;rname=bytestar;mapq=0;flag|=0x4;}
		
//		assert(false) : "\npos="+pos+"\ncigar='"+cigar+"'\nVERSION="+VERSION+"\na1="+a1+", b1="+b1+"\n\n"+r.toString();
		
//		rnext=(r2==null ? stringstar : (r.mapped() && !r2.mapped()) ? "chr"+Gene.chromCodes[r.chrom] : "chr"+Gene.chromCodes[r2.chrom]);
		rnext=((r2==null || (!r.mapped() && !r2.mapped())) ? bytestar : (r.mapped() && r2.mapped()) ? (idx1==idx2 ? byteequals : name2) : byteequals);
		
		if(Data.scaffoldPrefixes){
			 if(rname!=null && rname!=bytestar){
				 int k=Tools.indexOf(rname, (byte)'$');
				 rname=Arrays.copyOfRange(rname, k+1, rname.length);
			 }
			 if(rnext!=null && rnext!=bytestar){
				 int k=Tools.indexOf(rnext, (byte)'$');
				 rnext=Arrays.copyOfRange(rnext, k+1, rnext.length);
			 }
		}
		
		if(r2==null){
			pnext=0;
		}else if(r2.mapped()){
			pnext=Tools.max(1, a2+1);
		}else if(r.mapped()){
			pnext=Tools.max(1, a1+1);
		}else{
			pnext=0;
		}
		tlen=(r2==null || !r.mapped() || !r2.mapped() || idx1!=idx2) ? 0 : 1+(Data.max(r.stop, r2.stop)-Data.min(r.start, r2.start));
//		if(r2==null || r.stop<=r2.start){
//			//plus sign
//		}else if(r2.stop<=r.start){
//			//minus sign
//			tlen=-tlen;
//		}else{
//			//They overlap... a lot.  Physically shorter than read length.
//			if(r.start<=r2.start){
//				
//			}else{
//				tlen=-tlen;
//			}
//		}
		//This version is less technically correct (does not account for very short insert reads) but probably more inline with what is expected
		if(r2==null || r.start<r2.start || (r.start==r2.start && r.pairnum()==0)){
			//plus sign
		}else{
			//minus sign
			tlen=-tlen;
		}
		
//		if(r.secondary()){
////			seq=qual=stringstar;
//			seq=qual=bytestar;
//		}else{
//			if(r.strand()==Gene.PLUS){
////				seq=new String(r.bases);
//				seq=r.bases.clone();
//				if(r.quality==null){
////					qual=stringstar;
//					qual=bytestar;
//				}else{
////					StringBuilder q=new StringBuilder(r.quality.length);
////					for(byte b : r.quality){
////						q.append((char)(b+33));
////					}
////					qual=q.toString();
//					qual=new byte[r.quality.length];
//					for(int i=0, j=qual.length-1; i<qual.length; i++, j--){
//						qual[i]=(byte)(r.quality[j]+33);
//					}
//				}
//			}else{
////				seq=new String(AminoAcid.reverseComplementBases(r.bases));
//				seq=AminoAcid.reverseComplementBases(r.bases);
//				if(r.quality==null){
////					qual=stringstar;
//					qual=bytestar;
//				}else{
////					StringBuilder q=new StringBuilder(r.quality.length);
////					for(int i=r.quality.length-1; i>=0; i--){
////						q.append((char)(r.quality[i]+33));
////					}
////					qual=q.toString();
//					qual=new byte[r.quality.length];
//					for(int i=0, j=qual.length-1; i<qual.length; i++, j--){
//						qual[i]=(byte)(r.quality[j]+33);
//					}
//				}
//			}
//		}
		
		if(r.secondary()){
//			seq=qual=bytestar;
			seq=qual=null;
		}else{
			seq=r.bases;
			if(r.quality==null){
//				qual=bytestar;
				qual=null;
			}else{
				qual=r.quality;
			}
		}
		
		
		if(r.mapped()){
			optional=new ArrayList<String>(8);

//			if(!r.secondary()){optional.add(r.ambiguous() ? "XT:A:R" : "XT:A:U");} //Not sure what do do for secondary alignments
			if(!r.secondary() && r.ambiguous()){optional.add("XT:A:R");} //Not sure what do do for secondary alignments
			
			int nm=r.bases.length;
			int dels=0;
			if(perfect){nm=0;}
			else if(r.match!=null){
				int delsCurrent=0;
				for(byte b : r.match){
					if(b=='m' || b=='C'){nm--;}
					if(b=='D'){delsCurrent++;}
					else{
						if(delsCurrent<=INTRON_LIMIT){dels+=delsCurrent;}
						delsCurrent=0;
					}
				}
				if(delsCurrent<=INTRON_LIMIT){dels+=delsCurrent;}
//				assert(false) : nm+",  "+dels+", "+delsCurrent+", "+r.bases.length+", "+r.match.length;
			}
			
			//Samtools puts nm tag in wrong place for deletions.  This "if" block is not really necessary except for exact text match to samtools.
			if(false && dels>0){
				if(MAKE_SM_TAG){optional.add("SM:i:"+mapq);}
//				optional.add("AM:i:"+Data.min(mapq, r2==null ? mapq : (r2.mapped() ? Data.max(1, -20+(int)((r2.mapScore*60f/(100*r2.mapLength)))) : 0)));
				optional.add("AM:i:"+Data.min(mapq, r2==null ? mapq : (r2.mapped() ? Data.max(1, r2.mapScore/r2.mapLength) : 0)));
				if(r.match!=null){optional.add("NM:i:"+(nm+dels));}
			}else{
				if(perfect){optional.add("NM:i:0");}
				else if(r.match!=null){optional.add("NM:i:"+(nm+dels));}
				if(MAKE_SM_TAG){optional.add("SM:i:"+mapq);}
//				optional.add("AM:i:"+Data.min(mapq, r2==null ? mapq : (r2.mapped() ? Data.max(1, -20+(int)((r2.mapScore*60f/(100*r2.mapLength)))) : 0)));
				optional.add("AM:i:"+Data.min(mapq, r2==null ? mapq : (r2.mapped() ? Data.max(1, r2.mapScore/r2.mapLength) : 0)));
			}
			
			if(MAKE_TOPHAT_TAGS){
				optional.add("AS:i:0");
				if(cigar==null || cigar.indexOf('N')<0){
					optional.add("XN:i:0");
				}else{
				}
				optional.add("XM:i:0");
				optional.add("XO:i:0");
				optional.add("XG:i:0");
				if(cigar==null || cigar.indexOf('N')<0){
					optional.add("YT:Z:UU");
				}else{
				}
				optional.add("NH:i:1");
			}else if(MAKE_XM_TAG){//XM tag.  For bowtie compatibility; unfortunately it is poorly defined.
				int x=0;
				if(r.discarded() || (!r.ambiguous() && !r.mapped())){
					x=0;//TODO: See if the flag needs to be present in this case.
				}else if(r.mapped()){
					x=1;
					if(r.numSites()>0 && r.numSites()>0){
						int z=r.topSite().score;
						for(int i=1; i<r.sites.size(); i++){
							SiteScore ss=r.sites.get(i);
							if(ss!=null && ss.score==z){x++;}
						}
					}
					if(r.ambiguous()){x=Tools.max(x, 2);}
				}
				if(x>=0){optional.add("XM:i:"+x);}
			}
			
			//XS tag
			if(MAKE_XS_TAG){
				String xs=makeXSTag(r);
				if(xs!=null){
					optional.add(xs);
					assert(r2==null || r.pairnum()!=r2.pairnum());
//					assert(r2==null || !r2.mapped() || r.strand()==r2.strand() || makeXSTag(r2)==xs) : 
//						"XS problem:\n"+r+"\n"+r2+"\n"+xs+"\n"+makeXSTag(r2)+"\n";
				}
			}
			
			if(MAKE_MD_TAG){
				String md=makeMdTag(r.chrom, r.start, r.match, r.bases, r.colorspace());
				if(md!=null){optional.add(md);}
			}
			
			if(r.mapped() && MAKE_NH_TAG){
				if(ReadStreamWriter.OUTPUT_SAM_SECONDARY_ALIGNMENTS && r.numSites()>1){
					optional.add("NH:i:"+r.sites.size());
				}else{
					optional.add("NH:i:1");
				}
			}
			
			if(MAKE_STOP_TAG && (perfect || r.match!=null)){optional.add(makeStopTag(pos, seq.length, cigar, perfect));}
			
			if(MAKE_IDENTITY_TAG && (perfect || r.match!=null)){optional.add(makeIdentityTag(r.match, perfect));}
			
			if(MAKE_INSERT_TAG && r2!=null){
				if(r.mapped() ||r.originalSite!=null){
					optional.add("X8:Z:"+r.insertSizeMapped(false)+(r.originalSite==null ? "" : ","+r.insertSizeOriginalSite()));
//					assert(r.originalSite==null || r2.originalSite==null || r.insertSizeOriginalSite()==r2.insertSizeOriginalSite());
//					assert(!r.mapped() || !r2.mapped() || r.insertSizeMapped()==r2.insertSizeMapped());
				}
			}
			if(MAKE_CORRECTNESS_TAG){
				final SiteScore ss0=r.originalSite;
				if(ss0!=null){
					optional.add("X9:Z:"+(ss0.isCorrect(r.chrom, r.strand(), r.start, r.stop, 0) ? "T" : "F"));
				}
			}
			
			if(MAKE_CUSTOM_TAGS){
				int sites=r.numSites() + (r.originalSite==null ? 0 : 1);
				if(sites>0){
					StringBuilder sb=new StringBuilder();
					sb.append("X1:Z:");
					if(r.sites!=null){
						for(SiteScore ss : r.sites){
							sb.append('$');
							sb.append(ss.toText());
						}
					}
					if(r.originalSite!=null){
						sb.append('$');
						sb.append('*');
						sb.append(r.originalSite.toText());
					}
					optional.add(sb.toString());
				}

				if(r.match!=null){
					byte[] match=r.match;
					if(!r.shortmatch()){
						match=Read.toShortMatchString(match);
					}
					optional.add("X2:Z:"+new String(match));
				}

				optional.add("X3:i:"+r.mapScore);
				optional.add("X4:i:"+r.mapLength);
				optional.add("X5:Z:"+r.numericID);
				optional.add("X6:i:"+(r.flags|(r.match==null ? 0 : Read.SHORTMATCHMASK)));
				if(r.copies>1){optional.add("X7:i:"+r.copies);}
			}
			
		}
//		assert(r.pairnum()==1) : "\n"+r.toText(false)+"\n"+this+"\n"+r2;
	}
	
	public SamLine(String[] s){
		assert(!s[0].startsWith("@")) : "Tried to make a SamLine from a header: "+s[0];
		assert(s.length>=11) : "\nNot all required fields are present: "+s.length+"\nline='"+Arrays.toString(s)+"'\n";
		if(s.length<11){
			System.err.println("Invalid SamLine: "+Arrays.toString(s));
			return;
		}
		qname=s[0];
		flag=Integer.parseInt(s[1]);
		rname=s[2].getBytes();
		pos=Integer.parseInt(s[3]);
//		try {
//			Integer.parseInt(s[4]);
//		} catch (NumberFormatException e) {
//			System.err.println(Arrays.toString(s));
//		}
		mapq=Character.isDigit(s[4].charAt(0)) ? Integer.parseInt(s[4]) : 99; //Added for non-compliant mappers that put * here
		cigar=s[5];
		rnext=s[6].getBytes();
		pnext=Integer.parseInt(s[7]);
		tlen=Character.isDigit(s[8].charAt(0)) ? Integer.parseInt(s[8]) : 0; //Added for non-compliant mappers that put * here
//		seq=s[9];
//		qual=s[10];
		seq=(s[9].equals(stringstar) ? null : s[9].getBytes());
		qual=(s[10].equals(stringstar) ? null : s[10].getBytes());
		
		if(mapped() && strand()==Gene.MINUS){
			if(seq!=bytestar){AminoAcid.reverseComplementBasesInPlace(seq);}
			if(qual!=bytestar){Tools.reverseInPlace(qual);}
		}
		
		if(qual!=null && qual!=bytestar){
			for(int i=0; i<qual.length; i++){qual[i]-=33;}
		}
		
		if(s.length>11){
			optional=new ArrayList<String>(s.length-11);
			for(int i=11; i<s.length; i++){
				optional.add(s[i]);
			}
		}
	}
	
	public SamLine(byte[] s){
		assert(s[0]!='@') : "Tried to make a SamLine from a header: "+new String(s);
		
		int a=0, b=0;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(s);
		qname=(b==a+1 && s[a]=='*' ? null : new String(s, a, b-a));
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(s);
		flag=Tools.parseInt(s, a, b);
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(s);
		rname=(b==a+1 && s[a]=='*' ? null : Arrays.copyOfRange(s, a, b));
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(s);
		pos=Tools.parseInt(s, a, b);
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(s);
		mapq=Tools.parseInt(s, a, b);
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 5: "+new String(s);
		cigar=(b==a+1 && s[a]=='*' ? null : new String(s, a, b-a));
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 6: "+new String(s);
		rnext=(b==a+1 && s[a]=='*' ? null : Arrays.copyOfRange(s, a, b));
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 7: "+new String(s);
		pnext=Tools.parseInt(s, a, b);
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 8: "+new String(s);
		tlen=Tools.parseInt(s, a, b);
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 9: "+new String(s);
//		seq=new String(s, a, b-a);
		seq=(b==a+1 && s[a]=='*' ? null : Arrays.copyOfRange(s, a, b));
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 10: "+new String(s);
//		qual=new String(s, a, b-a);
		qual=(b==a+1 && s[a]=='*' ? null : Arrays.copyOfRange(s, a, b));
		b++;
		a=b;

		assert((seq==bytestar)==(Tools.equals(seq, bytestar)));
		assert((qual==bytestar)==(Tools.equals(qual, bytestar)));
		
		if(mapped() && strand()==Gene.MINUS){
			if(seq!=bytestar){AminoAcid.reverseComplementBasesInPlace(seq);}
			if(qual!=bytestar){Tools.reverseInPlace(qual);}
		}
		
		if(qual!=null && qual!=bytestar){
			for(int i=0; i<qual.length; i++){qual[i]-=33;}
		}
		
		if(b<s.length){
			optional=new ArrayList<String>(4);
			while(b<s.length){
				while(b<s.length && s[b]!='\t'){b++;}
				if(b>a){
					String x=new String(s, a, b-a);
					optional.add(x);
				}else{
					//Empty field
				}
				b++;
				a=b;
			}
		}
	}
	
	
	public Read parseName(){
		try {
			String[] answer=qname.split("_");
			long id=Long.parseLong(answer[0]);
			byte trueChrom=Gene.toChromosome(answer[1]);
			byte trueStrand=Byte.parseByte(answer[2]);
			int trueLoc=Integer.parseInt(answer[3]);
			int trueStop=Integer.parseInt(answer[4]);
//			byte[] quals=qual.getBytes();
//			byte[] quals=qual;
//			for(int i=0; i<quals.length; i++){quals[i]-=33;}
//			Read r=new Read(seq.getBytes(), trueChrom, trueStrand, trueLoc, trueStop, qname, quals, false, id);
			Read r=new Read(seq, trueChrom, trueStrand, trueLoc, trueStop, qname, qual, false, id);
			return r;
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}
	
	public long parseNumericId(){
		return Long.parseLong(qname.substring(0, qname.indexOf('_')));
	}
	
	public Read toRead(boolean parseCustom){
		
		SiteScore originalSite=null;
		long numericId_=0;
		boolean synthetic=false;
		
		if(parseCustom){
			try {
				String[] answer=qname.split("_");
				numericId_=Long.parseLong(answer[0]);
				byte trueChrom=Gene.toChromosome(answer[1]);
				byte trueStrand=Byte.parseByte(answer[2]);
				int trueLoc=Integer.parseInt(answer[3]);
				int trueStop=Integer.parseInt(answer[4]);
				
				originalSite=new SiteScore(trueChrom, trueStrand, trueLoc, trueStop, 0, 0);
				synthetic=true;
				
			} catch (NumberFormatException e) {
				System.err.println("Failed to parse "+qname);
			} catch (NullPointerException e) {
				System.err.println("Bad read with no name.");
				return null;
			}
		}
//		assert(false) : originalSite;
		
		
		if(Data.GENOME_BUILD>=0){
			
		}
		
		int chrom_=-1; 
		byte strand_=strand();
		int start_=start();
		int stop_=stop();
		assert(start_<=stop_) : start_+", "+stop_;
		boolean cs_=colorspace();
		
		if(Data.GENOME_BUILD>=0 && rname!=null && (rname.length!=1 || rname[0]!='*')){
			ScafLoc sc=Data.getScafLoc(rname);
			assert(sc!=null) : "Can't find scaffold in reference with name "+new String(rname)+"\n"+this;
			if(sc!=null){
				chrom_=sc.chrom;
				start_+=sc.loc;
				stop_+=sc.loc;
			}
		}

////		byte[] quals=(qual==null || (qual.length()==1 && qual.charAt(0)=='*')) ? null : qual.getBytes();
////		byte[] quals=(qual==null || (qual.length==1 && qual[0]=='*')) ? null : qual.clone();
//		byte[] quals=(qual==null || (qual.length==1 && qual[0]=='*')) ? null : qual;
//		byte[] bases=seq==null ? null : seq.clone();
//		if(strand_==Gene.MINUS){//Minus-mapped SAM lines have bases and quals reversed
//			AminoAcid.reverseComplementBasesInPlace(bases);
//			Tools.reverseInPlace(quals);
//		}
//		Read r=new Read(bases, chrom_, strand_, start_, stop_, qname, quals, cs_, numericId_);
		
		final Read r;
		{
			byte[] seqX=(seq==null || (seq.length==1 && seq[0]=='*')) ? null : seq;
			byte[] qualX=(qual==null || (qual.length==1 && qual[0]=='*')) ? null : qual;
			String qnameX=(qname==null || qname.equals(stringstar)) ? null : qname;
			r=new Read(seqX, chrom_, strand_, start_, stop_, qnameX, qualX, cs_, numericId_);
		}
		
		r.setMapped(mapped());
		r.setSynthetic(synthetic);
//		r.setPairnum(pairnum()); //TODO:  Enable after fixing assertions that this will break in read input streams.
		if(originalSite!=null){
			r.originalSite=originalSite;
		}
		
		r.mapScore=mapq;
		r.setSecondary(!primary());
		
//		if(mapped()){
//			r.list=new ArrayList<SiteScore>(1);
//			r.list.add(new SiteScore(r.chrom, r.strand(), r.start, r.stop, 0));
//		}
		
//		System.out.println(optional);
		if(optional!=null){
			for(String s : optional){
				if(s.equals("XT:A:R")){
					r.setAmbiguous(true);
				}else if(s.startsWith("X1:Z:")){
//					System.err.println("Found X1 tag!\t"+s);
					String[] split=s.split("\\$");
//					assert(false) : Arrays.toString(split);
					ArrayList<SiteScore> list=new ArrayList<SiteScore>(3);
					
					for(int i=1; i<split.length; i++){
//						System.err.println("Processing ss\t"+split[i]);
						String s2=split[i];
						SiteScore ss=SiteScore.fromText(s2);
						if(s2.charAt(0)=='*'){
							r.originalSite=ss;
						}else{
							list.add(ss);
						}
					}
//					System.err.println("List size = "+list.size());
					if(list.size()>0){r.sites=list;}
				}else if(s.startsWith("X2:Z:")){
					String s2=s.substring(5);
					r.match=s2.getBytes();
				}else if(s.startsWith("X3:i:")){
					String s2=s.substring(5);
//					r.mapScore=Integer.parseInt(s2); //Messes up generation of ROC curve
				}else if(s.startsWith("X4:i:")){
					String s2=s.substring(5);
					r.mapLength=Integer.parseInt(s2);
				}else if(s.startsWith("X5:Z:")){
					String s2=s.substring(5);
					r.numericID=Long.parseLong(s2);
				}else if(s.startsWith("X6:i:")){
					String s2=s.substring(5);
					r.flags=Integer.parseInt(s2);
				}else if(s.startsWith("X7:i:")){
					String s2=s.substring(5);
					r.copies=Integer.parseInt(s2);
				}else{
//					System.err.println("Unknown SAM field:"+s);
				}
			}
		}
		
		if(r.match==null && cigar!=null && (CONVERT_CIGAR_TO_MATCH || cigar.indexOf('M')<0)){
			r.match=cigarToShortMatch(cigar, true);
			
			if(r.match!=null){
				r.setShortMatch(true);
				if(Tools.indexOf(r.match, (byte)'B')>=0){
					boolean success=r.fixMatchB();
//					if(!success){r.match=null;}
//					assert(false) : new String(r.match);
				}
//				assert(false) : new String(r.match);
			}
//			assert(false) : new String(r.match);
//			System.err.println(">\n"+cigar+"\n"+(r.match==null ? "null" : new String(r.match)));
		}
//		assert(false) : new String(r.match);
		
//		System.err.println("Resulting read: "+r.toText());
		
		return r;
		
	}
	
	/** Aproximate length of result of SamLine.toText() */
	public int textLength(){
		int len=11; //11 tabs
		len+=(3+9+3+9);
		len+=(tlen>999 ? 9 : 3);
		
		len+=(qname==null ? 1 : qname.length());
		len+=(rname==null ? 1 : rname.length);
		len+=(rnext==null ? 1 : rnext.length);
		len+=(cigar==null ? 1 : cigar.length());
		len+=(seq==null ? 1 : seq.length);
		len+=(qual==null ? 1 : qual.length);
		
		if(optional!=null){
			len+=optional.size();
			for(String s : optional){len+=s.length();}
		}
		return len;
	}
	
	public ByteBuilder toBytes(ByteBuilder bb){
		
		final int buflen=Tools.max((rname==null ? 1 : rname.length), (rnext==null ? 1 : rnext.length), (seq==null ? 1 : seq.length), (qual==null ? 1 : qual.length));
		
		if(bb==null){bb=new ByteBuilder(textLength()+4);}
		if(qname==null){bb.append('*').append('\t');}else{bb.append(qname).append('\t');}
		bb.append(flag).append('\t');
		append(bb, rname).append('\t');
		bb.append(pos).append('\t');
		bb.append(mapq).append('\t');
		if(cigar==null){bb.append('*').append('\t');}else{bb.append(cigar).append('\t');}
		append(bb, rnext).append('\t');
		bb.append(pnext).append('\t');
		bb.append(tlen).append('\t');
		
		if(mapped() && strand()==Gene.MINUS){
			appendReverseComplimented(bb, seq).append('\t');
			appendQualReversed(bb, qual);
		}else{
			append(bb, seq).append('\t');
			appendQual(bb, qual);
		}

//		assert(seq.getClass()==String.class);
//		assert(qual.getClass()==String.class);
//		sb.append(seq).append('\t');
//		sb.append(qual);
		
		if(optional!=null){
			for(String s : optional){
				bb.append('\t').append(s);
			}
		}
		return bb;
	}
	
	public StringBuilder toText(){
		
		final int buflen=Tools.max((rname==null ? 1 : rname.length), (rnext==null ? 1 : rnext.length), (seq==null ? 1 : seq.length), (qual==null ? 1 : qual.length));
		final char[] buffer=Shared.getTLCB(buflen);
		
		StringBuilder sb=new StringBuilder(textLength()+4);
		if(qname==null){sb.append('*').append('\t');}else{sb.append(qname).append('\t');}
		sb.append(flag).append('\t');
		append(sb, rname, buffer).append('\t');
		sb.append(pos).append('\t');
		sb.append(mapq).append('\t');
		if(cigar==null){sb.append('*').append('\t');}else{sb.append(cigar).append('\t');}
		append(sb, rnext, buffer).append('\t');
		sb.append(pnext).append('\t');
		sb.append(tlen).append('\t');
		
		if(mapped() && strand()==Gene.MINUS){
			appendReverseComplimented(sb, seq, buffer).append('\t');
			appendQualReversed(sb, qual, buffer);
		}else{
			append(sb, seq, buffer).append('\t');
			appendQual(sb, qual, buffer);
		}

//		assert(seq.getClass()==String.class);
//		assert(qual.getClass()==String.class);
//		sb.append(seq).append('\t');
//		sb.append(qual);
		
		if(optional!=null){
			for(String s : optional){
				sb.append('\t').append(s);
			}
		}
		return sb;
	}
	
	public String toString(){return toText().toString();}
	
//	Bit Description
//	0x1 template having multiple fragments in sequencing
//	0x2 each fragment properly aligned according to the aligner
//	0x4 fragment unmapped
//	0x8 next fragment in the template unmapped
//	0x10 SEQ being reverse complemented
//	0x20 SEQ of the next fragment in the template being reversed
//	0x40 the first fragment in the template
//	0x80 the last fragment in the template
//	0x100 secondary alignment
//	0x200 not passing quality controls
//	0x400 PCR or optical duplicate
	
	public boolean hasMate(){
		return (flag&0x1)==0x1;
	}
	
	public boolean properPair(){
		return (flag&0x2)==0x2;
	}
	
	public boolean mapped(){
		return (flag&0x4)!=0x4;
//		0x4 fragment unmapped
//		0x8 next fragment in the template unmapped
	}
	
	public boolean nextMapped(){
		return (flag&0x8)!=0x8;
//		0x4 fragment unmapped
//		0x8 next fragment in the template unmapped
	}

	public byte strand(){
		return ((flag&0x10)==0x10 ? (byte)1 : (byte)0);
	}
	
	public byte nextStrand(){
		return ((flag&0x20)==0x20 ? (byte)1 : (byte)0);
	}
	
	public boolean firstFragment(){
		return (flag&0x40)==0x40;
	}
	
	public boolean lastFragment(){
		return (flag&0x80)==0x80;
	}
	
	public int pairnum(){
		return firstFragment() ? 0 : lastFragment() ? 1 : 0;
	}

	public boolean primary(){return (flag&0x100)==0;}
	public void setPrimary(boolean b){
		if(b){
			flag=flag|0x100;
		}else{
			flag=flag&~0x100;
		}
	}
	
	public boolean discarded(){
		return (flag&0x200)==0x200;
	}
	
//	/** Assumes rname is an integer. */
//	public int chrom(){
//		if(Data.GENOME_BUILD<0){return -1;}
//		HashMap sc
//	}
	
	/** Assumes rname is an integer. */
	public int chrom_old(){
		assert(false);
		if(!Character.isDigit(rname[0]) && !Character.isDigit(rname[rname.length-1])){
			if(warning){
				warning=false;
				System.err.println("Warning - sam lines need a chrom field.");
			}
			return -1;
		}
		assert(Shared.anomaly || '*'==rname[0] || (Character.isDigit(rname[0]) && Character.isDigit(rname[rname.length-1]))) : 
			"This is no longer correct, considering that sam lines are named by scaffold.  They need a chrom field.\n"+new String(rname);
		if(rname==null || Arrays.equals(rname, bytestar) || !(Character.isDigit(rname[0]) && Character.isDigit(rname[rname.length-1]))){return -1;}
		//return Gene.toChromosome(new String(rname));
		//return Integer.parseInt(new String(rname)));
		final byte z='0';
		int x=rname[0]-z;
		for(int i=1; i<rname.length; i++){
			x=(x*10)+(rname[i]-z);
		}
		return x;
	}
	
	public int start(){
		int x=countLeadingClip(cigar);
		return pos-1-x;
	}
	
	public int stop(){
		return stop(start());
	}
	
	public int stop(int start){
		if(!mapped() || cigar==null || cigar.charAt(0)=='*'){
//			return -1;
			return start+(seq==null ? 0 : Tools.max(0, seq.length-1));
		}
		int r=start+calcCigarLength(cigar)-1;
//		System.err.println("start= "+start+", stop="+r);
		return r;
	}
	
	public int stop2(){
		if(mapped() && cigar!=null && cigar.charAt(0)!='*'){return stop();}
//		return (seq==null ? -1 : start()+seq.length());
		return (seq==null ? -1 : start()+seq.length);
	}
	
	public long numericId(){
		return 0;
	}
	
	public boolean colorspace(){
//		for(int i=0; i<seq.length(); i++){
//			char c=seq.charAt(i);
		for(int i=0; i<seq.length; i++){
			char c=(char)seq[i];
			if(c=='0' || c=='1' || c=='2' || c=='3'){return true;}
			if(c=='A' || c=='C' || c=='G' || c=='T'){return false;}
			if(Character.isLetter(c)){
				if(c!='N'){return false;}
			}else{
				assert(false) : "\n"+c+"\n"+this.toText()+"\n";
			}
		}
		return false; //Can't tell...  but it doesn't matter.
	}
	
	public boolean pairedOnSameChrom(){
		return Tools.equals(rnext, byteequals) || Arrays.equals(rname, rnext);
	}
	
	public static int calcCigarLength(String cigar){
		if(cigar==null){return 0;}
		int len=0;
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='M' || c=='D' || c=='N' || c=='S' || c=='X' || c=='='){
					len+=current;
				}else if (c=='H'){ 
					//In this case, the base string is the wrong length since letters were truncated.
					//Therefore, the bases cannot be used for calling variations after mapping.
					//Hard clipping messes up original location verification.
					//Therefore...  len+=current would be best in practice, but for GRADING purposes, leaving it disabled is best.

					//len+=current;
				}else if(c=='I'){
					//do nothing
				}else if(c=='P'){
					throw new RuntimeException("Unhandled cigar symbol: "+c+"\n"+cigar+"\n");
					//'P' is currently poorly defined
				}else{
					throw new RuntimeException("Unhandled cigar symbol: "+c+"\n"+cigar+"\n");
				}
				current=0;
			}
		}
		return len;
	}
	
	/** Length of clipped initial bases.  Used to calculate correct start location of clipped reads. */
	public static int countLeadingClip(String cigar){
		if(cigar==null){return 0;}
		int len=0;
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Character.isLetter(c) || c=='='){
				if(c=='H' || (SUBTRACT_LEADING_SOFT_CLIP && c=='S')){
					len+=current;
				}else{
					break;
				}
				current=0;
			}else{
				current=(current*10)+(c-'0');
			}
		}
		return len;
	}
	
	/** Length of clipped (out of bounds) initial insertions and deletions. */
	public static int countLeadingIndels(int rloc, byte[] match){
		if(match==null || rloc>=0){return 0;}
		int dels=0;
		int inss=0;
		int cloc=0;
		for(int mloc=0; mloc<match.length && rloc<0; mloc++){
			byte b=match[mloc];
			assert(!Character.isDigit(b));
			if(b=='D'){
				dels++;
				rloc++;
			}else if(b=='I'){
				inss++;
				cloc++;
			}else{
				rloc++;
				cloc++;
			}
		}
		return dels-inss;
	}
	
	/**
	 * @param cigar
	 * @return Max consecutive match, sub, del, ins, or clip symbols
	 */
	public static final int[] cigarToMdsiMax(String cigar) {
		if(cigar==null){return null;}
		int[] msdic=new int[5];
		
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='M' || c=='='){
					msdic[0]=Tools.max(msdic[0], current);
				}else if(c=='X'){
					msdic[1]=Tools.max(msdic[1], current);
				}else if(c=='D' || c=='N'){
					msdic[2]=Tools.max(msdic[2], current);
				}else if(c=='I'){
					msdic[3]=Tools.max(msdic[3], current);
				}else if(c=='S' || c=='H' || c=='P'){
					msdic[4]=Tools.max(msdic[4], current);
				}
				current=0;
			}
		}
		return msdic;
	}
	
	/**
	 * @param cigar
	 * @return Total number of match, sub, del, ins, or clip symbols
	 */
	public static final int[] cigarToMsdic(String cigar) {
		if(cigar==null){return null;}
		int[] msdic=new int[5];
		
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='M' || c=='='){
					msdic[0]+=current;
				}else if(c=='X'){
					msdic[1]+=current;
				}else if(c=='D' || c=='N'){
					msdic[2]+=current;
				}else if(c=='I'){
					msdic[3]+=current;
				}else if(c=='S' || c=='H' || c=='P'){
					msdic[4]+=current;
				}
				current=0;
			}
		}
		return msdic;
	}
	
	/**
	 * @param cigar
	 * @return Match string of this cigar string when possible, otherwise null
	 */
	public static final byte[] cigarToShortMatch(String cigar, boolean allowM) {
		if(cigar==null || cigar.equals(stringstar)){return null;}
		
		int total=0;
		int current=0;
		
//		int totalLen=0;
//		int currentLen=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				
				if(c=='M'){
					if(!allowM){return null;} //Loss of information
				}else if(c=='H'){
					current=0; //Information destroyed
				}else if(c=='P'){
					return null; //Undefined symbol
				}

				total+=current;
				current=0;
			}
		}
		
		if(total<1){return null;}
		
		StringBuilder sb=new StringBuilder(cigar.length());
		
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='='){
					sb.append('m');
					if(current>1){sb.append(current);}
				}else if(c=='X'){
					sb.append('S');
					if(current>1){sb.append(current);}
				}else if(c=='D' || c=='N'){
					sb.append('D');
					if(current>1){sb.append(current);}
				}else if(c=='I'){
					sb.append('I');
					if(current>1){sb.append(current);}
				}else if(c=='S'){
					sb.append('C');
					if(current>1){sb.append(current);}
				}else if(c=='M'){
					sb.append('B');
					if(current>1){sb.append(current);}
				}
				current=0;
			}
		}
		
		byte[] match=new byte[sb.length()];
		for(int i=0; i<match.length; i++){
			match[i]=(byte)sb.charAt(i);
		}
		return match;
	}
	
	private static StringBuilder append(StringBuilder sb, byte[] a, char[] buffer){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		{//This is actually faster
			assert(buffer.length>=a.length);
			for(int i=0; i<a.length; i++){
				buffer[i]=(char)a[i];
			}
			sb.append(buffer, 0, a.length);
		}
//		for(byte b : a){
//			sb.append((char)b);
//		}
		return sb;
	}
	
	private static StringBuilder appendReverseComplimented(StringBuilder sb, byte[] a, char[] buffer){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		{//This is actually faster
			assert(buffer.length>=a.length);
			for(int i=0, j=a.length-1; j>=0; i++, j--){buffer[i]=(char)AminoAcid.baseToComplementExtended[a[j]];}
			sb.append(buffer, 0, a.length);
		}
//		for(int i=a.length-1; i>=0; i--){
//			sb.append((char)AminoAcid.baseToComplementEbuffertended[a[i]]);
//		}
		return sb;
	}
	
	private static StringBuilder appendQual(StringBuilder sb, byte[] a, char[] buffer){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		{//This is actually faster
			assert(buffer.length>=a.length);
			for(int i=0; i<a.length; i++){buffer[i]=(char)(a[i]+33);}
			sb.append(buffer, 0, a.length);
		}
//		for(byte b : a){
//			sb.append((char)(b+33));
//		}
		return sb;
	}
	
	private static StringBuilder appendQualReversed(StringBuilder sb, byte[] a, char[] buffer){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		{//This is actually faster
			assert(buffer.length>=a.length);
			for(int i=0, j=a.length-1; j>=0; i++, j--){buffer[i]=(char)(a[j]+33);}
			sb.append(buffer, 0, a.length);
		}
//		for(int i=a.length-1; i>=0; i--){
//			sb.append((char)(a[i]+33));
//		}
		return sb;
	}
	
	private static ByteBuilder append(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		return sb.append(a);
	}
	
	private static ByteBuilder appendReverseComplimented(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}

		sb.ensureExtra(a.length);
		byte[] buffer=sb.array;
		int i=sb.length;
		for(int j=a.length-1; j>=0; i++, j--){buffer[i]=AminoAcid.baseToComplementExtended[a[j]];}
		sb.length+=a.length;

		return sb;
	}
	
	private static ByteBuilder appendQual(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}

		sb.ensureExtra(a.length);
		byte[] buffer=sb.array;
		int i=sb.length;
		for(int j=0; j<a.length; i++, j++){buffer[i]=(byte)(a[j]+33);}
		sb.length+=a.length;

		return sb;
	}
	
	private static ByteBuilder appendQualReversed(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		
		sb.ensureExtra(a.length);
		byte[] buffer=sb.array;
		int i=sb.length;
		for(int j=a.length-1; j>=0; i++, j--){buffer[i]=(byte)(a[j]+33);}
		sb.length+=a.length;
		
		return sb;
	}
	
	/** Assumes a custom name including original location */
	public byte[] originalContig(){
//		assert(PARSE_CUSTOM);
		int loc=-1;
		int count=0;
		for(int i=0; i<qname.length() && loc==-1; i++){
			if(qname.charAt(i)=='_'){
				count++;
				if(count==6){loc=i;}
			}
		}
		if(loc==-1){
			return null;
		}
		return qname.substring(loc+1).getBytes();
	}
	
	/** Assumes a custom name including original location */
	public int originalContigStart(){
//		assert(PARSE_CUSTOM);
		int loc=-1;
		int count=0;
		for(int i=0; i<qname.length() && loc==-1; i++){
			if(qname.charAt(i)=='_'){
				count++;
				if(count==5){loc=i;}
			}
		}
		if(loc==-1){
			return -1;
		}
		
		int sum=0;
		int mult=1;
		for(int i=loc+1; i<qname.length(); i++){
			char c=qname.charAt(i);
			if(!Character.isDigit(c)){
				if(i==loc+1 && c=='-'){mult=-1;}
				else{break;}
			}else{
				sum=(sum*10)+(c-'0');
			}
		}
		return sum*mult;
	}
	
	public String matchTag(){
		if(optional==null){return null;}
		for(String s : optional){
			if(s.startsWith("X2:Z:")){
				return s;
			}
		}
		return null;
	}
	
	private void setFrom(SamLine sl){
		qname=sl.qname;
		flag=sl.flag;
		rname=sl.rname;
		pos=sl.pos;
		mapq=sl.mapq;
		cigar=sl.cigar;
		rnext=sl.rnext;
		pnext=sl.pnext;
		tlen=sl.tlen;
		seq=sl.seq;
		qual=sl.qual;
		optional=sl.optional;
	}
	
	private String makeXSTag(Read r){
		if(r.mapped() && cigar!=null && cigar.indexOf('N')>=0){
//			System.err.println("For read "+r.pairnum()+" mapped to strand "+r.strand());
			boolean plus=(r.strand()==Gene.PLUS); //Assumes secondstrand=false
//			System.err.println("plus="+plus);
			if(r.pairnum()!=0){plus=!plus;}
//			System.err.println("plus="+plus);
			if(XS_SECONDSTRAND){plus=!plus;}
//			System.err.println("plus="+plus);
			return (plus ? XSPLUS : XSMINUS);
		}else{
			return null;
		}
	}
	
	public byte[] rname(){return rname;}
	public byte[] rnext(){return rnext;}
	
	public String qname;
	public int flag;
	private byte[] rname;
	public int pos;
	public int mapq;
	public String cigar;
	private byte[] rnext;
	public int pnext;
	public int tlen;
	public byte[] seq;
	public byte[] qual;
	public ArrayList<String> optional;
	
	public Object obj;
	
	/** Turn this off for RNAseq */
	public static boolean MAKE_MD_TAG=false;
	public static boolean MAKE_SM_TAG=false;
	public static boolean MAKE_XM_TAG=false;
	public static boolean MAKE_XS_TAG=false;
	public static boolean MAKE_AS_TAG=false; //TODO: Alignment score from aligner
	public static boolean MAKE_NH_TAG=false;
	public static boolean MAKE_TOPHAT_TAGS=false;
	public static boolean XS_SECONDSTRAND=false;
	public static boolean MAKE_IDENTITY_TAG=false;
	public static boolean MAKE_STOP_TAG=false;
	public static boolean MAKE_CUSTOM_TAGS=false;
	public static boolean MAKE_INSERT_TAG=false;
	public static boolean MAKE_CORRECTNESS_TAG=false;
	public static boolean CONVERT_CIGAR_TO_MATCH=false;
	public static boolean SOFT_CLIP=true;
	/** OK to use the "setFrom" function which uses the old SamLine instead of translating the read, if a genome is not loaded. Should be false when processing occurs. */
	public static boolean SET_FROM_OK=false;
	/** For paired reads, keep original names rather than changing read2's name to match read1 */
	public static boolean KEEP_NAMES=false;
	public static float VERSION=1.3f;
	/** Tells program when to use 'N' rather than 'D' in cigar strings */
	public static int INTRON_LIMIT=Integer.MAX_VALUE;

	private static boolean warning=System.getProperty("user.dir").contains("/bushnell/");
	
	/** SSAHA2 incorrectly calculates the start position of reads with soft-clipped starts, and needs this enabled. */
	public static boolean SUBTRACT_LEADING_SOFT_CLIP=true;
	/** Sort header scaffolds in alphabetical order to be more compatible with Tophat */
	public static boolean SORT_SCAFFOLDS=false;
	public static boolean verbose=false;
	
	private static final String stringstar="*";
	private static final byte[] bytestar=new byte[] {(byte)'*'};
	private static final byte[] byteequals=new byte[] {(byte)'='};
	private static final String XSPLUS="XS:A:+", XSMINUS="XS:A:-";
	
}
