package stream;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;

import align2.Shared;
import align2.Tools;

import dna.Data;
import dna.Gene;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.TextFile;


public class FASTQ {
	
	public static void writeFASTQ(Read[] reads, String fname){
		StringBuilder sb=new StringBuilder();
		for(Read r : reads){
			String[] quad=toFASTQ(r);
			for(int i=0; i<quad.length; i++){
				sb.append(quad[i]);
				sb.append('\n');
			}
		}
		ReadWrite.writeString(sb, fname);
	}
	
//	public static boolean isInterleaved(String fname){
//		if(!TEST_INTERLEAVED && !FORCE_INTERLEAVED){return false;}
//		assert(tf.is!=System.in && !tf.name.equals("stdin") && !tf.name.startsWith("stdin."));
//		if(tf.is!=System.in && !tf.name.equals("stdin") && !tf.name.startsWith("stdin.")){return FORCE_INTERLEAVED;}
//		String s=null;
//		
//		String[] oct=new String[8];
//	}
	
//	public static boolean isInterleaved_old(String fname){
////		assert(false) : TEST_INTERLEAVED+", "+FORCE_INTERLEAVED;
//		if(!TEST_INTERLEAVED && !FORCE_INTERLEAVED){
//			testQuality(fname);
//			return false;
//		}
//		assert(!fname.equals("stdin") && !fname.startsWith("stdin."));
//		if(fname.equals("stdin") || fname.startsWith("stdin.")){return FORCE_INTERLEAVED;}
////		TextFile tf=new TextFile(fname);
////		assert(tf.is!=System.in);
////		if(tf.is==System.in){return FORCE_INTERLEAVED;}
//		
//		String[] oct=new String[8];
//		int cntr=0;
//		
//		{
//			InputStream is=ReadWrite.getInputStream(fname, false, false);
//			BufferedReader br=new BufferedReader(new InputStreamReader(is));
//			try {
//				for(String s=br.readLine(); s!=null && cntr<8; s=br.readLine()){
//					oct[cntr]=s;
//					cntr++;
//				}
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
//		
//		if(oct[0]==null){return false;}
//		
//		testQuality(oct);
//		
//		if(cntr<8){return false;}
////		assert(false);
//		assert(oct[0]==null || oct[0].startsWith("@")) : "Does not appear to be a valid FASTQ file:\n"+new String(oct[0]);
//		assert(oct[2]==null || oct[2].startsWith("+")) : "Does not appear to be a valid FASTQ file:\n"+new String(oct[2]);
//		assert(oct[4]==null || oct[4].startsWith("@")) : "Does not appear to be a valid FASTQ file:\n"+new String(oct[4]);
//		assert(oct[6]==null || oct[6].startsWith("+")) : "Does not appear to be a valid FASTQ file:\n"+new String(oct[6]);
//		
//		if(FORCE_INTERLEAVED){return true;}
//		if(PARSE_CUSTOM && fname.contains("_interleaved.f")){
//			return true;
//		}
//		
//		return testPairNames(oct[0], oct[4]);
//	}
	
	public static boolean isInterleaved(final String fname){
		
		if(!DETECT_QUALITY && !TEST_INTERLEAVED){return FORCE_INTERLEAVED;}
		final String[] oct=getFirstOctet(fname);
		if(oct==null){return FORCE_INTERLEAVED;}
		
		if(DETECT_QUALITY){testQuality(oct);}
		if(TEST_INTERLEAVED){return testInterleaved(oct, fname);}
		return FORCE_INTERLEAVED;
	}
	
	private static String[] getFirstOctet(String fname){
		if(fname==null){return null;}
		if(fname.equalsIgnoreCase("stdin") || fname.toLowerCase().startsWith("stdin.")){return null;}
		
		String[] oct=new String[8];
		int cntr=0;
		
		{
			InputStream is=ReadWrite.getInputStream(fname, false, false);
			BufferedReader br=new BufferedReader(new InputStreamReader(is));
			try {
				for(String s=br.readLine(); s!=null && cntr<8; s=br.readLine()){
					oct[cntr]=s;
					cntr++;
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return oct;
	}
	
	public static byte testQuality(String fname){
		if(fname==null){return ASCII_OFFSET;}
		if(!DETECT_QUALITY || fname.equalsIgnoreCase("stdin") || fname.toLowerCase().startsWith("stdin.")){return ASCII_OFFSET;}
		
		String[] oct=getFirstOctet(fname);
		return testQuality(oct);
	}
	
	private static boolean testInterleaved(String[] oct, String fname){
		if(oct==null || oct.length<8){return false;}
		for(int i=0; i<oct.length; i++){
			if(oct[i]==null){return false;}
		}
		
		assert(oct[0]==null || oct[0].startsWith("@")) : "File "+fname+"\ndoes not appear to be a valid FASTQ file:\n"+new String(oct[0]);
		assert(oct[2]==null || oct[2].startsWith("+")) : "File "+fname+"\ndoes not appear to be a valid FASTQ file:\n"+new String(oct[2]);
		assert(oct[4]==null || oct[4].startsWith("@")) : "File "+fname+"\ndoes not appear to be a valid FASTQ file:\n"+new String(oct[4]);
		assert(oct[6]==null || oct[6].startsWith("+")) : "File "+fname+"\ndoes not appear to be a valid FASTQ file:\n"+new String(oct[6]);
		
		if(FORCE_INTERLEAVED){return true;}
		if(PARSE_CUSTOM && fname.contains("_interleaved.")){return true;}
		
		return testPairNames(oct[0], oct[4]);
	}
	
	private static byte testQuality(String[] oct){
		if(oct==null || oct[0]==null){return ASCII_OFFSET;}
		
		int qflips=0;
		for(int k=0; k<2; k++){
			int a=1+4*k, b=3+4*k;
			if(oct.length<b || oct[a]==null || oct[b]==null){break;}
			byte[] bases=oct[a].getBytes();
			byte[] quals=oct[b].getBytes();
			//		assert(false) : Arrays.toString(quals);
			for(int i=0; i<quals.length; i++){
				quals[i]-=ASCII_OFFSET; //Convert from ASCII33 to native.
				if(DETECT_QUALITY){
					if(ASCII_OFFSET==33 && (quals[i]>QUAL_THRESH /*|| (bases[i]=='N' && quals[i]>20)*/)){
						System.err.println("Changed from ASCII-33 to ASCII-64 on input quality "+(quals[i]+ASCII_OFFSET)+" while prescanning.");
						qflips++;
						ASCII_OFFSET=64;
						if(DETECT_QUALITY_OUT){ASCII_OFFSET_OUT=64;}
						for(int j=0; j<=i; j++){
							quals[j]=(byte)(quals[j]-31);
						}
					}else if(ASCII_OFFSET==64 && (quals[i]<-5)){
						System.err.println("Changed from ASCII-64 to ASCII-33 on input quality "+(quals[i]+ASCII_OFFSET)+" while prescanning.");
						ASCII_OFFSET=33;
						if(DETECT_QUALITY_OUT){ASCII_OFFSET_OUT=33;}
						qflips++;
						for(int j=0; j<=i; j++){
							quals[j]=(byte)(quals[j]+31);
						}
					}
				}
				assert(quals[i]>=-5) : "ASCII encoding for quality (currently ASCII-"+ASCII_OFFSET+") appears to be wrong.\n"
					+oct[k]+"\n"+oct[k+3]+"\n"+Arrays.toString(oct[k+3].getBytes());
				assert(qflips<2) : "Failed to auto-detect quality coding; quitting.";
			}
		}
		
		return ASCII_OFFSET;
	}
	
	public static boolean testPairNames(Read r1, Read r2){
		if(r1==null || r2==null){return false;}
		return testPairNames(r1.id, r2.id);
	}
	
	public static boolean testPairNames(String id1, String id2){
		
		if(id1==null || id2==null){return false;}
		
		final int idxSlash1=id1.lastIndexOf('/');
		final int idxSlash2=id2.lastIndexOf('/');
		final int idxSpace1=id1.indexOf(' ');
		final int idxSpace2=id2.indexOf(' ');
		//			System.out.println("idxSlash1="+idxSlash1+", idxSlash2="+idxSlash2+", idxSpace1="+idxSpace1+", idxSpace2="+idxSpace2);
		if(idxSlash1==idxSlash2 && idxSlash1>1){
			//				System.out.println("A");
			String[] split1=id1.split("/");
			String[] split2=id2.split("/");
			//				System.out.println(Arrays.toString(split1));
			//				System.out.println(Arrays.toString(split2));

			if(split1.length>1 && split2.length>1 && split1[0].equals(split2[0])){
				//					System.out.println("B");
				if(split1[split1.length-1].contains(" ")){
					split1[split1.length-1]=split1[split1.length-1].split(" ")[0];
					//						System.out.println("B1: "+Arrays.toString(split1));
				}
				if(split2[split2.length-1].contains(" ")){
					split2[split2.length-1]=split2[split2.length-1].split(" ")[0];
					//						System.out.println("B2: "+Arrays.toString(split2));
				}
				if(split1[split1.length-1].equals("1") && split2[split2.length-1].equals("2")){
					//						System.out.println("B3");
					return true;
				}
			}
		}

		if(idxSpace1==idxSpace2 && idxSpace1>=0){
			//				System.out.println("C");
			if(idxSpace1==idxSpace2 && idxSpace1>1){
				//					System.out.println("D");
				String[] split1=id1.split(" ");
				String[] split2=id2.split(" ");
				//					System.out.println(Arrays.toString(split1));
				//					System.out.println(Arrays.toString(split2));

				if(split1.length>1 && split2.length>1 && split1[0].equals(split2[0])){
					//						System.out.println("E");
					if(split1[1].startsWith("1:") && split2[1].startsWith("2:")){return true;}
				}
			}
		}
		return false;
	}
	
	public static String[] toFASTQ(Read r){
		String id=customID(r);
		return toFASTQ(r.bases, r.quality, id==null ? ""+r.numericID : id);
	}
	
	public static String customID(Read r){
		if(PARSE_CUSTOM && (r.chrom>-1 && r.stop>-1)){
			if(Data.GENOME_BUILD>=0){
				final int chrom1=r.chrom;
				final int start1=r.start;
				final int stop1=r.stop;
				int idx1=Data.scaffoldIndex(chrom1, (start1+stop1)/2);
				byte[] name1=Data.scaffoldNames[chrom1][idx1];
				int a1=Data.scaffoldRelativeLoc(chrom1, start1, idx1);
				if(r.mate==null || !ADD_PAIRNUM_TO_CUSTOM_ID){
					return (r.id==null ? ""+r.numericID : r.id)+"_chr"+r.chrom+"_"+r.strand()+"_"+r.start+"_"+r.stop+"_"+a1+"_"+new String(name1);
				}else{
					return (r.id==null ? ""+r.numericID : r.id)+"_chr"+r.chrom+"_"+r.strand()+"_"+r.start+"_"+r.stop+"_"+a1+"_"+new String(name1)+" /"+(r.pairnum()+1);
				}
			}else{
				if(r.mate==null || !ADD_PAIRNUM_TO_CUSTOM_ID){
					return (r.id==null ? ""+r.numericID : r.id)+"_chr"+r.chrom+"_"+r.strand()+"_"+r.start+"_"+r.stop;
				}else{
					return (r.id==null ? ""+r.numericID : r.id)+"_chr"+r.chrom+"_"+r.strand()+"_"+r.start+"_"+r.stop+" /"+(r.pairnum()+1);
				}
			}
		}
		return r.id;
	}
	
	private static int fastqLength(Read r){
		int len=6; //newlines, @, +
		len+=(r.id==null ? Tools.stringLength(r.numericID) : r.id.length());
		len+=(r.bases==null ? 0 : r.bases.length);
		len+=(r.quality==null ? 0 : r.quality.length);
		return len;
	}
	
	public static ByteBuilder toFASTQ(Read r, ByteBuilder bb){
		int len=fastqLength(r);
		final String id;
		final byte[] bases=r.bases, quals=r.quality;
		if(PARSE_CUSTOM && (r.chrom>-1 && r.stop>-1)){
			id=customID(r);
			if(id!=null){len+=id.length();}
		}else{
			id=r.id;
		}
		if(bb==null){bb=new ByteBuilder(len);}
		else{bb.ensureExtra(len);}
		
		bb.append('@');
		if(id==null){bb.append(r.numericID);}
		else{bb.append(id);}
		bb.append('\n');
		
//		if(bases!=null){for(byte b : bases){sb.append((char)b);}}
//		sb.append('\n');
//		sb.append('+');
//		sb.append('\n');
//		if(quals!=null){for(byte b : quals){sb.append((char)(b+ASCII_OFFSET_OUT));}}
		
		if(bases==null){
			bb.append('\n').append('+').append('\n');
			if(verbose){System.err.println("A:\n"+bb);}
		}else{
			bb.append(bases);
			bb.append('\n').append('+').append('\n');
			if(verbose){System.err.println("B:\n"+bb);}
			if(quals==null){
				final byte q=(byte)(30+ASCII_OFFSET_OUT);
				final int blen=bases.length;
				bb.ensureExtra(blen);
				for(int i=0, j=bb.length; i<blen; i++, j++){bb.array[j]=q;}
				bb.length+=blen;
				if(verbose){System.err.println("C:\n"+bb);}
			}else{
				bb.ensureExtra(quals.length);
				for(int i=0, j=bb.length; i<quals.length; i++, j++){
					byte q=quals[i];
					bb.array[j]=(byte)(q+ASCII_OFFSET_OUT);
				}
				bb.length+=quals.length;
				if(verbose){System.err.println("D:\n"+bb);}
			}
		}
		if(verbose){System.err.println("E:\n"+bb);}
		
//		sb.append('\n');
		return bb;
	}
	
	public static StringBuilder toFASTQ(Read r, StringBuilder sb){
		int len=fastqLength(r);
		final String id;
		final byte[] bases=r.bases, quals=r.quality;
		if(PARSE_CUSTOM && (r.chrom>-1 && r.stop>-1)){
			id=customID(r);
			if(id!=null){len+=id.length();}
		}else{
			id=r.id;
		}
		if(sb==null){sb=new StringBuilder(len);}
		else{sb.ensureCapacity(len);}
		
		sb.append('@');
		if(id==null){sb.append(r.numericID);}
		else{sb.append(id);}
		sb.append('\n');
		
//		if(bases!=null){for(byte b : bases){sb.append((char)b);}}
//		sb.append('\n');
//		sb.append('+');
//		sb.append('\n');
//		if(quals!=null){for(byte b : quals){sb.append((char)(b+ASCII_OFFSET_OUT));}}
		
		if(bases==null){
			sb.append('\n').append('+').append('\n');
		}else{
			char[] buffer=Shared.getTLCB(bases.length);
			for(int i=0; i<bases.length; i++){buffer[i]=(char)bases[i];}
			sb.append(buffer, 0, bases.length);
			sb.append('\n').append('+').append('\n');
			if(quals==null){
				final char q=(char)(30+ASCII_OFFSET_OUT);
				final int blen=bases.length;
				for(int i=0; i<blen; i++){buffer[i]=q;}
				sb.append(buffer, 0, blen);
			}else{
				for(int i=0; i<quals.length; i++){buffer[i]=(char)(quals[i]+ASCII_OFFSET_OUT);}
				sb.append(buffer, 0, quals.length);
			}
		}
		
//		sb.append('\n');
		return sb;
	}
	
	public static String[] toFASTQ(byte[] bases, byte[] quality, String identifier){
		String[] out=new String[4];
		
		identifier=(identifier==null ? ""+incr() : identifier);
		if(quality==null){
			byte[] x=new byte[bases.length];
			for(int i=0; i<bases.length; i++){
				x[i]=30;
			}
			quality=x;
		}
		
		byte[] q2=new byte[quality.length];
		for(int i=0; i<q2.length; i++){q2[i]=(byte)(quality[i]+ASCII_OFFSET_OUT);}
		
		assert(quality.length==bases.length);
		
		out[0]="@"+identifier;
		out[1]=new String(bases);
		out[2]="+"/*+identifier*/;
		out[3]=new String(q2);
		
		return out;
	}
	
	
	public static Read[] toReads(TextFile tf, int maxReadsToReturn, boolean colorspace, long numericID, boolean interleaved){
		ArrayList<Read> list=toReadList(tf, maxReadsToReturn, colorspace, numericID, interleaved);
		assert(list.size()<=maxReadsToReturn);
		return list.toArray(new Read[list.size()]);
	}
	
	public static final String makeId(String s){
		if(s==null || s.length()<1){return null;}
		char c=s.charAt(0);
		int start=0, stop=s.length();
		if(c=='@' || c=='>'){start=1;}
		if(Shared.TRIM_READ_COMMENTS){
			for(int i=start; i<stop; i++){
				if(Character.isWhitespace(s.charAt(i))){
					stop=i;
					break;
				}
			}
		}
		return stop<=start ? null : start==0 && stop==s.length() ? s : s.substring(start, stop);
	}
	
	public static final String makeId(byte[] s){
		if(s==null || s.length<1){return null;}
		byte c=s[0];
		int start=0, stop=s.length;
		if(c=='@' || c=='>'){start=1;}
		if(Shared.TRIM_READ_COMMENTS){
			for(int i=start; i<stop; i++){
				if(Character.isWhitespace(s[i])){
					stop=i;
					break;
				}
			}
		}
		return stop<=start ? null : new String(s, start, stop-start);
	}
	
	public static ArrayList<Read> toReadList(TextFile tf, int maxReadsToReturn, boolean colorspace, long numericID, boolean interleaved){
		String s=null;
		ArrayList<Read> list=new ArrayList<Read>(Data.min(16384, maxReadsToReturn));
		
		String[] quad=new String[4];
		
		int cntr=0;
		int added=0;
		
		Read prev=null;
		
		for(s=tf.nextLine(); s!=null && added<maxReadsToReturn; s=tf.nextLine()){
			quad[cntr]=s;
			cntr++;
			if(cntr==4){
				assert(quad[0].startsWith("@")) : "\nError in "+tf.name+", line "+tf.lineNum+"\n"+quad[0]+"\n"+quad[1]+"\n"+quad[2]+"\n"+quad[3]+"\n";
				assert(quad[2].startsWith("+")) : "\nError in "+tf.name+", line "+tf.lineNum+"\n"+quad[0]+"\n"+quad[1]+"\n"+quad[2]+"\n"+quad[3]+"\n";
				
//				if(quad[0].startsWith("@HW") || quad[0].startsWith("@FC")){ascii_offset=66;} //TODO: clumsy
				
				final String id=makeId(quad[0]);
				
				Read r=null;
				
				byte[] bases=quad[1].getBytes();
				byte[] quals=quad[3].getBytes();
//				assert(false) : Arrays.toString(quals);
				for(int i=0; i<quals.length; i++){
					quals[i]-=ASCII_OFFSET; //Convert from ASCII33 to native.
					if(DETECT_QUALITY && ASCII_OFFSET==33 && (quals[i]>QUAL_THRESH /*|| (bases[i]=='N' && quals[i]>20)*/)){
						if(numericID<1){
							System.err.println("Changed from ASCII-33 to ASCII-64 on input "+((char)quals[i])+": "+quals[i]+" -> "+(quals[i]-31));
						}else{
							System.err.println("Warning! Changed from ASCII-33 to ASCII-64 on input "+((char)quals[i])+": "+quals[i]+" -> "+(quals[i]-31));
							System.err.println("Up to "+numericID+" prior reads may have been generated with incorrect qualities.");
							System.err.println("If this is a problem you may wish to re-run with the flag 'qin=64'.");
						}
						ASCII_OFFSET=64;
						for(int j=0; j<=i; j++){
							quals[j]=(byte)(quals[j]-31);
						}
					}
					assert(quals[i]>=-5) : "\n"+quad[0]+"\n"+quad[3];
				}
//				assert(false) : Arrays.toString(quals);
//				assert(false) : new String(quad[0]);
				if(PARSE_CUSTOM && quad[0]!=null && quad[0].indexOf('_')>0){
					String[] answer=quad[0].split("_");
					if(answer.length>=5){
						try {
							byte trueChrom=Gene.toChromosome(answer[1]);
							byte trueStrand=Byte.parseByte(answer[2]);
							int trueLoc=Integer.parseInt(answer[3]);
							int trueStop=Integer.parseInt(answer[4]);
							r=new Read(bases, trueChrom, trueStrand, trueLoc, trueStop, id, quals, colorspace, numericID);
							r.setSynthetic(true);
						} catch (NumberFormatException e) {}
					}
				}
				if(r==null){
					r=new Read(bases, 0, (byte)0, 0, 0, id, quals, colorspace, numericID);
				}
				
				cntr=0;
				
				if(interleaved){
					if(prev==null){prev=r;}
					else{
						prev.mate=r;
						r.mate=prev;
						r.setPairnum(1);
						list.add(prev);
						added++;
						numericID++;
						prev=null;
					}
				}else{
					list.add(r);
					added++;
					numericID++;
				}
				
				
				if(added>=maxReadsToReturn){break;}
				
//				System.out.println(r.chrom+", "+r.strand+", "+r.loc);
//				assert(false);
			}
		}
		assert(list.size()<=maxReadsToReturn);
		return list;
	}
	
	public static ArrayList<Read> toReadList(ByteFile tf, int maxReadsToReturn, boolean colorspace, long numericID, boolean interleaved){
		byte[] s=null;
		ArrayList<Read> list=new ArrayList<Read>(Data.min(16384, maxReadsToReturn));
		
		byte[][] quad=new byte[4][];
		
		int cntr=0;
		int added=0;
		
		Read prev=null;
		
		for(s=tf.nextLine(); s!=null && added<maxReadsToReturn; s=tf.nextLine()){
			quad[cntr]=s;
			cntr++;
			if(cntr==4){
//				
//				if(verbose){
//					System.err.println("\nASCII offset is "+ASCII_OFFSET);
//					System.err.println("quad:");
//					System.err.println(new String(quad[0]));
//					System.err.println(new String(quad[1]));
//					System.err.println(new String(quad[2]));
//					System.err.println(new String(quad[3]));
//				}
//				
//				assert(quad[0][0]==(byte)'@') : "\nError in "+tf.name+", line "+tf.lineNum()+"\n"+
//					new String(quad[0])+"\n"+new String(quad[1])+"\n"+new String(quad[2])+"\n"+new String(quad[3])+"\n";
//				assert(quad[2][0]==(byte)'+') : "\nError in "+tf.name+", line "+tf.lineNum()+"\n"+
//					new String(quad[0])+"\n"+new String(quad[1])+"\n"+new String(quad[2])+"\n"+new String(quad[3])+"\n";
//
////				if(quad[0].startsWith("@HW") || quad[0].startsWith("@FC")){ascii_offset=66;} //TODO: clumsy
//				
//				final String id=makeId(quad[0]);
//				
//				Read r=null;
//				
//				byte[] bases=quad[1];
//				byte[] quals=quad[3];
////				assert(false) : Arrays.toString(quals);
//				for(int i=0; i<quals.length; i++){
//					quals[i]-=ASCII_OFFSET; //Convert from ASCII33 to native.
//					if(DETECT_QUALITY && ASCII_OFFSET==33 && (quals[i]>QUAL_THRESH /*|| (bases[i]=='N' && quals[i]>20)*/)){
//						if(numericID<1){
//							System.err.println("Changed from ASCII-33 to ASCII-64 on input "+((char)quals[i])+": "+quals[i]+" -> "+(quals[i]-31));
//						}else{
//							System.err.println("Warning! Changed from ASCII-33 to ASCII-64 on input "+((char)quals[i])+": "+quals[i]+" -> "+(quals[i]-31));
//							System.err.println("Up to "+numericID+" prior reads may have been generated with incorrect qualities.");
//							System.err.println("If this is a problem you may wish to re-run with the flag 'qin=64'.");
//							errorState=true;
//						}
//						ASCII_OFFSET=64;
//						for(int j=0; j<=i; j++){
//							quals[j]=(byte)(quals[j]-31);
//						}
//					}
//					if(quals[i]<-5){
//						if(!negativeFive){
//							for(int j=0; j<quals.length; j++){quals[j]=Tools.max(quals[j], (byte)33);}
//							System.err.println("\nThe ASCII quality encoding offset ("+ASCII_OFFSET+") is not set correctly; quality value below -5.\n" +
//									"Please re-run with the flag 'qin=33'.\nProblematic read number "+numericID+":\n" +
//						
//							"\n"+new String(quad[0])+"\n"+new String(quad[1])+"\n"+new String(quad[2])+"\n"+new String(quad[3])+"\n");
//						}
//						errorState=true;
//						negativeFive=true;
//						return null;
//					}
//					assert(quals[i]>=-5);
////					assert(quals[i]>=-5) : "The ASCII quality encoding level is not set correctly.  Quality value below -5:" +
////							"\n"+new String(quad[0])+"\n"+new String(quad[1])+"\n"+new String(quad[2])+"\n"+new String(quad[3]);
//				}
////				assert(false) : Arrays.toString(quals);
////				assert(false) : PARSE_CUSTOM+"\n"+new String(quad[0]);
//				if(PARSE_CUSTOM){
//					if(quad[0]!=null && Tools.indexOf(quad[0], (byte)'_')>0){
//						String temp=new String(quad[0]);
//						if(temp.endsWith(" /1") || temp.endsWith(" /2")){temp=temp.substring(0, temp.length()-3);}
//						String[] answer=temp.split("_");
//
//						if(answer.length>=5){
//							try {
//								byte trueChrom=Gene.toChromosome(answer[1]);
//								byte trueStrand=Byte.parseByte(answer[2]);
//								int trueLoc=Integer.parseInt(answer[3]);
//								int trueStop=Integer.parseInt(answer[4]);
//								r=new Read(bases, trueChrom, trueStrand, trueLoc, trueStop, id, quals, colorspace, numericID);
//								r.setSynthetic(true);
//							} catch (NumberFormatException e) {
//								PARSE_CUSTOM=false;
//								System.err.println("Turned off PARSE_CUSTOM because could not parse "+new String(quad[0]));
//							}
//						}else{
//							PARSE_CUSTOM=false;
//							System.err.println("Turned off PARSE_CUSTOM because answer="+Arrays.toString(answer));
//						}
//					}else{
//						PARSE_CUSTOM=false;
//						System.err.println("Turned off PARSE_CUSTOM because quad[0]="+new String(quad[0])+", index="+Tools.indexOf(quad[0], (byte)'_'));
//					}
//				}
//				if(r==null){
//					r=new Read(bases, 0, (byte)0, 0, 0, id, quals, colorspace, numericID);
//				}
				
				
				Read r=quadToRead(quad, true, false, colorspace, tf, numericID);
				cntr=0;
				
				if(interleaved){
					if(prev==null){prev=r;}
					else{
						prev.mate=r;
						r.mate=prev;
						r.setPairnum(1);
						list.add(prev);
						added++;
						numericID++;
						prev=null;
					}
				}else{
					list.add(r);
					added++;
					numericID++;
				}
				
				
				if(added>=maxReadsToReturn){break;}
				
//				System.out.println(r.chrom+", "+r.strand+", "+r.loc);
//				assert(false);
			}
		}
		assert(list.size()<=maxReadsToReturn);
		return list;
	}
	
	public static byte[][] scarfToQuad(final byte[] scarf, byte[][] quad){
		
		int a=-1, b=-1;
		final byte colon=':';
		for(int i=scarf.length-1; i>=0; i--){
			if(scarf[i]==colon){
				if(b<0){b=i;}
				else{
					assert(a<0);
					a=i;
					break;
				}
			}
		}
		if(a<0 || b<0){
			throw new RuntimeException("Misformatted scarf line: "+new String(scarf));
		}
		if(quad==null){quad=new byte[4][];}
		quad[0]=Arrays.copyOfRange(scarf, 0, a);
		quad[1]=Arrays.copyOfRange(scarf, a+1, b);
		quad[3]=Arrays.copyOfRange(scarf, b+1, scarf.length);
		return quad;
	}
	
	public static Read quadToRead(final byte[][] quad, boolean fastq, boolean scarf, boolean colorspace, ByteFile tf, long numericID){

		if(verbose){
			System.err.println("\nASCII offset is "+ASCII_OFFSET);
			System.err.println("quad:");
			System.err.println(new String(quad[0]));
			System.err.println(new String(quad[1]));
			System.err.println(new String(quad[2]));
			System.err.println(new String(quad[3]));
		}

		assert(scarf || quad[0][0]==(byte)'@') : "\nError in "+tf.name()+", line "+tf.lineNum()+"\n"+
		new String(quad[0])+"\n"+new String(quad[1])+"\n"+new String(quad[2])+"\n"+new String(quad[3])+"\n";
		assert(scarf || quad[2][0]==(byte)'+') : "\nError in "+tf.name()+", line "+tf.lineNum()+"\n"+
		new String(quad[0])+"\n"+new String(quad[1])+"\n"+new String(quad[2])+"\n"+new String(quad[3])+"\n";

		//			if(quad[0].startsWith("@HW") || quad[0].startsWith("@FC")){ascii_offset=66;} //TODO: clumsy

		final String id=makeId(quad[0]);

		Read r=null;

		byte[] bases=quad[1];
		byte[] quals=quad[3];
		//			assert(false) : Arrays.toString(quals);
		for(int i=0; i<quals.length; i++){
			quals[i]-=ASCII_OFFSET; //Convert from ASCII33 to native.
			if(DETECT_QUALITY && ASCII_OFFSET==33 && (quals[i]>QUAL_THRESH /*|| (bases[i]=='N' && quals[i]>20)*/)){
				if(numericID<1){
					System.err.println("Changed from ASCII-33 to ASCII-64 on input "+((char)quals[i])+": "+quals[i]+" -> "+(quals[i]-31));
				}else{
					System.err.println("Warning! Changed from ASCII-33 to ASCII-64 on input "+((char)quals[i])+": "+quals[i]+" -> "+(quals[i]-31));
					System.err.println("Up to "+numericID+" prior reads may have been generated with incorrect qualities.");
					System.err.println("If this is a problem you may wish to re-run with the flag 'qin=64'.");
					errorState=true;
				}
				ASCII_OFFSET=64;
				for(int j=0; j<=i; j++){
					quals[j]=(byte)(quals[j]-31);
				}
			}
			if(quals[i]<-5){
				if(!negativeFive){
					for(int j=0; j<quals.length; j++){quals[j]=Tools.max(quals[j], (byte)33);}
					System.err.println("\nThe ASCII quality encoding offset ("+ASCII_OFFSET+") is not set correctly; quality value below -5.\n" +
							"Please re-run with the flag 'qin=33'.\nProblematic read number "+numericID+":\n" +

						"\n"+new String(quad[0])+"\n"+new String(quad[1])+"\n"+new String(quad[2])+"\n"+new String(quad[3])+"\n");
				}
				errorState=true;
				negativeFive=true;
				return null;
			}
			assert(quals[i]>=-5);
			//				assert(quals[i]>=-5) : "The ASCII quality encoding level is not set correctly.  Quality value below -5:" +
			//						"\n"+new String(quad[0])+"\n"+new String(quad[1])+"\n"+new String(quad[2])+"\n"+new String(quad[3]);
		}
		//			assert(false) : Arrays.toString(quals);
		//			assert(false) : PARSE_CUSTOM+"\n"+new String(quad[0]);
		if(PARSE_CUSTOM){
			if(quad[0]!=null && Tools.indexOf(quad[0], (byte)'_')>0){
				String temp=new String(quad[0]);
				if(temp.endsWith(" /1") || temp.endsWith(" /2")){temp=temp.substring(0, temp.length()-3);}
				String[] answer=temp.split("_");

				if(answer.length>=5){
					try {
						byte trueChrom=Gene.toChromosome(answer[1]);
						byte trueStrand=Byte.parseByte(answer[2]);
						int trueLoc=Integer.parseInt(answer[3]);
						int trueStop=Integer.parseInt(answer[4]);
						r=new Read(bases, trueChrom, trueStrand, trueLoc, trueStop, id, quals, colorspace, numericID);
						r.setSynthetic(true);
					} catch (NumberFormatException e) {
						PARSE_CUSTOM=false;
						System.err.println("Turned off PARSE_CUSTOM because could not parse "+new String(quad[0]));
					}
				}else{
					PARSE_CUSTOM=false;
					System.err.println("Turned off PARSE_CUSTOM because answer="+Arrays.toString(answer));
				}
			}else{
				PARSE_CUSTOM=false;
				System.err.println("Turned off PARSE_CUSTOM because quad[0]="+new String(quad[0])+", index="+Tools.indexOf(quad[0], (byte)'_'));
			}
		}
		if(r==null){
			r=new Read(bases, 0, (byte)0, 0, 0, id, quals, colorspace, numericID);
		}
		return r;
	}
	
	public static ArrayList<Read> toScarfReadList(ByteFile tf, int maxReadsToReturn, boolean colorspace, long numericID, boolean interleaved){
		byte[] s=null;
		ArrayList<Read> list=new ArrayList<Read>(Data.min(16384, maxReadsToReturn));
		
		byte[][] quad=new byte[4][];
		
		int added=0;
		
		Read prev=null;
		
		for(s=tf.nextLine(); s!=null && added<maxReadsToReturn; s=tf.nextLine()){
			scarfToQuad(s, quad);
			Read r=quadToRead(quad, false, true, colorspace, tf, numericID);

			if(interleaved){
				if(prev==null){prev=r;}
				else{
					prev.mate=r;
					r.mate=prev;
					r.setPairnum(1);
					list.add(prev);
					added++;
					numericID++;
					prev=null;
				}
			}else{
				list.add(r);
				added++;
				numericID++;
			}


			if(added>=maxReadsToReturn){break;}
		}
		assert(list.size()<=maxReadsToReturn);
		return list;
	}
	
	public static String qualToString(byte[] quals){
		byte[] q2=new byte[quals.length];
		for(int i=0; i<quals.length; i++){
			q2[i]=(byte)(quals[i]+ASCII_OFFSET);
		}
		return new String(q2);
	}
	
	/** Return true if this has detected an error */
	public static boolean errorState(){return errorState;}
	/** TODO */
	private static boolean errorState=false;
	private static boolean negativeFive=false;
	
	private static synchronized long incr(){return incr++;}
	private static long incr=10000000000L;
	
	public static boolean PARSE_CUSTOM=false;
	public static byte ASCII_OFFSET=33;
	public static byte ASCII_OFFSET_OUT=33;
	public static boolean TEST_INTERLEAVED=true;
	public static boolean FORCE_INTERLEAVED=false;
	public static boolean DETECT_QUALITY=true;
	public static boolean DETECT_QUALITY_OUT=true;
	public static boolean ADD_PAIRNUM_TO_CUSTOM_ID=true;
	
	public static final int QUAL_THRESH=54;
	public static boolean verbose=false;
	
}
