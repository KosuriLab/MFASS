package align2;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Random;

import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;



import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.FastaToChromArrays;
import dna.Gene;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.SummaryFile;
import fileIO.TextStreamWriter;

public final class RandomReads3 {
	
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		t.start();

		FASTQ.ADD_PAIRNUM_TO_CUSTOM_ID=false;
		FastaReadInputStream.SPLIT_READS=false;
		FastaReadInputStream.MIN_READ_LEN=1;
		Data.GENOME_BUILD=-1;
		int build=-1;
		String ref=null;
		String out=null;
		
		int number=0;
		int minlen=100;
		int maxlen=100;

		int minInsLen=1;
		int minSubLen=2;
		int minDelLen=1;
		int minNLen=1;

		int maxInsLen=12;
		int maxSubLen=12;
		int maxDelLen=400;
		int maxNLen=1;
		
		int minChrom=-1;
		int maxChrom=-1;
		
		int maxSnps=3;
		int maxInss=2;
		int maxDels=2;
		int maxSubs=2;
		int maxNs=0;

		float snpRate=0;
		float insRate=0;
		float delRate=0;
		float subRate=0;
		float nRate=0;
		PERFECT_READ_RATIO=0;

//		float snpRate=0.4f;
//		float insRate=0.2f;
//		float delRate=0.2f;
//		float subRate=0.2f;
//		float nRate=0.2f;
//		PERFECT_READ_RATIO=0.5f;
		
		String adapter=null;
		
		long seed2=Long.MIN_VALUE;
		
		int minQuality=28;
		int midQuality=32;
		int maxQuality=36;
		
		int minInsert=-1, maxInsert=-1;
		
		boolean paired=false;
		boolean colorspace=false;
		
		ReadWrite.USE_PIGZ=true;
		
		for(int i=0; i<args.length; i++){
//			assert(s.contains("=")) : "All arguments must be of the form word=number, e.g., reads=10000";
			final String arg=args[i];
			final String[] split=arg.split("=");
			assert(split.length<=2);
//			assert(split.length==2);
			final String a=split[0].toLowerCase();
			final String b=(split.length<2 ? "true" : split[1]);
			
			int x=-1;
			try {
				x=Integer.parseInt(b);
			} catch (NumberFormatException e) {}

			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("reads")){
				number=x;
			}else if(a.equals("len") || a.equals("length") || a.equals("readlen")){
				minlen=maxlen=x;
			}else if(a.equals("overwrite") || a.equals("ow")){
				OVERWRITE=Tools.parseBoolean(b);
			}else if(a.startsWith("minlen")){
				minlen=x;
				maxlen=Tools.max(minlen, maxlen);
			}else if(a.startsWith("maxlen")){
				maxlen=x;
				minlen=Tools.min(minlen, maxlen);
			}else if(a.equals("adapter")){
				adapter=b;
			}else if(a.equals("amp")){
				AMP=Integer.parseInt(b);
			}else if(a.equals("slashes") || a.equals("addslashes") || a.equals("slash") || a.equals("addpairnum") || a.equals("pairnum")){
				FASTQ.ADD_PAIRNUM_TO_CUSTOM_ID=Tools.parseBoolean(b);
			}else if(a.equals("snprate")){
				snpRate=Float.parseFloat(b);
			}else if(a.equals("subrate")){
				subRate=Float.parseFloat(b);
			}else if(a.equals("delrate")){
				delRate=Float.parseFloat(b);
			}else if(a.equals("insrate")){
				insRate=Float.parseFloat(b);
			}else if(a.equals("nrate")){
				nRate=Float.parseFloat(b);
			}else if(a.equals("maxsnps")){
				maxSnps=Integer.parseInt(b);
			}else if(a.equals("maxdels")){
				maxDels=Integer.parseInt(b);
			}else if(a.equals("maxsubs")){
				maxSubs=Integer.parseInt(b);
			}else if(a.equals("maxinss")){
				maxInss=Integer.parseInt(b);
			}else if(a.equals("maxns")){
				maxNs=Integer.parseInt(b);
			}else if(a.startsWith("maxdellen")){
				maxDelLen=Integer.parseInt(b);
			}else if(a.startsWith("maxsublen")){
				maxSubLen=Integer.parseInt(b);
			}else if(a.startsWith("maxinslen")){
				maxInsLen=Integer.parseInt(b);
			}else if(a.startsWith("maxnlen")){
				maxNLen=Integer.parseInt(b);
			}else if(a.startsWith("mindellen")){
				minDelLen=Integer.parseInt(b);
			}else if(a.startsWith("minsublen")){
				minSubLen=Integer.parseInt(b);
			}else if(a.startsWith("mininslen")){
				minInsLen=Integer.parseInt(b);
			}else if(a.startsWith("minnlen")){
				minNLen=Integer.parseInt(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
			}else if(a.startsWith("seed")){
				seed2=Long.parseLong(b);
			}else if(a.equals("ref") || a.equals("reference")){
				ref=b;
			}else if(a.equals("s") || a.startsWith("snp")){
				maxSnps=x;
				snpRate=1;
			}else if(a.equals("i") || a.startsWith("ins")){
				maxInss=(x>0 ? 1 : 0);
				maxInsLen=x;
				insRate=1;
			}else if(a.equals("d") || a.startsWith("del")){
				maxDels=(x>0 ? 1 : 0);
				maxDelLen=x;
				delRate=1;
			}else if(a.equals("u") || a.startsWith("sub")){
				maxSubs=(x>0 ? 1 : 0);
				maxSubLen=x;
				subRate=1;
			}else if(a.equals("n")){
				maxNs=x;
				nRate=1;
				minNLen=maxNLen=1;
			}else if(a.startsWith("minchrom")){
				minChrom=x;
			}else if(a.equals("int") || a.equals("interleaved") || a.equals("interleave")){
				OUTPUT_INTERLEAVED=Tools.parseBoolean(b);
				if(OUTPUT_INTERLEAVED){paired=true;}
			}else if(a.equals("biasedsnps")){
				BIASED_SNPS=Tools.parseBoolean(b);
			}else if(a.startsWith("maxchrom")){
				maxChrom=x;
			}else if(a.startsWith("build") || a.startsWith("genome")){
				build=x;
//				assert(false) : "Set genome to "+x;
			}else if(a.startsWith("minq")){
				minQuality=x;
				midQuality=Tools.max(midQuality,  minQuality);
				maxQuality=Tools.max(maxQuality,  minQuality);
			}else if(a.startsWith("midq")){
				midQuality=x;
			}else if(a.startsWith("maxq")){
				maxQuality=x;
				midQuality=Tools.min(midQuality,  maxQuality);
				minQuality=Tools.min(minQuality,  maxQuality);
			}else if(a.startsWith("qual") || a.equals("q")){
				minQuality=midQuality=maxQuality=x;
			}else if(a.equals("mininsert")){
				minInsert=x;
			}else if(a.equals("maxinsert")){
				maxInsert=x;
			}else if(a.startsWith("minmid")){
				mateMiddleMin=x;
			}else if(a.startsWith("maxmid")){
				mateMiddleMax=x;
			}else if(a.startsWith("paired")){
				paired=Tools.parseBoolean(b);
			}else if(a.startsWith("superflat")){
				SUPERFLAT_DIST=Tools.parseBoolean(b);
			}else if(a.startsWith("flat")){
				FLAT_DIST=Tools.parseBoolean(b);
			}else if(a.startsWith("bell") || a.startsWith("gauss") || a.startsWith("round")){
				BELL_DIST=Tools.parseBoolean(b);
			}else if(a.startsWith("unique")){
				USE_UNIQUE_SNPS=Tools.parseBoolean(b);
			}else if(a.startsWith("adderrors") || a.startsWith("usequality")){
				ADD_ERRORS_FROM_QUALITY=Tools.parseBoolean(b);
			}else if(a.startsWith("replacenoref")){
				REPLACE_NOREF=Tools.parseBoolean(b);
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("colorspace")){
				colorspace=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("ziplevel") || a.equals("zl")){
				if(x>=0){
					ReadWrite.ZIPLEVEL=Tools.min(x, 9);
				}
			}else if(a.equals("ext") || a.equals("extension")){
				fileExt=b;
				if(fileExt==null){fileExt=".fq.gz";}
				if(!fileExt.startsWith(".")){fileExt="."+fileExt;}
			}else if(a.equals("perfect")){
				PERFECT_READ_RATIO=Float.parseFloat(b);
			}else if(a.equals("singlescaffold")){
				FORCE_SINGLE_SCAFFOLD=Tools.parseBoolean(b);
			}else if(a.equals("minoverlap") || a.equals("overlap")){
				MIN_SCAFFOLD_OVERLAP=Integer.parseInt(b);
			}else if(a.equals("usegzip") || a.equals("gzip")){
				ReadWrite.USE_GZIP=Tools.parseBoolean(b);
			}else if(a.equals("usepigz") || a.equals("pigz")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					int zt=Integer.parseInt(b);
					if(zt<1){ReadWrite.USE_PIGZ=false;}
					else{
						ReadWrite.USE_PIGZ=true;
						if(zt>1){
							ReadWrite.MAX_ZIP_THREADS=zt;
							ReadWrite.ZIP_THREAD_DIVISOR=1;
						}
					}
				}else{ReadWrite.USE_PIGZ=Tools.parseBoolean(b);}
			}else if(a.equals("usegunzip") || a.equals("gunzip")){
				ReadWrite.USE_GUNZIP=Tools.parseBoolean(b);
			}else if(a.equals("useunpigz") || a.equals("unpigz")){
				ReadWrite.USE_UNPIGZ=Tools.parseBoolean(b);
			}else{throw new RuntimeException("Unknown parameter "+args[i]);}

		}
//		assert(false) : OUTPUT_INTERLEAVED;
		assert(build>=0) : "Please specify a genome.";

		if(minInsert>-1){mateMiddleMin=minInsert-2*maxlen;}
		if(maxInsert>-1){mateMiddleMax=maxInsert-2*minlen;}
		
		ArrayList<ChromosomeArray> chromlist=null;
		if(ref!=null){
			chromlist=writeRef(ref, build);
		}
		
		Data.setGenome(build);
		if(minChrom<1){minChrom=1;}
		if(maxChrom<1){maxChrom=Data.numChroms;}
		
		if(chromlist==null){
			Data.loadChromosomes(minChrom, maxChrom);
		}else{
			assert(chromlist.size()==maxChrom-minChrom+1) : chromlist.size()+", "+minChrom+", "+maxChrom;
			for(ChromosomeArray cha : chromlist){
				Data.chromosomePlusMatrix[cha.chromosome]=cha;
			}
		}
		if(Shared.TRIM_READ_COMMENTS){Data.trimScaffoldNames();}
		
		if(number<1){
			Data.sysout.println("No reads to generate; quitting.");
			return;
		}
		
		RandomReads3 rr=(seed2==Long.MIN_VALUE ? new RandomReads3(paired) : 
			new RandomReads3((seed2==-1 ? System.nanoTime() : seed2), paired));
		if(adapter!=null){
			rr.adapter1=adapter.getBytes();
			rr.adapter2=AminoAcid.reverseComplementBases(rr.adapter1);
			rr.adapter1=rr.adapter2; //For PacBio, since adapters never appear in plus configuration
		}
		
		if(REPLACE_NOREF){
			for(int chrom=minChrom; chrom<=maxChrom; chrom++){
				ChromosomeArray cha=Data.getChromosome(chrom);
				final byte[] array=cha.array;
				final byte n='N';
				for(int i=0; i<array.length; i++){
					if(array[i]==n){
						array[i]=AminoAcid.numberToBase[rr.randyNoref.nextInt()&3];
					}
				}
			}
		}
		
		if(PERFECT_READ_RATIO>=1){
			snpRate=insRate=delRate=subRate=0;
			maxSnps=maxInss=maxDels=maxSubs=maxNs=0;
		}

		if(delRate<=0 || maxDelLen<=0 || maxDels<=0){
			delRate=0;
			maxDelLen=minDelLen=maxDels=0;
		}
		if(insRate<=0 || maxInsLen<=0 || maxInss<=0){
			insRate=0;
			maxInsLen=minInsLen=maxInss=0;
		}
		if(subRate<=0 || maxSubLen<=0 || maxSubs<=0){
			subRate=0;
			maxSubLen=minSubLen=maxSubs=0;
		}
		if(snpRate<=0 || maxSnps<=0){
			snpRate=0;
			maxSnps=0;
		}
		if(nRate<=0 || maxNLen<=0 || maxNs<=0){
			nRate=0;
			maxNLen=minNLen=maxNs=0;
		}
		
		System.err.println("snpRate="+snpRate+", max="+maxSnps+", unique="+USE_UNIQUE_SNPS);
		System.err.println("insRate="+insRate+", max="+maxInss+", len=("+minInsLen+"-"+maxInsLen+")");
		System.err.println("delRate="+delRate+", max="+maxDels+", len=("+minDelLen+"-"+maxDelLen+")");
		System.err.println("subRate="+subRate+", max="+maxSubs+", len=("+minSubLen+"-"+maxSubLen+")");
		System.err.println("nRate  ="+nRate+", max="+maxNs+", len=("+minNLen+"-"+maxNLen+")");
		System.err.println("genome="+Data.GENOME_BUILD);
		System.err.println("PERFECT_READ_RATIO="+PERFECT_READ_RATIO);
		System.err.println("ADD_ERRORS_FROM_QUALITY="+ADD_ERRORS_FROM_QUALITY);
		System.err.println("REPLACE_NOREF="+REPLACE_NOREF);
		System.err.println("paired="+paired);
		System.err.println("read length="+(minlen==maxlen ? ""+minlen : minlen+"-"+maxlen));
		if(paired){
			System.err.println("insert size="+(mateMiddleMin+2*minlen)+"-"+(mateMiddleMax+2*maxlen));
		}
		
//		assert(false) : OUTPUT_INTERLEAVED;
		String fname1="reads_B"+Data.GENOME_BUILD+"_"+number+"x"+maxlen+"bp_"
			+(maxSnps==0 || snpRate==0 ? 0 : maxSnps)+"S_"+(maxInss==0 || insRate==0 ? 0 : +maxInsLen)+"I_"+(maxDels==0 || delRate==0 ? 0 : maxDelLen)+"D_"+
			(maxSubs==0 || subRate==0 ? 0 : maxSubLen)+"U_"+
			(maxNs==0 || nRate==0 ? 0 : maxNs)+"N"/*+"_chr"+minChrom+"-"+maxChrom*/+(paired ? (OUTPUT_INTERLEAVED ? "_interleaved" : "_1") : "")+fileExt;
		
		String fname2=(!paired || OUTPUT_INTERLEAVED) ? null : "reads_B"+Data.GENOME_BUILD+"_"+number+"x"+maxlen+"bp_"
			+(maxSnps==0 || snpRate==0 ? 0 : maxSnps)+"S_"+(maxInss==0 || insRate==0 ? 0 : +maxInsLen)+"I_"+(maxDels==0 || delRate==0 ? 0 : maxDelLen)+"D_"+
			(maxSubs==0 || subRate==0 ? 0 : maxSubLen)+"U_"+
			(maxNs==0 || nRate==0 ? 0 : maxNs)+"N"/*+"_chr"+minChrom+"-"+maxChrom*/+"_2"+fileExt;
		
		if(out!=null){
			fname1=out.replaceFirst("#", "1");
			fname2=(!out.contains("#") || !paired || OUTPUT_INTERLEAVED) ? null : out.replaceFirst("#", "2");
		}
		
		rr.writeRandomReadsX(number, minlen, maxlen,
				maxSnps, maxInss, maxDels, maxSubs, maxNs,
				snpRate, insRate, delRate, subRate, nRate,
				minInsLen, minDelLen, minSubLen, minNLen,
				maxInsLen, maxDelLen, maxSubLen, maxNLen,
				minChrom, maxChrom, colorspace, minQuality, midQuality, maxQuality, fname1, fname2);
		
		t.stop();
		Data.sysout.println("Wrote "+fname1);
		if(fname2!=null){Data.sysout.println("Wrote "+fname2);}
		Data.sysout.println("Time: \t"+t);
		
	}
	
	
	private static ArrayList<ChromosomeArray> writeRef(String reference, int build){
		ArrayList<ChromosomeArray> chromlist=null;
		if(reference!=null){
			{
				File f=new File(reference);
				if(!f.exists() || !f.isFile() || !f.canRead()){throw new RuntimeException("Cannot read file "+f.getAbsolutePath());}
			}
			{
				String s=align2.IndexMaker4.fname(1, 1, 13, 1, false);
				String dir=new File(s).getParent();
				dir=dir.replace('\\', '/');
				dir=dir.replace("ref/index/", "ref/genome/");
				String sf=dir+"/summary.txt";
				if(!NODISK && new File(sf).exists() && SummaryFile.compare(sf, reference)){
					//do nothing
					System.err.println("NOTE:\tIgnoring reference file because it already appears to have been processed.");
					System.err.println("NOTE:\tIf you wish to regenerate the index, please manually delete "+dir+"/summary.txt");
					return null;
				}
				File f=new File(dir);
				if(f.exists()){
					File[] f2=f.listFiles();
					if(f2!=null && f2.length>0){
						if(OVERWRITE){
							Data.sysout.println("NOTE:\tDeleting contents of "+dir+" because reference is specified and overwrite="+OVERWRITE);
							for(File f3 : f2){
								if(f3.isFile()){
									String f3n=f3.getName();
									if((f3n.contains(".chrom") || f3n.endsWith(".txt") || f3n.endsWith(".txt.gz")) && !f3n.endsWith("list.txt")){
										f3.delete();
									}
								}
							}
						}else{
							Data.sysout.println(Arrays.toString(f2));
							throw new RuntimeException("\nThere is already a reference at location '"+f.getAbsolutePath()+"'.  " +
									"Please delete it (and the associated index), or use a different build ID, " +
									"or remove the 'reference=' parameter from the command line, or set overwrite=true.");
						}
					}
				}
				dir=dir.replace("ref/genome/", "ref/index/");
				f=new File(dir);
				if(f.exists()){
					File[] f2=f.listFiles();
					if(f2!=null && f2.length>0){
						if(OVERWRITE){
							Data.sysout.println("NOTE:\tDeleting contents of "+dir+" because reference is specified and overwrite="+OVERWRITE);
							for(File f3 : f2){
								if(f3.isFile()){f3.delete();}
							}
						}else{
							throw new RuntimeException("\nThere is already an index at location '"+f.getAbsolutePath()+"'.  " +
									"Please delete it, or use a different build ID, or remove the 'reference=' parameter from the command line.");
						}
					}
				}
			}
			
			Data.sysout.println("Writing reference.");
			
			int oldzl=ReadWrite.ZIPLEVEL;
			ReadWrite.ZIPLEVEL=Tools.max(4, ReadWrite.ZIPLEVEL);
			
			int minScaf=-1;
			int midPad=500;
			int maxChromLen=-1;
			boolean genScaffoldInfo=true;
			
			maxChromLen=maxChromLen>0 ? maxChromLen : FastaToChromArrays.MAX_LENGTH;
			minScaf=minScaf>-1 ? minScaf : FastaToChromArrays.MIN_SCAFFOLD;
			midPad=midPad>-1 ? midPad : FastaToChromArrays.MID_PADDING;
			String[] ftcaArgs=new String[] {reference, ""+build, "writeinthread=false", "genscaffoldinfo="+genScaffoldInfo, "retain", "waitforwriting=false",
					"gzip="+(Data.CHROMGZ), "chromc="+Data.CHROMC, "maxlen="+maxChromLen,
							"writechroms="+(!NODISK), "minscaf="+minScaf, "midpad="+midPad, "nodisk="+NODISK};
			
			chromlist=FastaToChromArrays.main2(ftcaArgs);
			
			ReadWrite.ZIPLEVEL=oldzl;
		}
		return chromlist;
	}
	
	
	public RandomReads3(boolean paired_){
		this(getSeed(), paired_);
	}
	
	public RandomReads3(long seed, boolean paired_){
		if(randomChrom==null){
			synchronized(getClass()){
				if(randomChrom==null){
					randomChrom=fillRandomChrom();
				}
			}
		}
		randy=new Random(seed+1);
		randy2=new Random(seed+2);
		randyMutationType=new Random(seed+3);
		randyCSError=new Random(seed+4);
		randyQual=new Random(seed+5);
		randyAdapter=new Random(seed+25);
		paired=paired_;

		randyPerfectRead=new Random(seed+20);
		randyLength=new Random(seed+21);
		randyAmp=new Random(seed+22);
		
		if(paired){
			randyMate=new Random(seed+6);
			randy2Mate=new Random(seed+7);
			randyMutationTypeMate=new Random(seed+8);
			randyCSErrorMate=new Random(seed+9);
			randyQualMate=new Random(seed+10);
			randyAdapterMate=new Random(seed+30);
		}else{
			randyMate=null;
			randy2Mate=null;
			randyMutationTypeMate=null;
			randyCSErrorMate=null;
			randyQualMate=null;
			randyAdapterMate=null;
		}
		
		if(REPLACE_NOREF){
			randyNoref=new Random(seed+31);
		}else{
			randyNoref=null;
		}
	}
	
	private final void addErrorsFromQuality(Read r, Random randy){
		for(int i=0; i<r.quality.length; i++){
			if(r.bases[i]!='N'){
				float prob=QualityTools.PROB_ERROR[r.quality[i]];
				if(randy.nextFloat()<prob){
					byte old=r.bases[i];
					if(!r.colorspace()){old=AminoAcid.baseToNumber[old];}
					byte b=(byte)((old+randy.nextInt(3)+1)%4);
					if(!r.colorspace()){b=AminoAcid.numberToBase[b];}
					assert(old!=b);
					r.bases[i]=b;
				}
			}
		}
	}

	public byte[] addAdapter(byte[] bases, int[] locs, int readlen, Random rand, byte[] adapter){
//		Data.sysout.println("Adding adapter "+new String(adapter));
		assert(readlen<=bases.length);
		int mod=Tools.max((readlen+1)/2, readlen-30-adapter.length);
		int index=rand.nextInt(mod);
		index-=adapter.length/2;
		for(int i=0, j=index; i<adapter.length; i++, j++){
			if(j>=0 && j<bases.length){bases[j]=adapter[i];}
		}
		return bases;
	}

	public byte[] addSNP(byte[] bases, int[] locs, int readlen, Random rand){
		assert(readlen<=bases.length);
		int index=rand.nextInt(readlen);
		byte old=bases[index];
		byte oldNum=AminoAcid.baseToNumber[old];
		if(oldNum<0){oldNum=0;}
		int num;
		if(BIASED_SNPS && rand.nextInt(3)>0){
			num=(oldNum^3);
		}else{
			num=(oldNum+rand.nextInt(3)+1)%4;
		}
		assert(num>=0 && num<=3 && num!=oldNum);
		bases[index]=AminoAcid.numberToBase[num];
		return bases;
	}

	public byte[] addSNP(byte[] bases, int[] locs, int readlen, Random rand, BitSet bits){
		assert(readlen<=bases.length);
		int index=rand.nextInt(readlen);
		
		while(bits.get(index)){
			index=rand.nextInt(readlen);
		}
		bits.set(index);
		
		byte old=bases[index];
		byte oldNum=AminoAcid.baseToNumber[old];
		if(oldNum<0){oldNum=0;}
		int num;
		if(BIASED_SNPS && rand.nextInt(3)>0){
			num=(oldNum^3);
		}else{
			num=(oldNum+rand.nextInt(3)+1)%4;
		}
		assert(num>=0 && num<=3 && num!=oldNum) : num+", "+oldNum;
		bases[index]=AminoAcid.numberToBase[num];
		return bases;
	}
	

	public byte[] addSUB(byte[] bases, int[] locs, int minlen, int maxlen, int readlen, Random rand){
		assert(readlen<=bases.length) : readlen+", "+bases.length;
		assert(minlen>=1);
		assert(maxlen>=minlen);
		
//		int len=minlen+randy2.nextInt(maxlen-minlen+1);
		int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));
//		int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));

		assert(len>=minlen);
		assert(len<=maxlen);
		
//		System.err.println(minlen+", "+maxlen+", "+readlen+", "+s.length());
		
		int index=rand.nextInt(readlen-len+1);
		
		int lim=index+len-1;
		
		{//Change first and last to anything except old
			int i=index;
			
			byte old=bases[i];
			if(AminoAcid.isFullyDefined(old)){
				byte oldNum=AminoAcid.baseToNumber[old];
				int num=(oldNum+rand.nextInt(4))%4;
				assert(num>=0 && num<=3);
				byte base=AminoAcid.numberToBase[num];
				bases[i]=base;
			}
			
			i=lim;
			old=bases[i];
			if(AminoAcid.isFullyDefined(old)){
				byte oldNum=AminoAcid.baseToNumber[old];
				int num=(oldNum+rand.nextInt(4))%4;
				assert(num>=0 && num<=3);
				byte base=AminoAcid.numberToBase[num];
				bases[i]=base;
			}
		}
		
		for(int i=index+1; i<lim; i++){ //Change middles to anything
			byte old=bases[i];
			if(AminoAcid.isFullyDefined(old)){
				byte oldNum=AminoAcid.baseToNumber[old];
				int num=(oldNum+rand.nextInt(3)+1)%4;
				assert(num>=0 && num<=3 && num!=oldNum);
				byte base=AminoAcid.numberToBase[num];
				bases[i]=base;
			}
		}
		return bases;
	}
	

	public byte[] addN(byte[] bases, int[] locs, int minlen, int maxlen, int readlen, Random rand, BitSet bits){
		assert(readlen<=bases.length) : readlen+", "+bases.length;
		assert(minlen>=1);
		assert(maxlen>=minlen);
		
//		int len=minlen+randy2.nextInt(maxlen-minlen+1);
		int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));

		assert(len>=minlen);
		assert(len<=maxlen);
		
//		System.err.println(minlen+", "+maxlen+", "+readlen+", "+s.length());
		
		int index=rand.nextInt(readlen-len+1);
		if(bits!=null){
			int trials=40;
			while(bits.get(index) && (trials--)>0){
				index=rand.nextInt(readlen-len+1);
			}
			bits.set(index);
		}
		
		int lim=index+len-1;
		
		for(int i=index; i<=lim; i++){bases[i]='N';}

		return bases;
	}
	
	public byte[] addInsertion(byte[] bases, int[] locs, int minlen, int maxlen, int readlen, int[] dif, Random rand){
		assert(readlen<=bases.length) : readlen+", "+bases.length;
		assert(minlen>0);
		assert(maxlen>=minlen);
		
//		int len=minlen+randy2.nextInt(maxlen-minlen+1);
		int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));
//		int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));
		
		len=Tools.min(len, readlen-dif[1]-2);
		if(len<1){return bases;}
		
		dif[0]-=len;
		dif[1]+=len;
		
		int index=rand.nextInt(readlen-len+1); //Assures that all inserted bases will be within the read
		
//		System.err.println("Added insertion "+len+" at "+index);
		
		byte[] bases2=new byte[bases.length+len];
		for(int i=0; i<index; i++){bases2[i]=bases[i];}
		
		for(int i=bases.length-1, j=bases2.length-1; i>=index; i--, j--){
//			if(verbose){
//				System.err.println("i="+i+", bases.length="+bases.length+", j="+j+", bases2.length="+bases2.length+", locs.length="+locs.length+"\n"+Arrays.toString(locs));
//			}
			if(j<locs.length){locs[j]=locs[i];}
			bases2[j]=bases[i];
		}
		
//		for(int i=bases.length-1; i>=index; i--){
//			bases2[i+len]=bases[i];
////			locs[i+len]=locs[i];
//		}
//		final int locfill=locs[(index==0 ? 0 : index-1)];
		for(int i=index; i<index+len; i++){
			int x=rand.nextInt(4);
			byte b=AminoAcid.numberToBase[x];
			bases2[i]=b;
//			locs[i]=locfill;
			locs[i]=-1;
		}
		
		return bases2;
	}
	
	public int[] makeDelsa(int dels, int minlen, int maxlen, Random rand){
		if(dels<1){return null;}
		assert(minlen>0);
		assert(maxlen>=minlen);
		int[] delsa=new int[dels];
		for(int i=0; i<delsa.length; i++){
//			int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));
			int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));
			delsa[i]=len;
		}
		return delsa;
	}
	
	public byte[] addDeletion(byte[] bases, int[] locs, int len, int readlen, int[] dif, Random rand){
		assert(bases.length>=readlen+len) : "bases.len="+bases.length+", readlen="+readlen+", len="+len+", dif="+Arrays.toString(dif);
		assert(len>0);
		
		dif[0]+=len;
		
//		int index=randy2.nextInt(s.length()-len);
		int index=1+rand.nextInt(readlen-1); //Assures there will never be a deletion of the first base, which would not technically be a deletion.
		
//		System.err.println("Added deletion "+len+" at "+index);
		
		byte[] bases2=new byte[bases.length-len];
		for(int i=0; i<index; i++){bases2[i]=bases[i];}
		for(int i=index; i<bases2.length; i++){
			bases2[i]=bases[i+len];
			locs[i]=locs[i+len];
		}
		
		return bases2;
	}
	
	public int randomChrom(Read r0, int minChrom, int maxChrom){
		if(r0!=null){return r0.chrom;}
		
		int x=-1;
		int chrom=-1;
		
		assert(minChrom<=maxChrom) : minChrom+", "+maxChrom;
		while(chrom<minChrom || chrom>maxChrom){
			x=randy.nextInt();
			chrom=randomChrom[(x&0x7FFFFFFF)%randomChrom.length];
		}
		return chrom;
	}
	
	public int randomStrand(Read r0, int minChrom, int maxChrom, boolean sameStrandMate){
		if(r0!=null){
			return sameStrandMate ? r0.strand() : r0.strand()^1;
		}
		return randy.nextInt()&1;
	}
	
	public int randomLoc(Read r0, int chrom, int readlen, int minMiddle, int maxMiddle, int strand){
		
		if(r0!=null){
			return randomLocPaired(r0, chrom, readlen, minMiddle, maxMiddle, strand);
		}
		return randomLocSingle(chrom, readlen);
	}
	
	public int randomLocPaired(Read r0, int chrom, int readlen, int minMiddle, int maxMiddle, int strand){
		assert(r0!=null);
		
		final int midRange=maxMiddle-minMiddle+1;
		int middle=(maxMiddle+minMiddle)/2;
		if(SUPERFLAT_DIST){
			//		Data.sysout.print(other.numericID);
			middle=((int)(r0.numericID%midRange))+minMiddle;
			//		Data.sysout.println("\t"+middle);
		}else if(FLAT_DIST){
			middle=randyMate.nextInt(midRange)+minMiddle;
		}else if(BELL_DIST){
			assert(false) : "TODO";
		}else{
			middle=(randyMate.nextInt(midRange)+randyMate.nextInt(midRange))/2+minMiddle;
		}

		//	Data.sysout.println(sameStrand+": "+other.strand+" -> "+strand);
		int x;
		if(r0.strand()==Gene.PLUS){
			x=r0.stop+middle+1;
		}else{
			x=r0.start-middle-readlen;
		}
		return x;
	}
	
	public int randomLocSingle(int chrom, int readlen){
		ChromosomeArray cha=Data.getChromosome(chrom);
		byte[] array=cha.array;
		if(readlen>=(cha.maxIndex-cha.minIndex)){return -1;}
		
		int loc=-1;
		for(int i=0; loc<0 && i<24; i++){
			loc=randy.nextInt(cha.maxIndex-readlen);
			for(int j=0; j<readlen; j++){
				if(!AminoAcid.isFullyDefined(array[j+loc])){
					loc=-1;
					break;
				}
			}
		}
		return loc;
	}
	
	public void writeRandomReadsX(int numReads, int minlen, int maxlen,
			int maxSnps, int maxInss, int maxDels, int maxSubs, int maxNs,
			float snpRate, float insRate, float delRate, float subRate, float nRate,
			int minInsLen, int minDelLen, int minSubLen, int minNLen,
			int maxInsLen, int maxDelLen, int maxSubLen, int maxNLen, 
			int minChrom, int maxChrom, boolean colorspace,
			int minQual, int midQual, int maxQual, String fname1, String fname2){
		FASTQ.PARSE_CUSTOM=true;
		
		TextStreamWriter tsw1=new TextStreamWriter(fname1, OVERWRITE, false, true);
		tsw1.start();
		TextStreamWriter tsw2=null;
		if(fname2!=null){
			assert(!fname2.equalsIgnoreCase(fname1));
			tsw2=new TextStreamWriter(fname2, OVERWRITE, false, true);
			tsw2.start();
		}
		
		assert(minQual<=midQual);
		assert(midQual<=maxQual);
		assert(minQual>=0 && maxQual<60);
		
		final int maxQualP=maxQual;//Tools.max(35, maxQual);
		final int midQualP=midQual;//30;
		final int minQualP=minQual;//Tools.min(25, maxQual);
		
		final BitSet bits=new BitSet(maxlen+1);
		final int[] locs=new int[(int)Tools.min(300000000, maxlen+(maxDelLen*(long)maxDels))];
		
		Read lastRead=null;
		int ampLevel=0;
		int ampLength=2000;
		
		for(int i=0; i<numReads; i++){
			
			final boolean perfect=randyPerfectRead.nextFloat()<PERFECT_READ_RATIO;
			
			final byte baseQuality;
			final byte slant;
			{
				byte baseSlant=(perfect ? (byte)5 : (byte)(maxQual-minQual+1));
				slant=(byte)((randyQual.nextInt(baseSlant)+randyQual.nextInt(baseSlant)+1)/2);
				if(randyQual.nextBoolean()){
					int range=(perfect ? maxQualP-midQualP+1 : maxQual-midQual+1);
					int delta=Tools.min(randyQual.nextInt(range), randyQual.nextInt(range));
					baseQuality=(byte)((perfect ? midQualP : midQual)+delta);
				}else{
					int range=perfect ? midQualP-minQualP+1 : midQual-minQual+1;
					int delta=randyQual.nextInt(range);
					baseQuality=(byte)((perfect ? midQualP : midQual)-delta);
				}
			}
			
			int forceChrom=-1, forceLoc=-1;
			if(AMP>1 && lastRead!=null){
				if(ampLevel>0){
					forceChrom=lastRead.chrom;
//					forceLoc=lastRead.start+4-randyAmp.nextInt(9);
//					forceLoc=lastRead.start+10-randyAmp.nextInt(21);
					
					int a=ampLength;
					int b=a*2+1;
					if(randyAmp.nextBoolean()){
						forceLoc=lastRead.start+a-randyAmp.nextInt(b);
					}else{
//						if(randyAmp.nextBoolean()){
//							forceLoc=lastRead.start+a-(randyAmp.nextInt(b)+randyAmp.nextInt(b))/2;
//						}else{
							forceLoc=lastRead.start+a-(randyAmp.nextInt(b)+randyAmp.nextInt(b)+randyAmp.nextInt(b))/3;
//						}
					}
				}else{
					
					int a1=AMP;
					if(randyAmp.nextInt(30)==0){a1*=7;}
					
					if(randyAmp.nextInt(3)>0){
						ampLevel=Tools.min(randyAmp.nextInt(a1), randyAmp.nextInt(a1));
					}else{
						double log=Math.log10(a1*7);
						ampLevel=(int)Math.round(Math.pow(10, randyAmp.nextDouble()*log));
					}
					
					ampLength=500+randyAmp.nextInt(3001);
//					ampLevel=randyAmp.nextInt(AMP);
				}
			}
			
			Read r1=makeRead(null, minlen, maxlen, minChrom, maxChrom,
					maxSnps, maxInss, maxDels, maxSubs, maxNs,
					snpRate, insRate, delRate, subRate, nRate,
					minInsLen, minDelLen, minSubLen, minNLen,
					maxInsLen, maxDelLen, maxSubLen, maxNLen, 
					mateMiddleMin, mateMiddleMax, mateSameStrand, 
					minQual, midQual, maxQual, baseQuality, slant, 
					perfect, colorspace, nextReadID, locs, bits, forceChrom, forceLoc);

//			assert(false) : r1;
			if(paired && r1!=null){

				Read r2=null;
				for(int tries=0; r2==null && tries<100; tries++){
					r2=makeRead(r1, minlen, maxlen, minChrom, maxChrom,
							maxSnps, maxInss, maxDels, maxSubs, maxNs,
							snpRate, insRate, delRate, subRate, nRate,
							minInsLen, minDelLen, minSubLen, minNLen,
							maxInsLen, maxDelLen, maxSubLen, maxNLen, 
							mateMiddleMin, mateMiddleMax, mateSameStrand, 
							minQual, midQual, maxQual, baseQuality, slant, 
							perfect, colorspace, nextReadID, locs, bits, -1, -1);
				}
				
				if(r2!=null){
					r1.mate=r2;
					r2.mate=r1;
				}else{
					r1=null;
				}
				
//				Data.sysout.println(r.strand()+"\t"+r.insertSize());
			}
			if(r1!=null){
//				assert(false) : r1;
				tsw1.println(r1);
				if(r1.mate!=null){
					r1.mate.setPairnum(1);
					if(tsw2!=null){tsw2.println(r1.mate);}
					else{tsw1.println(r1.mate);}
						
				}
				nextReadID++;
			}else{
				i--;
			}
			ampLevel=Tools.max(0, ampLevel-1);
			if(ampLevel==0){lastRead=null;}

			if(lastRead==null){lastRead=r1;}
//			System.err.println("Made "+r.start+" ~ "+r.stop+" = "+(r.stop-r.start));
		}
		tsw1.poison();
		if(tsw2!=null){tsw2.poison();}
	}
	
	public Read makeRead(Read r0, int minlen, int maxlen, int minChrom, int maxChrom,
			int maxSnps, int maxInss, int maxDels, int maxSubs, int maxNs,
			float snpRate, float insRate, float delRate, float subRate, float nRate,
			int minInsLen, int minDelLen, int minSubLen, int minNLen,
			int maxInsLen, int maxDelLen, int maxSubLen, int maxNLen, 
			int minMiddle, int maxMiddle, boolean sameStrand, 
			int minQual, int midQual, int maxQual, byte baseQuality, byte slant, 
			boolean perfect, boolean colorspace, long rid, int[] locs, BitSet bits,
			int FORCE_CHROM, int FORCE_LOC){
		
//		verbose=(rid==3860);
		
		int SNPs=0;
		int INSs=0;
		int DELs=0;
		int SUBs=0;
		int Ns=0;
		int adapters=0;
			
		while(SNPs<maxSnps && randyMutationType.nextFloat()<snpRate){SNPs++;}
		while(INSs<maxInss && randyMutationType.nextFloat()<insRate){INSs++;}
		while(DELs<maxDels && randyMutationType.nextFloat()<delRate){DELs++;}
		while(SUBs<maxSubs && randyMutationType.nextFloat()<subRate){SUBs++;}
		while(Ns<maxNs && randyMutationType.nextFloat()<nRate){Ns++;}
		
//		final boolean perfect=randyPerfectRead.nextFloat()<PERFECT_READ_RATIO;
		if(perfect){SNPs=INSs=DELs=SUBs=Ns=0;}
		
		if(verbose){
			System.err.println("\nMaking read with snps="+SNPs+", inss="+INSs+", dels="+DELs+", subs="+SUBs+", Ns="+Ns);
			System.err.println("perfect="+perfect);
		}
		
		int[] delsa=makeDelsa(DELs, minDelLen, maxDelLen, randy2);
		
		final int readlen=(minlen==maxlen ? maxlen : minlen+randyLength.nextInt(maxlen-minlen+1));
		int inititallen0=readlen+(delsa==null ? 0 : (int)Tools.sum(delsa));
		
		if(verbose){
			System.err.println("delsa="+Arrays.toString(delsa));
			System.err.println("readlen="+readlen+", inititallen0="+inititallen0);
		}
		
		final int chrom=(FORCE_CHROM>=0 ? FORCE_CHROM : randomChrom(r0, minChrom, maxChrom));
		if(chrom<0){return null;}
		final int strand=randomStrand(r0, minChrom, maxChrom, sameStrand);
		int loc=(FORCE_LOC>=0 ? FORCE_LOC : randomLoc(r0, chrom, inititallen0, minMiddle, maxMiddle, strand));
		
		if(verbose){
			System.err.println("chrom="+chrom+", loc="+loc+"~"+(loc+inititallen0-1)+", strand="+strand+", chalen="+Data.getChromosome(chrom).maxIndex);
		}
		
		if(r0!=null){
			int y=loc+inititallen0-1;
			if(loc<0){loc=0; y=readlen-1; maxDels=0; delsa=null; inititallen0=readlen;}
			ChromosomeArray cha=Data.getChromosome(chrom);
			if(y>cha.maxIndex){y=cha.maxIndex; loc=y-readlen+1; maxDels=0; delsa=null; inititallen0=readlen;}
			if(verbose){
				System.err.println("After pair compensation:");
				System.err.println("delsa="+Arrays.toString(delsa));
				System.err.println("readlen="+readlen+", inititallen0="+inititallen0);
				System.err.println("chrom="+chrom+", loc="+loc+", strand="+strand);
			}
			assert(y<=cha.maxIndex) : y+", "+cha.maxIndex;
			assert(cha.get(y)>0) : cha.get(y);
		}
		
		if(loc<0){
			if(verbose){
				System.err.println("Bad values; returning null.");
			}
			return null;
		}
		
		final ChromosomeArray cha=Data.getChromosome(chrom);
		if(readlen>=(cha.maxIndex-cha.minIndex)){
			if(verbose){
				System.err.println("Too long; returning null.");
			}
			return null;
		}
		if(loc>=cha.maxIndex || loc<0){return null;}
		byte[] bases=cha.getBytes(loc, loc+inititallen0-1);
		assert(bases[0]>0 && bases[bases.length-1]>0) : Arrays.toString(bases);
		assert(strand==Gene.MINUS || strand==Gene.PLUS);
		
		for(int i=0; i<bases.length; i++){
			locs[i]=i+loc;
		}
		if(verbose){
			System.err.println(new String(bases));
			System.err.println(Arrays.toString(Arrays.copyOf(locs, bases.length)));
		}
		
		if(adapter1!=null && (rid&3)==0){
			bases=addAdapter(bases, locs, readlen, randyAdapter, adapter1);
			adapters++;
		}else if(adapter2!=null && (rid&3)==1){
			bases=addAdapter(bases, locs, readlen, randyAdapter, adapter2);
			adapters++;
		}
		
		int[] dif=new int[] {0, 0};
		
		for(int j=0; delsa!=null && j<delsa.length; j++){
			bases=addDeletion(bases, locs, delsa[j], readlen, dif, randy2);
			if(verbose){
				System.err.println("After adding del "+delsa[j]+": ");
				System.err.println(new String(bases));
				System.err.println(Arrays.toString(Arrays.copyOf(locs, bases.length)));
			}
		}
		if(bases.length>readlen){bases=Arrays.copyOf(bases, readlen);}
		assert(bases.length==readlen);
		
		for(int j=0; j<INSs; j++){
			bases=addInsertion(bases, locs, minInsLen, maxInsLen, readlen, dif, randy2);
			if(verbose){
				System.err.println("After adding ins: ");
				System.err.println("'"+new String(bases)+"'");
				System.err.println(Arrays.toString(Arrays.copyOf(locs, Tools.min(locs.length, bases.length))));
			}
		}
		if(bases.length!=readlen){bases=Arrays.copyOf(bases, readlen);}
		assert(bases.length==readlen);
		
		if(USE_UNIQUE_SNPS){
			bits.clear();
			for(int j=0; j<SNPs; j++){bases=addSNP(bases, locs, readlen, randy2, bits);}
		}else{
			for(int j=0; j<SNPs; j++){bases=addSNP(bases, locs, readlen, randy2);}
		}
		
		for(int j=0; j<SUBs; j++){bases=addSUB(bases, locs, minSubLen, maxSubLen, readlen, randy2);}
		
		if(USE_UNIQUE_SNPS){
			bits.clear();
			for(int j=0; j<Ns; j++){bases=addN(bases, locs, minNLen, maxNLen, readlen, randy2, bits);}
		}else{
			for(int j=0; j<Ns; j++){bases=addN(bases, locs, minNLen, maxNLen, readlen, randy2, null);}
		}
		
		//Fill insertions in loc array
		for(int i=1; i<bases.length; i++){
			if(locs[i]<0){locs[i]=locs[i-1];}
		}
		for(int i=bases.length-2; i>=0; i--){
			if(locs[i]<0){locs[i]=locs[i+1];}
		}
		final int x=locs[0], y=locs[bases.length-1];
		if(verbose){
			System.err.println("After adding SNPs, SUBs, Ns, and fixing locs: ");
			System.err.println("'"+new String(bases)+"'");
			System.err.println(Arrays.toString(Arrays.copyOf(locs, Tools.min(locs.length, bases.length))));
		}
		
//		if(FORCE_LOC>=0 || FORCE_CHROM>=0){
//			if(y<0 || y+readlen>)
//		}
		assert(FORCE_LOC>=0 || FORCE_CHROM>=0 || y<=cha.maxIndex) : y+", "+r0;
		assert(FORCE_LOC>=0 || FORCE_CHROM>=0 || cha.get(y)>0) : cha.get(y);
		
		if(strand==Gene.MINUS){
			AminoAcid.reverseComplementBasesInPlace(bases);
			//Reverse loc array; not really necessary
			for(int i=0, lim=bases.length/2; i<lim; i++){
				int tmp=locs[i];
				locs[i]=locs[bases.length-i-1];
				locs[bases.length-i-1]=tmp;
			}
			if(verbose){
				System.err.println("After reverse-complement: ");
				System.err.println(new String(bases));
				System.err.println(Arrays.toString(Arrays.copyOf(locs, bases.length)));
			}
		}
		
		if(verbose){
			System.err.println("Final lineup: ");
			System.err.println(new String(bases));
			for(int i=0; i<bases.length; i++){
				byte c=cha.get(locs[i]);
				if(strand==1){c=AminoAcid.baseToComplementExtended[c];}
				System.err.print((char)c);
			}
			System.err.println();
		}
		
		byte[] quals=null;
		if(USE_FIXED_QUALITY){
			quals=getFixedQualityRead(bases.length);
		}else{
			if(perfect){
				quals=QualityTools.makeQualityArray(bases.length, randyQual, 30, 40, baseQuality, slant);
			}else{
				quals=QualityTools.makeQualityArray(bases.length, randyQual, minQual, maxQual, baseQuality, slant);
			}
		}
		for(int j=0; j<quals.length; j++){
			if(!AminoAcid.isFullyDefined(bases[j])){quals[j]=0;}
		}
		
		
//		Read r=new Read(bases, chrom, (byte)strand, loc, loc+bases.length-1, rid, quals, false);
		Read r=new Read(bases, chrom, (byte)strand, x, y, rid, quals, false);
		r.setSynthetic(true);
		assert(r.bases.length==readlen);
		
		if(ADD_ERRORS_FROM_QUALITY && !perfect){addErrorsFromQuality(r, randyQual);}
		
		assert(r.bases.length==readlen);
		
//		r.stop=r.start+readlen+dif[0]-1;
		
		assert(r.stop>r.start) : r;
		
		if(colorspace){
			r=r.translateToColorspace(true);
			r.obj=new String(r.bases); //TODO - for testing
		}
		r.mapLength=r.bases.length;
		if(adapters>0){r.setHasAdapter(true);}
		if(FORCE_SINGLE_SCAFFOLD && !Data.isSingleScaffold(r.chrom, r.start, r.stop)){return null;}
		if(MIN_SCAFFOLD_OVERLAP>0 && Data.scaffoldOverlapLength(r.chrom, r.start, r.stop)<MIN_SCAFFOLD_OVERLAP){return null;}
		return r;
	}
	
	
	public void addColorspaceErrors(Read r, int errors){
		assert(r.colorspace());
		for(int i=0; i<errors; i++){
			int loc=randyCSError.nextInt(r.bases.length);
			if(r.bases[loc]!='N'){
				assert(r.quality[loc]>=4);
				r.quality[loc]=(byte) (4+randyCSError.nextInt(r.quality[loc]));
				r.bases[loc]=(byte)((r.bases[loc]+randyCSError.nextInt(3)+1)&3);
			}
		}
	}
	
	private static int[] randomChrom;

	private static int[] fillRandomChrom(){

		int[] in=Arrays.copyOf(Data.chromLengths, Data.chromLengths.length);
		long total=Tools.sum(in);
		int div=(int)(total/8192);
		for(int i=0; i<in.length; i++){in[i]=((in[i]+div-1)/div);}


		int sum=0;
		for(int i=0; i<in.length; i++){sum+=in[i];}
		int[] out=new int[sum];
		sum=0;
		for(int chrom=0; chrom<in.length; chrom++){
			int size=in[chrom];
			for(int j=0; j<size; j++){
				out[sum+j]=chrom;
			}
			sum+=size;
		}
		return out;
	}
	
	public static final byte[] getFixedQualityRead(int bases){
		if(fixedQuality[bases]==null){
			fixedQuality[bases]=new byte[bases];
			Arrays.fill(fixedQuality[bases], FIXED_QUALITY_VALUE);
		}
		return fixedQuality[bases];
	}
	
	private static synchronized long getSeed(){
		long r=seed*1000;
		seed++;
		return r;
	}
	private static long seed=0;

	private final Random randy;
	private final Random randy2;
	private final Random randyMutationType;
	private final Random randyCSError;
	private final Random randyQual;
	private final Random randyAdapter;
	
	private final Random randyMate;
	private final Random randy2Mate;
	private final Random randyMutationTypeMate;
	private final Random randyCSErrorMate;
	private final Random randyQualMate;
	private final Random randyAdapterMate;

	private final Random randyPerfectRead;
	private final Random randyNoref;

	private final Random randyLength;
	
	private final Random randyAmp;
	
	public final boolean paired;
	
	private long nextReadID=0;
	private byte[] adapter1=null;
	private byte[] adapter2=null;
	
	private static final byte[][] fixedQuality=new byte[1000][];
	
	public static final boolean USE_FIXED_QUALITY=false;
	public static final byte FIXED_QUALITY_VALUE=24;
	public static boolean ADD_ERRORS_FROM_QUALITY=true;
	public static boolean REPLACE_NOREF=false;
	public static boolean OUTPUT_INTERLEAVED=false;
	public static String fileExt=".fq.gz";
	public static boolean verbose=false;
	
//	public static int mateLen=35;
	public static boolean mateSameStrand=false;
	public static int mateMiddleMin=-100; //default -25
	public static int mateMiddleMax=100; //default 475
	public static boolean SUPERFLAT_DIST=false;
	public static boolean FLAT_DIST=false;
	public static boolean BELL_DIST=false;
	public static boolean BIASED_SNPS=false;
	
	public static boolean NODISK=false;
	
	public static int AMP=1;
	
	public static float PERFECT_READ_RATIO=0f;

	public static boolean USE_UNIQUE_SNPS=true;
	public static boolean FORCE_SINGLE_SCAFFOLD=true;
	public static int MIN_SCAFFOLD_OVERLAP=1;
	public static boolean OVERWRITE=true;
	
}
