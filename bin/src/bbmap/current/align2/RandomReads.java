package align2;

import java.util.Arrays;
import java.util.Random;

import stream.FASTQ;
import stream.Read;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;

public final class RandomReads {
	
	
	public static void main(String[] args){
		
		Data.GENOME_BUILD=-1;
		
		int number=100000;
		int readlen=100;
		
		boolean paired=false;

//		final int maxInsertionLen=12;
//		final int maxSubLen=12;
//		final int maxDeletionLen=400;

		final int maxInsertionLen=20;
		final int maxSubLen=20;
		final int maxDeletionLen=20;
		
		int minChrom=-1;
		int maxChrom=-1;
		
		int maxSnps=3;
		int maxInss=2;
		int maxDels=2;
		int maxSubs=2;

		float snpRate=0.4f;
		float insRate=0.2f;
		float delRate=0.2f;
		float subRate=0.2f;
		PERFECT_READ_RATIO=0.4f;
		
		int minQuality=30;
		int midQuality=33;
		int maxQuality=35;
		
		
		minQuality=5;
		midQuality=20;
		maxQuality=35;
		
		
		maxSnps=3;//4;
		maxInss=3;//2;
		maxDels=3;
		maxSubs=3;//2;
		
		snpRate=0.25f;
		insRate=0.25f;
		delRate=0.25f;
		subRate=0.25f;//0.3f;
		PERFECT_READ_RATIO=0.2f;//0.2f;//0.8f
		
		
		
		
		for(int i=0; i<args.length; i++){
			String arg=args[i].toLowerCase();
			assert(arg.contains("=")) : "All arguments must be of the form word=number, e.g., reads=10000";
			String[] split=arg.split("=");
			assert(split.length==2);
			final String a=split[0], b=split[1];
			
			int x=-1;
			try {
				x=Integer.parseInt(b);
			} catch (NumberFormatException e) {}

			if(a.equals("reads")){number=x;}
			else if(a.startsWith("len")){readlen=x;}
			else if(a.equals("s") || a.startsWith("snp")){
				maxSnps=x;
				snpRate=1;
			}else if(a.equals("i") || a.startsWith("ins")){
				maxInss=x;
				insRate=1;
			}else if(a.equals("d") || a.startsWith("del")){
				maxDels=x;
				delRate=1;
			}else if(a.equals("u") || a.startsWith("sub")){
				maxSubs=x;
				subRate=1;
			}else if(a.startsWith("minchrom")){
				minChrom=x;
			}else if(a.startsWith("maxchrom")){
				maxChrom=x;
			}else if(a.startsWith("build")){
				Data.setGenome(x);
			}else if(a.startsWith("minq")){
				minQuality=x;
			}else if(a.startsWith("midq")){
				midQuality=x;
			}else if(a.startsWith("maxq")){
				maxQuality=x;
			}else if(a.startsWith("paired")){
				paired=Tools.parseBoolean(b);
			}else if(a.equals("perfect")){PERFECT_READ_RATIO=Float.parseFloat(b);}

		}
		
		assert(Data.GENOME_BUILD>=0);
		if(minChrom<1){minChrom=1;}
		if(maxChrom<1){maxChrom=Data.numChroms;}
		
		RandomReads rr=new RandomReads(paired);
		mateLen=readlen;
		
		System.err.println("snpRate =  \t"+snpRate);
		System.err.println("insRate =  \t"+insRate);
		System.err.println("delRate =  \t"+delRate);
		System.err.println("subRate =  \t"+subRate);
		System.err.println("Reads   =  \t"+number);
		System.err.println("Readlen =  \t"+readlen);
		System.err.println("Paired  =  \t"+paired);
		System.err.println("Genome  =  \t"+Data.GENOME_BUILD);
		System.err.println("PERFECT_READ_RATIO="+PERFECT_READ_RATIO);
		
		String fname="reads_B"+Data.GENOME_BUILD+"_"+number+"x"+readlen+"bp_"
			+maxSnps+"S_"+maxInss+"I_"+maxDels+"D_"+maxSubs+"U_chr"+minChrom+"-"+maxChrom+(paired ? "_#" : "")+".fq";
		
		Read[] reads=rr.makeRandomReadsX(number, readlen, 
				maxSnps, maxInss, maxDels, maxSubs,
				snpRate, insRate, delRate, subRate,
				maxInsertionLen, maxDeletionLen,  maxSubLen,
				minChrom, maxChrom, false, minQuality, midQuality, maxQuality);

		FASTQ.writeFASTQ(reads, fname.replace("#", "1"));
		if(paired){
			for(int i=0; i<reads.length; i++){
				reads[i]=reads[i].mate;
			}
			FASTQ.writeFASTQ(reads, fname.replace("#", "2"));
		}
	}
	
	public RandomReads(boolean paired_){
		this(getSeed(), paired_);
	}
	
	public RandomReads(long seed, boolean paired_){
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
		paired=paired_;
		
		randyPerfectRead=new Random(seed+20);
		
		if(paired){
			randyMate=new Random(seed+6);
			randy2Mate=new Random(seed+7);
			randyMutationTypeMate=new Random(seed+8);
			randyCSErrorMate=new Random(seed+9);
			randyQualMate=new Random(seed+10);
		}else{
			randyMate=null;
			randy2Mate=null;
			randyMutationTypeMate=null;
			randyCSErrorMate=null;
			randyQualMate=null;
		}
	}
	
	private final void addErrorsFromQuality(Read r, Random randy){
		for(int i=0; i<r.quality.length; i++){
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

	public String addSNP(String s, int readlen, Random rand){
		assert(readlen<=s.length());
		int index=rand.nextInt(readlen);
		char old=s.charAt(index);
		byte oldNum=AminoAcid.baseToNumber[old];
		int num=(oldNum+rand.nextInt(3)+1)%4;
		assert(num>=0 && num<=3 && num!=oldNum);
		byte[] bytes=s.getBytes();
		bytes[index]=AminoAcid.numberToBase[num];
		return new String(bytes);
	}
	

	public String addSUB(String s, int minlen, int maxlen, int readlen, Random rand){
		assert(readlen<=s.length()) : readlen+", "+s.length();
		assert(minlen>1);
		assert(maxlen>=minlen);
		
//		int len=minlen+rand.nextInt(maxlen-minlen+1);
		int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));

		assert(len>=minlen);
		assert(len<=maxlen);
		
		assert(readlen<=s.length());
		int index=rand.nextInt(readlen-len+1);
		byte[] bytes=s.getBytes();
		
		int lim=index+len-1;
		
		{//Change first and last to anything except old
			int i=index;
			
			byte old=bytes[i];
			if(AminoAcid.isFullyDefined(old)){
				byte oldNum=AminoAcid.baseToNumber[old];
				int num=(oldNum+rand.nextInt(4))%4;
				assert(num>=0 && num<=3);
				byte base=AminoAcid.numberToBase[num];
				bytes[i]=base;
			}
			
			i=lim;
			old=bytes[i];
			if(AminoAcid.isFullyDefined(old)){
				byte oldNum=AminoAcid.baseToNumber[old];
				int num=(oldNum+rand.nextInt(4))%4;
				assert(num>=0 && num<=3);
				byte base=AminoAcid.numberToBase[num];
				bytes[i]=base;
			}
		}
		
		for(int i=index+1; i<lim; i++){ //Change middles to anything
			byte old=bytes[i];
			if(AminoAcid.isFullyDefined(old)){
				byte oldNum=AminoAcid.baseToNumber[old];
				int num=(oldNum+rand.nextInt(3)+1)%4;
				assert(num>=0 && num<=3 && num!=oldNum);
				byte base=AminoAcid.numberToBase[num];
				bytes[i]=base;
			}
		}
		
		return new String(bytes);
	}
	
	
	public String addInsertion(String s, int minlen, int maxlen, int readlen, int[] dif, Random rand){
		
		
//		assert(false) : minlen+","+maxlen;
		assert(readlen<=s.length()) : readlen+", "+s.length();
		assert(minlen>0);
		assert(maxlen>=minlen);
		
//		int len=minlen+rand.nextInt(maxlen-minlen+1);
		int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));
		
		len=Tools.min(len, readlen-dif[1]-2);
//		assert(false) : len+", "+readlen+", "+dif[1];
		if(len<1){return s;}
		
		if(verbose){System.err.println("\nAdding insertion of len="+len+", dif="+dif[0]);}
		
		dif[0]-=len;
		dif[1]+=len;
		
		int index=rand.nextInt(readlen-len+1); //Assures that all inserted bases will be within the read
		
		String mid="";
		for(int i=0; i<len; i++){
			int x=rand.nextInt(4);
			mid=mid+(char)AminoAcid.numberToBase[x];
		}
		
		if(verbose){System.err.println("dif="+dif[0]+", index="+index+", mid="+mid);}
		
//		System.out.println("Added "+mid.length()+"bp insertion");
//		System.err.println("Added insertion "+len+" at "+index);
		
		String s2=s.substring(0, index)+mid+s.substring(index);
		
		if(verbose){System.err.println(s+"\n->\n"+s2);}
		
		return s2;
	}
	
	
	public String addDeletion(String s, int minlen, int maxlen, int readlen, int[] dif, Random rand){
		assert(s.length()>=readlen+maxlen);
		assert(minlen>0);
		assert(maxlen>=minlen);

//		int len=maxlen;
//		int len=minlen+rand.nextInt(maxlen-minlen+1);
		int len=minlen+Tools.min(rand.nextInt(maxlen-minlen+1), rand.nextInt(maxlen-minlen+1));
//		System.err.println("Made del len "+len);
		dif[0]+=len;
		
//		int index=rand.nextInt(s.length()-len);
		int index=1+rand.nextInt(readlen-1); //Assures there will never be a deletion of the first base, which would not technically be a deletion.
		
//		System.err.println("Added deletion "+len+" at "+index);
		
		String s2=s.substring(0, index)+s.substring(index+len);
		return s2;
	}
	
	
	public Read[] makeRandomReads(int readlen, int number, int minChrom, int maxChrom){
		Read[] out=new Read[number];
		for(int i=0; i<number; i++){
			out[i]=makeRandomRead(readlen, minChrom, maxChrom);
		}
		return out;
	}
	
	
	public Read makeRandomRead(int readlen, int minChrom, int maxChrom){
		int x=-1;
		int chrom=-1;
		
//		System.err.println("Making random read of length "+readlen);
		
		assert(minChrom<=maxChrom);
		while(chrom<minChrom || chrom>maxChrom){
			x=randy.nextInt();
			chrom=randomChrom[(x&0x7FFFFFFF)%randomChrom.length];
//			if(chrom>25 && Data.GENOME_BUILD==36){chrom=-1;}
		}
		byte strand=(byte) (x>=0 ? 0 : 1);
//		strand=0; //TODO
//		System.err.println("Chose chrom "+chrom+", strand "+strand);
		return makeRandomRead2(readlen, chrom, strand);
	}
	
	
	public Read makeRandomRead2(int readlen, int chrom, byte strand){
		byte[] s=null;
		ChromosomeArray cha=Data.getChromosome(chrom);
		
		int loc=-1;
		while(s==null){
//			loc=randy.nextInt(cha.maxIndex-40000);
			loc=randy.nextInt(cha.maxIndex-readlen);
//			loc=10180206;
			s=cha.getBytes(loc, loc+readlen-1);
			assert(s.length==readlen);
			
//			System.out.println(new String(s));
			
			if(AminoAcid.countUndefined(s)>5){
				s=null;
//				System.out.println("Tossed out string.");
			}
		}
//		System.err.println("Chose loc="+loc);
		assert(strand==Gene.MINUS || strand==Gene.PLUS);
		
		if(strand==Gene.MINUS){
			s=AminoAcid.reverseComplementBases(s);
		}
		long id=nextReadID;
		nextReadID++;
		Read r=new Read(s, chrom, strand, loc, loc+s.length-1, id, null, false);
		r.setSynthetic(true);
//		System.err.println("Made read "+r.start+", "+r.stop);
//		assert(readlen==20);
		assert(r.bases.length==readlen);
//		assert(false) : r.start+", "+r.stop;
		return r;
	}
	
	
	
	public Read[] makeRandomReadsX(int numReads, int readlen, 
			int maxSnps, int maxInss, int maxDels, int maxSubs, 
			float snpRate, float insRate, float delRate, float subRate,
			int maxInsertionLen, int maxDeletionLen, int maxSubLen, 
			int minChrom, int maxChrom, boolean colorspace,
			int minQual, int midQual, int maxQual){
		assert(minQual<=midQual);
		assert(midQual<=maxQual);
		assert(minQual>=0 && maxQual<48);
//		System.err.println("Called makeRandomReadsX("+numReads+", "+readlen+", "+maxSnps+", "+maxDels+", "+maxInss+", "+
//			snpRate+", "+delRate+", "+insRate+", "+maxInsertionLen+", "+maxDeletionLen+", "+minChrom+", "+maxChrom+")");
		
//		if(colorspace){readlen++;}
		
		Read[] reads=new Read[numReads];
		
//		assert(Index2.maxIndel==maxIndel); //Temporary
		

		final int maxQualP=Tools.max(35, maxQual);
		final int midQualP=30;
		final int minQualP=Tools.min(25, maxQual);
		
		for(int i=0; i<reads.length; i++){
			
//			System.err.println("Making read "+i);
			Read r=makeRandomRead(readlen+(maxDeletionLen*maxDels), minChrom, maxChrom);
			
			while(r==null){
				r=makeRandomRead(readlen+(maxDeletionLen*maxDels), minChrom, maxChrom);
			}
			
//			assert(r.numericID<=70192);
//			verbose=r.numericID==70192;
			if(verbose){System.err.println(r.header()+"\n"+r+"\n");}
			
			reads[i]=r;
//			System.err.println("B");
//			String s=r.bases;
			String s=new String(r.bases);
			
			if(r.strand()==1){s=AminoAcid.reverseComplementBases(s);}
			if(verbose){System.err.println(s+"\n");}
			
			int SNPs=0;
			int INSs=0;
			int DELs=0;
			int SUBs=0;
				
			while(SNPs<maxSnps && randyMutationType.nextFloat()<snpRate){SNPs++;}
			while(INSs<maxInss && randyMutationType.nextFloat()<insRate){INSs++;}
			while(DELs<maxDels && randyMutationType.nextFloat()<delRate){DELs++;}
			while(SUBs<maxSubs && randyMutationType.nextFloat()<subRate){SUBs++;}

			final boolean perfect=randyPerfectRead.nextFloat()<PERFECT_READ_RATIO;
			if(perfect){SNPs=INSs=DELs=SUBs=0;}
			
			
//			System.err.println("Making read with "+SNPs+", "+INSs+", "+DELs+", "+SUBs);
			
			int[] dif=new int[] {0, 0};
			

			if(verbose){System.err.println("Before "+DELs+" DEL: dif="+dif[0]+"\n"+s+"\n");}
			
			for(int j=0; j<DELs; j++){
//				lengthDif-=s.length();
				s=addDeletion(s, 1, maxDeletionLen, readlen, dif, randy2);
//				lengthDif+=s.length();
			}
			if(verbose){System.err.println("After "+DELs+" DEL: dif="+dif[0]+"\n"+s+"\n");}
			
			if(s.length()>readlen){s=s.substring(0, readlen);}
			assert(s.length()==readlen);

			if(verbose){System.err.println("After length adjust 1 to "+readlen+": dif="+dif[0]+"\n"+s+"\n");}

//			int preInsertLength=s.length();
			for(int j=0; j<INSs; j++){
//				lengthDif-=s.length();
				s=addInsertion(s, 1, maxInsertionLen, readlen, dif, randy2);
//				lengthDif+=s.length();
			}
			if(verbose){System.err.println("After "+INSs+" INS: dif="+dif[0]+"\n"+s+"\n");}
//			int insertLength=s.length()-preInsertLength;
			
			if(s.length()!=readlen){
				assert(s.length()>readlen);
				s=s.substring(0, readlen);
			}
			if(verbose){System.err.println("After length adjust 2 to "+readlen+": dif="+dif[0]+"\n"+s+"\n");}
			
//			if(s.length()!=readlen){
//				assert(s.length()>readlen);
//				boolean start=randyCutPos.nextBoolean();
//				if(start){//Take first part of string
//					s=s.substring(0, readlen);
//				}else{//take last part of string
//					s=s.substring(s.length()-readlen);
//				}
//			}

			for(int j=0; j<SNPs; j++){s=addSNP(s, readlen, randy2);}
			
			for(int j=0; j<SUBs; j++){s=addSUB(s, 2, maxSubLen, readlen, randy2);}
//			for(int j=0; j<SUBs; j++){s=addSUB(s, 28, 28, readlen, randy2);}
			
			if(r.strand()==1){s=AminoAcid.reverseComplementBases(s);}
			r.bases=s.getBytes();
			
			final byte baseQuality;
			final byte slant;
			{
				byte baseSlant=(perfect ? (byte)5 : (byte)20);
//				slant=(byte)((randyQual.nextInt(baseSlant)+randyQual.nextInt(baseSlant)+1)/2);
				slant=(byte)Tools.min(randyQual.nextInt(baseSlant), randyQual.nextInt(baseSlant));
				if((randyQual.nextInt()&3)>0){
					int range=(perfect ? maxQualP-midQualP+1 : maxQual-midQual+1);
					int delta=Tools.min(randyQual.nextInt(range), randyQual.nextInt(range));
					baseQuality=(byte)((perfect ? midQualP : midQual)+delta);
				}else{
					int range=perfect ? midQualP-minQualP+1 : midQual-minQual+1;
					int delta=Tools.min(randyQual.nextInt(range), randyQual.nextInt(range));
					baseQuality=(byte)((perfect ? midQualP : midQual)-delta);
				}
			}
			
			if(USE_FIXED_QUALITY){
				r.quality=getFixedQualityRead(r.bases.length);
			}else{
				if(perfect){
					r.quality=QualityTools.makeQualityArray(
							r.bases.length, randyQual, minQualP, maxQualP, baseQuality, slant);
				}else{
					r.quality=QualityTools.makeQualityArray(
							r.bases.length, randyQual, minQual, maxQual, baseQuality, slant);
				}
			}
			for(int j=0; j<r.quality.length; j++){
				if(!AminoAcid.isFullyDefined(r.bases[j])){r.quality[j]=0;}
			}
			if(ADD_ERRORS_FROM_QUALITY && !perfect){addErrorsFromQuality(r, randyQual);}
//			System.err.println(Arrays.toString(r.quality));
			
			assert(r.bases.length==readlen);
			
			r.stop=r.start+readlen+dif[0]-1;
			
			assert(r.stop>r.start) : "\n"+Read.header()+"\n"+r+"\n"+SNPs+", "+SUBs+", "+INSs+", "+DELs+"\n"+s+"\n";
			
			if(colorspace){
				r=reads[i]=r.translateToColorspace(true);
				r.obj=new String(r.bases); //TODO - for testing
			}
			r.mapLength=r.bases.length;
			
			
			if(paired){

				Read r2=makeMate(r, mateLen, 
						maxSnps, maxInss, maxDels, maxSubs, 
						snpRate, insRate, delRate, subRate,
						maxInsertionLen, maxDeletionLen, maxSubLen,
						mateMiddleMin, mateMiddleMax, mateSameStrand,
						minQual, maxQual, baseQuality, slant, perfect);
				
				while(r2==null){
					r2=makeMate(r, mateLen, 
							maxSnps, maxInss, maxDels, maxSubs, 
							snpRate, insRate, delRate, subRate,
							maxInsertionLen, maxDeletionLen, maxSubLen,
							mateMiddleMin, mateMiddleMax, mateSameStrand,
							minQual, maxQual, baseQuality, slant, perfect);
				}
				
				r.mate=r2;
				r2.mate=r;
			}
//			System.err.println("Made "+r.start+" ~ "+r.stop+" = "+(r.stop-r.start));
		}
		
//		if(colorspace){
//			for(int i=0; i<reads.length; i++){
//				Read r=reads[i];
//				reads[i]=r.translateToColorspace(true);
//				reads[i].obj=new String(r.bases); //TODO - for testing
//			}
//		}
		
//		if(paired){
//			for(int i=0; i<reads.length; i++){
//				Read r=reads[i];
//				
//				Read r2=makeMate(r, mateLen, 
//						maxSnps, maxInss, maxDels, maxSubs, 
//						snpRate, insRate, delRate, subRate,
//						maxInsertionLen, maxDeletionLen, maxSubLen,
//						mateMiddleMin, mateMiddleMax, mateSameStrand,
//						baseQuality);
//				r.mate=r2;
//				r2.mate=r;
//			}
//		}
		
		return reads;
	}
	
	public Read makeMate(Read other, int readlen, 
			int maxSnps, int maxInss, int maxDels, int maxSubs, 
			float snpRate, float insRate, float delRate, float subRate,
			int maxInsertionLen, int maxDeletionLen, int maxSubLen,
			int minMiddle, int maxMiddle, boolean sameStrand, 
			int minQual, int maxQual, byte baseQuality, byte slant, boolean perfect){

		int x=-1, y=-1;
		int chrom=other.chrom;

		assert(maxMiddle>=minMiddle);
//		assert(minMiddle>=0);
		int midRange=maxMiddle-minMiddle+1;
		int middle=(randyMate.nextInt(midRange)+randyMate.nextInt(midRange))/2+minMiddle;
		byte strand=(byte) (sameStrand ? other.strand() : other.strand()^1);
		
//		System.out.println(sameStrand+": "+other.strand+" -> "+strand);
		
		if(other.strand()==Gene.PLUS){
			x=other.stop+middle;
		}else{
			x=other.start-middle-readlen;
		}
		y=x+readlen+(maxDeletionLen*maxDels);
		if(x<0){x=0; y=readlen-1; maxDels=0;}
		if(y>Data.getChromosome(chrom).maxIndex){y=Data.getChromosome(chrom).maxIndex; x=y-readlen+1; maxDels=0;}

		String s=Data.getChromosome(chrom).getString(x, y);
		
//		System.out.println("Making string length "+s.length()+" from "+x+"-"+y+" of "+Data.getChromosome(chrom).maxIndex);
		
		//I already do this later.
//		if(strand==Gene.MINUS){
//			s=AminoAcid.reverseComplementBases(s);
//		}
		
		long id=other.numericID;

		int SNPs=0;
		int INSs=0;
		int DELs=0;
		int SUBs=0;

//		assert(maxSnps==0 || (snpRate>.0001 && snpRate<=1)) : maxSnps+", "+snpRate;
//		assert(maxInss==0 || (snpRate>.0001 && insRate<=1)) : maxInss+", "+insRate;
//		assert(maxDels==0 || (snpRate>.0001 && delRate<=1)) : maxDels+", "+delRate;
//		assert(maxSubs==0 || (snpRate>.0001 && subRate<=1)) : maxSubs+", "+subRate;
		
		while(SNPs<maxSnps && randyMutationTypeMate.nextFloat()<snpRate){SNPs++;}
		while(INSs<maxInss && randyMutationTypeMate.nextFloat()<insRate){INSs++;}
		while(DELs<maxDels && randyMutationTypeMate.nextFloat()<delRate){DELs++;}// Note: SNPs are used instead of deletions.
		while(SUBs<maxSubs && randyMutationTypeMate.nextFloat()<subRate){SUBs++;}
		
		
		if(perfect){SNPs=INSs=DELs=SUBs=0;}

		int[] dif=new int[] {0, 0};

		//			int lengthDif=0;

		for(int j=0; j<DELs; j++){
			//				lengthDif-=s.length();
			s=addDeletion(s, 1, maxDeletionLen, readlen, dif, randy2Mate);
			//				lengthDif+=s.length();
		}
		if(s.length()>readlen){s=s.substring(0, readlen);}
		assert(s.length()==readlen);

		//			int preInsertLength=s.length();
		
		int insLen=0;
		for(int j=0; j<INSs; j++){

			int tempLimit=Tools.min(maxInsertionLen, readlen-4-insLen);
			if(tempLimit<1){break;}
			
			insLen-=s.length();
			s=addInsertion(s, 1, tempLimit, readlen, dif, randy2Mate);
			insLen+=s.length();
			assert(insLen<readlen);
		}
		//			int insertLength=s.length()-preInsertLength;

		if(s.length()!=readlen){
			assert(s.length()>readlen);
			s=s.substring(0, readlen);
		}

		//			if(s.length()!=readlen){
		//				assert(s.length()>readlen);
		//				boolean start=randyCutPos.nextBoolean();
		//				if(start){//Take first part of string
		//					s=s.substring(0, readlen);
		//				}else{//take last part of string
		//					s=s.substring(s.length()-readlen);
		//				}
		//			}

		for(int j=0; j<SNPs; j++){s=addSNP(s, readlen, randy2Mate);}
		for(int j=0; j<SUBs; j++){s=addSUB(s, 2, maxSubLen, readlen, randy2Mate);}

		Read r=new Read(s.getBytes(), chrom, strand, x, y, id, null, false);
		r.setSynthetic(true);
		r.setPairnum(1);
		
		assert(other!=null);
		assert(sameStrand == (r.strand()==other.strand())) : "\n"+r.toText(false)+"\n"+other.toText(false)+"\n\n"+
		sameStrand+", "+r.strand()+", "+other.strand()+"\n"+r.pairnum()+", "+other.pairnum()+"\n";

		if(r.strand()==Gene.MINUS){
			AminoAcid.reverseComplementBasesInPlace(r.bases);
//			r.start=r.stop-readlen-dif[0]+1;
		}else{
//			r.stop=r.start+readlen+dif[0]-1;
		}
		r.stop=r.start+readlen+dif[0]-1;

		if(USE_FIXED_QUALITY){
			r.quality=getFixedQualityRead(r.bases.length);
		}else{
			r.quality=QualityTools.makeQualityArray(r.bases.length, randyQualMate, minQual, maxQual, baseQuality, slant);
		}
		for(int j=0; j<r.quality.length; j++){
			if(!AminoAcid.isFullyDefined(r.bases[j])){r.quality[j]=0;}
		}
		if(ADD_ERRORS_FROM_QUALITY && !perfect){addErrorsFromQuality(r, randyQualMate);}

		assert(r.bases.length==readlen);

		r.stop=r.start+readlen+dif[0]-1;
		assert(r.stop>r.start) : "DELs="+DELs+", INSs="+INSs+", SUBs="+SUBs+", SNPs="+SNPs+
			", r.start="+r.start+", r.stop="+r.stop+", "+Data.getChromosome(chrom).maxIndex+"\n"+
			(other==null ? "" : "\n\n"+other.toText(false)+"\n\n");


		if(other.colorspace()){
			r=r.translateToColorspace(true);
			r.obj=new String(r.bases); //TODO - for testing
		}
		
		assert(sameStrand == (r.strand()==other.strand())) : "\n"+r.toText(false)+"\n"+other.toText(false)+"\n\n"+
			sameStrand+", "+r.strand()+", "+other.strand()+"\n"+r.pairnum()+", "+other.pairnum()+"\n";
		
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
	
	
//	public static int[] approxChromLengths=new int[] {
//		0,
//		600
//	};

//	public static final int[] randomChrom=fillRandomChrom(approxChromLengths);
	public static int[] randomChrom;
	
	private static int[] fillRandomChrom(){
		assert(Data.chromLengths!=null);
		int[] in=Arrays.copyOf(Data.chromLengths, Data.chromLengths.length);
		long total=Tools.sum(in);
		int div=(int)(total/1000);
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
//		return 345;
	}
	private static long seed=0;

	private final Random randy;
	private final Random randy2;
	private final Random randyMutationType;
	private final Random randyCSError;
	private final Random randyQual;
	
	private final Random randyMate;
	private final Random randy2Mate;
	private final Random randyMutationTypeMate;
	private final Random randyCSErrorMate;
	private final Random randyQualMate;

	private final Random randyPerfectRead;
	
	public final boolean paired;
	
	private long nextReadID=0;
	
	private static final byte[][] fixedQuality=new byte[1000][];
	
	public static final boolean USE_FIXED_QUALITY=false;
	public static final byte FIXED_QUALITY_VALUE=20;
	public static final boolean ADD_ERRORS_FROM_QUALITY=true;
	
	public static int mateLen=35;
	public static boolean mateSameStrand=false;
	public static int mateMiddleMin=-25; //default -25
	public static int mateMiddleMax=475; //default 475
	
	public static float PERFECT_READ_RATIO=0.8f;
	public static boolean verbose=false;
	
}
