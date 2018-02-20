package stream;

import java.util.ArrayList;

import align2.RandomReads;
import align2.Shared;
import align2.Tools;

public class RandomReadInputStream extends ReadInputStream {
	
	public RandomReadInputStream(long number_, int readlen_, 
			int maxSnps_, int maxInss_, int maxDels_, int maxSubs_,
			float snpRate_, float insRate_, float delRate_, float subRate_,
			int maxInsertionLen_, int maxDeletionLen_,  int maxSubLen_,
			int minChrom_, int maxChrom_, boolean colorspace_, boolean paired_,
			int minQual_, int midQual_, int maxQual_){
		
		number=number_;
		readlen=readlen_;

		maxInsertionLen=maxInsertionLen_;
		maxSubLen=maxSubLen_;
		maxDeletionLen=maxDeletionLen_;
		
		minChrom=minChrom_;
		maxChrom=maxChrom_;
		
		maxSnps=maxSnps_;
		maxInss=maxInss_;
		maxDels=maxDels_;
		maxSubs=maxSubs_;

		snpRate=snpRate_;
		insRate=insRate_;
		delRate=delRate_;
		subRate=subRate_;
		
		colorspace=colorspace_;
		paired=paired_;
		
		minQual=(byte) minQual_;
		midQual=(byte) midQual_;
		maxQual=(byte) maxQual_;
		
		restart();
	}
	
	@Override
	public void start() {}
	
	
	@Override
	public boolean hasMore() {
		return number>consumed;
	}

	@Override
	public Read next() {
		if(consumed>=number){return null;}
		if(buffer==null || next>=buffer.length){fillBuffer();}
		Read r=buffer[next];
		buffer[next]=null;
		next++;
		consumed++;
		return r;
	}

	@Override
	public synchronized Read[] nextBlock() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(consumed>=number){return null;}
		if(buffer==null || next>=buffer.length){fillBuffer();}
		Read[] r=buffer;
		buffer=null;
		if(r!=null && r.length==0){r=null;}
		consumed+=(r==null ? 0 : r.length);
		return r;
	}
	
	@Override
	public synchronized ArrayList<Read> nextList() {
		return toList(nextBlock());
	}
	public final boolean preferArrays(){return true;}
	
	private synchronized void fillBuffer(){
		buffer=null;
		next=0;
		
		long toMake=number-generated;
		if(toMake<1){return;}
		toMake=Tools.min(toMake, BUF_LEN);
		
		Read[] reads=rr.makeRandomReadsX((int)toMake, readlen, 
				maxSnps, maxInss, maxDels, maxSubs,
				snpRate, insRate, delRate, subRate,
				maxInsertionLen, maxDeletionLen, maxSubLen,
				minChrom, maxChrom, colorspace,
				minQual, midQual, maxQual);
		
		generated+=reads.length;
		assert(generated<=number);
		buffer=reads;
	}
	
	public synchronized void restart(){
		next=0;
		buffer=null;
		consumed=0;
		generated=0;
		rr=new RandomReads(1, paired);
		rr.mateLen=readlen;
	}

	@Override
	public boolean close() {return false;}

	@Override
	public boolean paired() {
		return paired;
	}
	
	private Read[] buffer=null;
	private int next=0;
	
	public static final int BUF_LEN=Shared.READ_BUFFER_LENGTH;

	public long generated=0;
	public long consumed=0;
	
	public long number=100000;
	public int readlen=50;

	public int maxInsertionLen=6;
	public int maxSubLen=6;
	public int maxDeletionLen=100;
	
	public int minChrom=1;
	public int maxChrom=22;
	
	public int maxSnps=2;
	public int maxInss=2;
	public int maxDels=2;
	public int maxSubs=2;

	public float snpRate=0.3f;
	public float insRate=0.15f;
	public float delRate=0.15f;
	public float subRate=0.10f;
	
	public final boolean colorspace;
	public final boolean paired;

	public final byte minQual;
	public final byte midQual;
	public final byte maxQual;
	
	private RandomReads rr;

}
