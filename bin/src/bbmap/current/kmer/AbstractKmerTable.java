package kmer;

import dna.AminoAcid;
import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Oct 23, 2013
 *
 */
public abstract class AbstractKmerTable {
	
	/** Returns count */
	abstract int increment(long kmer);
	
	/** Returns number of entries created */
	abstract int incrementAndReturnNumCreated(final long kmer);

	abstract int set(long kmer, int value);
	
	/** Returns number of kmers added */
	public abstract int setIfNotPresent(long kmer, int value);

	abstract Object get(long kmer);

	public abstract int getCount(long kmer);

	public abstract boolean contains(long kmer);

//	abstract boolean insert(KmerLink n);

	public abstract void rebalance();

	abstract void resize();

	public abstract long size();
	public abstract int arrayLength();
	abstract boolean canResize();
	public abstract boolean canRebalance();
	
	public abstract boolean dumpKmersAsText(TextStreamWriter tsw, int k);

	static final StringBuilder toText(long kmer, int count, int k){
		StringBuilder sb=new StringBuilder(k+10);
		return toText(kmer, count, k, sb);
	}

//	static final StringBuilder toText(long kmer, int count, int k, StringBuilder sb){
//		for(int i=0; i<k; i++){
//			int x=(int)(kmer&3);
//			sb.append((char)AminoAcid.numberToBase[x]);
//			kmer>>=2;
//		}
//		sb.reverse();
//		sb.append('\t');
//		sb.append(count);
//		return sb;
//	}
	
	static final StringBuilder toText(long kmer, int count, int k, StringBuilder sb){
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(2*i))&3);
			sb.append((char)AminoAcid.numberToBase[x]);
		}
		sb.append('\t');
		sb.append(count);
		return sb;
	}
	
	static void appendKmerText(long kmer, int count, int k, StringBuilder sb){
		sb.setLength(0);
		toText(kmer, count, k, sb);
		sb.append('\n');
	}
	
	
	/** For buffered tables. */
	long flush(){
		throw new RuntimeException("Not supported.");
	}

}
