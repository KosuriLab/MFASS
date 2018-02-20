package kmer;

import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Nov 22, 2013
 *
 */
public class HashBuffer extends AbstractKmerTable {
	
	public HashBuffer(AbstractKmerTable[] tables_, int buflen_){
		tables=tables_;
		buflen=buflen_;
		assert(buflen<Short.MAX_VALUE);
		ways=tables.length;
		sizes=new short[ways];
		buffers=new long[ways][buflen];
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#increment(long)
	 */
	@Override
	int increment(long kmer) {
		throw new RuntimeException("Not supported");
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#increment(long)
	 */
	@Override
	public
	int incrementAndReturnNumCreated(long kmer) {
		final int way=(int)(kmer%ways);
		long[] buffer=buffers[way];
		buffer[sizes[way]]=kmer;
		sizes[way]++;
		if(sizes[way]>=buflen){
			return incrementBuffer(way);
		}
		return 0;
	}
	
	private int incrementBuffer(final int way){
		final int size=sizes[way];
		if(size<1){return 0;}
		sizes[way]=0;
		final long[] buffer=buffers[way];
		int added=0;
		final AbstractKmerTable table=tables[way];
		synchronized(table){
			for(int i=0; i<size; i++){
				final long kmer=buffer[i];
				added+=table.incrementAndReturnNumCreated(kmer);
			}
		}
		return added;
	}
	
	public final long flush(){
		long added=0;
		for(int i=0; i<ways; i++){added+=incrementBuffer(i);}
		return added;
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#set(long, int)
	 */
	@Override
	int set(long kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#setIfNotPresent(long, int)
	 */
	@Override
	public int setIfNotPresent(long kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#get(long)
	 */
	@Override
	Object get(long kmer) {
		final int way=(int)(kmer%ways);
		return tables[way].get(kmer);
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#getCount(long)
	 */
	@Override
	public int getCount(long kmer) {
		final int way=(int)(kmer%ways);
		return tables[way].getCount(kmer);
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#contains(long)
	 */
	@Override
	public boolean contains(long kmer) {
		final int way=(int)(kmer%ways);
		return tables[way].contains(kmer);
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#rebalance()
	 */
	@Override
	public void rebalance() {
		throw new RuntimeException("Unimplemented method");
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#resize()
	 */
	@Override
	void resize() {
		throw new RuntimeException("Unimplemented method");
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#size()
	 */
	@Override
	public long size() {
		throw new RuntimeException("Unimplemented method");
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#arrayLength()
	 */
	@Override
	public int arrayLength() {
		throw new RuntimeException("Unimplemented method");
	}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#canResize()
	 */
	@Override
	final boolean canResize() {return false;}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#canRebalance()
	 */
	@Override
	public final boolean canRebalance() {return false;}
	
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k){
		throw new RuntimeException("Unimplemented method");
	}
	
	private final AbstractKmerTable[] tables;
	private final int buflen;
	private final int ways;	
	private final short[] sizes;
	private final long[][] buffers;

}
