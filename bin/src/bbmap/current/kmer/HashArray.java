package kmer;

import java.util.ArrayList;
import java.util.Arrays;


import fileIO.TextStreamWriter;

import align2.Tools;

/**
 * Stores kmers in a long[] and counts in an int[], with a victim cache.
 * @author Brian Bushnell
 * @date Oct 25, 2013
 *
 */
public final class HashArray extends AbstractKmerTable {
	
	public HashArray(int initialSize, boolean autoResize_){
		if(initialSize>1){
			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
		}else{
			initialSize=1;
		}
		prime=initialSize;
		sizeLimit=(long)(sizeLimit=(long)(maxLoadFactor*prime));
		array=new long[prime+extra];
		counts=new int[prime+extra];
		victims=new HashForest(Tools.max(10, initialSize/8), autoResize_);
		Arrays.fill(array, -1);
		autoResize=autoResize_;
	}
	


	int increment(final long kmer){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				counts[cell]++;
				if(counts[cell]<0){counts[cell]=Integer.MAX_VALUE;}
				return counts[cell];
			}else if(n==-1){
				array[cell]=kmer;
				size++;
				counts[cell]=1;
				if(autoResize && size>sizeLimit){resize();}
				return 1;
			}
		}
		return victims.increment(kmer);
	}
	
	int incrementAndReturnNumCreated(final long kmer){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				counts[cell]++;
				if(counts[cell]<0){counts[cell]=Integer.MAX_VALUE;}
				return 0;
			}else if(n==-1){
				array[cell]=kmer;
				size++;
				counts[cell]=1;
				if(autoResize && size>sizeLimit){resize();}
				return 1;
			}
		}
		return victims.incrementAndReturnNumCreated(kmer);
	}
	
	
	int set(long kmer, int value){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				counts[cell]=value;
				return 0;
			}else if(n==-1){
				array[cell]=kmer;
				counts[cell]=value;
				size++;
				if(autoResize && size>sizeLimit){resize();}
				return 1;
			}
		}
		return victims.set(kmer, value);
	}
	
	public int setIfNotPresent(long kmer, int value){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				return 0;
			}else if(n==-1){
				array[cell]=kmer;
				counts[cell]=value;
				size++;
				if(autoResize && size>sizeLimit){resize();}
				return 1;
			}
		}
//		System.err.println("size="+size+", prime="+prime+", limit="+sizeLimit);
		return victims.setIfNotPresent(kmer, value);
	}
	
	final Object get(long kmer){
		throw new RuntimeException("Unimplemented");
	}
	
	public final int getCount(long kmer){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				return counts[cell];
			}else if(n==-1){
				return 0;
			}
		}
		return victims.getCount(kmer);
	}
	
	public boolean contains(long kmer){
		int cell=(int)(kmer%prime);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				return true;
			}else if(n==-1){
				return false;
			}
		}
		return victims.contains(kmer);
	}
	
//	boolean insert(KmerNode n){
//		n.left=null;
//		n.right=null;
//		int cell=(int)(n.pivot%prime);
//		if(array[cell]==null){
//			array[cell]=n;
//			return true;
//		}
//		return array[cell].insert(n);
//	}
	
	public void rebalance(){
//		ArrayList<KmerNode> list=new ArrayList<KmerNode>(1000);
//		for(int i=0; i<array.length; i++){
//			if(array[i]!=null){array[i]=array[i].rebalance(list);}
//		}
		throw new RuntimeException("Unimplemented");
	}
	
//	boolean containsKey(LongM key){
//		return contains(key.value());
//	}
	
//	boolean put(LongM key, int number){
//		set(key.value(), number);
//		return true;
//	}
	
	synchronized void resize(){
//		assert(false);
//		System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
		if(prime>=maxPrime){
			sizeLimit=0xFFFFFFFFFFFFL;
			return;
		}

		final long maxAllowedByLoadFactor=(long)(size*minLoadMult);
		final long minAllowedByLoadFactor=(long)(size*maxLoadMult);

//		sizeLimit=Tools.min((long)(maxLoadFactor*prime), maxPrime);
		
		assert(maxAllowedByLoadFactor>=minAllowedByLoadFactor);
		if(maxAllowedByLoadFactor<prime){
			sizeLimit=(long)(maxLoadFactor*prime);
			return;
		}
		
		long x=10+(long)(prime*resizeMult);
		x=Tools.max(x, minAllowedByLoadFactor);
		x=Tools.min(x, maxAllowedByLoadFactor);
		
		int prime2=(int)Tools.min(maxPrime, Primes.primeAtLeast(x));
		
		if(prime2<=prime){
			sizeLimit=(long)(maxLoadFactor*prime);
			assert(prime2==prime) : "Resizing to smaller array? "+size+", "+prime+", "+x;
			return;
		}
		
		final long oldSize=size, oldVSize=victims.size;
		
		prime=prime2;
//		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
		long[] oldk=array;
		int[] oldc=counts;
		KmerNode[] oldv=victims.array;
		array=new long[prime2+extra];
		Arrays.fill(array, -1);
		counts=new int[prime2+extra];
		ArrayList<KmerNode> list=new ArrayList<KmerNode>((int)(victims.size)); //Can fail if more than Integer.MAX_VALUE
		for(int i=0; i<oldv.length; i++){
			if(oldv[i]!=null){oldv[i].traverseInfix(list);}
		}
		Arrays.fill(oldv, null);
		victims.size=0;
		size=0;
		sizeLimit=Long.MAX_VALUE;
		
		for(int i=0; i<oldk.length; i++){
			if(oldk[i]>-1){set(oldk[i], oldc[i]);}
		}

		for(KmerNode n : list){
			if(n.pivot>-1){set(n.pivot, n.count);}
		}
		
		assert(oldSize+oldVSize==size+victims.size) : oldSize+", "+oldVSize+" -> "+size+", "+victims.size;
		
		sizeLimit=(long)(maxLoadFactor*prime);
	}
	
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k){
		tsw.print("HashArray:\n");
		for(int i=0; i<array.length; i++){
			long kmer=array[i];
			if(kmer!=-1){
				tsw.print(toText(kmer, counts[i], k).append('\n'));
				tsw.println(""+kmer);
			}
		}
		if(victims!=null){
			victims.dumpKmersAsText(tsw, k);
		}
		return true;
	}
	
	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#size()
	 */
	@Override
	public long size() {return size;}
	
	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#arrayLength()
	 */
	@Override
	public int arrayLength() {return array.length;}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#canResize()
	 */
	@Override
	boolean canResize() {return true;}

	/* (non-Javadoc)
	 * @see jgi.AbstractKmerTable#canRebalance()
	 */
	@Override
	public boolean canRebalance() {return true;}
	
	long[] array;
	int[] counts;
	HashForest victims;
	int prime;
	long size=0;
	long sizeLimit;
	final boolean autoResize;
	
	final static int extra=19;
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE-extra);
	final static float resizeMult=2f; //Resize by a minimum of this much
	final static float minLoadFactor=0.6f; //Resize by enough to get the load above this factor
	final static float maxLoadFactor=0.95f; //Resize by enough to get the load under this factor
	final static float minLoadMult=1/minLoadFactor;
	final static float maxLoadMult=1/maxLoadFactor;
	

	
}
