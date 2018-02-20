package kmer;

import java.util.ArrayList;


import fileIO.TextStreamWriter;

import align2.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 23, 2013
 *
 */
public final class HashForest extends AbstractKmerTable {
	
	public HashForest(int initialSize, boolean autoResize_){
		if(initialSize>1){
			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
		}else{
			initialSize=1;
		}
		prime=initialSize;
		sizeLimit=(long) (initialSize*resizeMult);
		array=new KmerNode[prime];
		autoResize=autoResize_;
	}
	
	int increment(long kmer){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=(kmer<n.pivot ? n.left : n.right);
		}
		if(n==null){
			n=new KmerNode(kmer, 1);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				if(kmer<prev.pivot){
					prev.left=n;
				}else{
					prev.right=n;
				}
			}
			if(autoResize && size>sizeLimit){resize();}
		}else{
			n.count++;
			if(n.count<0){n.count=Integer.MAX_VALUE;}
		}
		return n.count;
	}
	
	int incrementAndReturnNumCreated(long kmer){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=(kmer<n.pivot ? n.left : n.right);
		}
		if(n==null){
			n=new KmerNode(kmer, 1);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				if(kmer<prev.pivot){
					prev.left=n;
				}else{
					prev.right=n;
				}
			}
			if(autoResize && size>sizeLimit){resize();}
			return 1;
		}else{
			n.count++;
			if(n.count<0){n.count=Integer.MAX_VALUE;}
			return 0;
		}
	}
	
	int set(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerNode n=array[cell];
		if(n==null){
			array[cell]=new KmerNode(kmer, value);
		}else{
			x=n.set(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	public int setIfNotPresent(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerNode n=array[cell];
		if(n==null){
			array[cell]=new KmerNode(kmer, value);
		}else{
			x=n.setIfNotPresent(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	final KmerNode get(long kmer){
//		int cell=(int)(kmer%prime);
//		KmerNode n=array[cell];
//		return n==null ? null : n.get(kmer);
		
		int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		while(n!=null && n.pivot!=kmer){
			n=(kmer<n.pivot ? n.left : n.right);
		}
		return n;
	}
	
	public final int getCount(long kmer){
//		KmerNode node=get(kmer);
//		return node==null ? 0 : node.count;
		
		int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		while(n!=null && n.pivot!=kmer){
			n=(kmer<n.pivot ? n.left : n.right);
		}
		return n==null ? 0 : n.count;
	}
	
//	int getCount(LongM kmer){
//		KmerNode node=get(kmer.value());
//		return node==null ? 0 : node.count;
//	}
	
	public boolean contains(long kmer){
		return get(kmer)!=null;
	}
	
	boolean insert(KmerNode n){
		n.left=null;
		n.right=null;
		int cell=(int)(n.pivot%prime);
		if(array[cell]==null){
			array[cell]=n;
			return true;
		}
		return array[cell].insert(n);
	}
	
	public void rebalance(){
		ArrayList<KmerNode> list=new ArrayList<KmerNode>(1000);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i]=array[i].rebalance(list);}
		}
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
		sizeLimit=Tools.max((long)(size*1.4), (long)(maxLoadFactor*prime));

		final long maxAllowedByLoadFactor=(long)(size*minLoadMult);
		final long minAllowedByLoadFactor=(long)(size*maxLoadMult);
		assert(maxAllowedByLoadFactor>=minAllowedByLoadFactor);
		if(maxAllowedByLoadFactor<prime){return;}
		
		long x=10+(long)(prime*resizeMult);
		x=Tools.max(x, minAllowedByLoadFactor);
		x=Tools.min(x, maxAllowedByLoadFactor);
		
		int prime2=(int)Tools.min(maxPrime, Primes.primeAtLeast(x));
		
		if(prime2<=prime){return;}
		
		prime=prime2;
//		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
		KmerNode[] old=array;
		array=new KmerNode[prime2];
		ArrayList<KmerNode> list=new ArrayList<KmerNode>(1000);
		for(int i=0; i<old.length; i++){
			if(old[i]!=null){
				old[i].traverseInfix(list);
				for(KmerNode n : list){insert(n);}
				list.clear();
			}
		}
		sizeLimit=Tools.max((long)(size*1.4), (long)(maxLoadFactor*prime));
	}
	
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k){
		tsw.print("HashForest:\n");
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null && node.count>0){
				StringBuilder sb=new StringBuilder();
				tsw.print(node.dumpKmersAsText(sb, k));
			}
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
	
	KmerNode[] array;
	int prime;
	long size=0;
	long sizeLimit;
	final boolean autoResize;
	
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE);
	final static float resizeMult=2.5f; //Resize by a minimum of this much
	final static float minLoadFactor=0.75f; //Resize by enough to get the load above this factor
	final static float maxLoadFactor=2.5f; //Resize by enough to get the load under this factor
	final static float minLoadMult=1/minLoadFactor;
	final static float maxLoadMult=1/maxLoadFactor;
	

	
}
