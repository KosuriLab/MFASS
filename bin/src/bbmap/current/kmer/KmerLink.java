package kmer;

import java.util.ArrayList;

/**
 * @author Brian Bushnell
 * @date Oct 22, 2013
 *
 */
public class KmerLink {

	public KmerLink(long pivot_){
		pivot=pivot_;
	}
	
	public KmerLink(long pivot_, int value_){
		pivot=pivot_;
		count=value_;
	}
	
	int increment(long kmer){
		if(pivot<0){pivot=kmer; return (count=1);} //Allows initializing empty nodes to -1
		if(kmer==pivot){
			if(count<Integer.MAX_VALUE){count++;}
			return count;
		}
		if(next==null){next=new KmerLink(kmer, 1); return 1;}
		return next.increment(kmer);
	}
	
	/** Returns number of nodes added */
	int set(long kmer, int value){
		if(pivot<0){pivot=kmer; count=value; return 1;} //Allows initializing empty nodes to -1
		if(kmer==pivot){count=value; return 0;}
		if(next==null){next=new KmerLink(kmer, value); return 1;}
		return next.set(kmer, value);
	}
	
	/** Returns number of nodes added */
	int setIfNotPresent(long kmer, int value){
		if(pivot<0){pivot=kmer; count=value; return 1;} //Allows initializing empty nodes to -1
		if(kmer==pivot){return 0;}
		if(next==null){next=new KmerLink(kmer, value); return 1;}
		return next.setIfNotPresent(kmer, value);
	}
	
	KmerLink get(long kmer){
		if(kmer==pivot){return this;}
		return next==null ? null : next.get(kmer);
	}
	
	boolean insert(KmerLink n){
		assert(pivot!=-1);
		if(pivot==n.pivot){return false;}
		if(next==null){next=n; return true;}
		return next.insert(n);
	}
	
	int getCount(long kmer){
		KmerLink node=get(kmer);
		return node==null ? 0 : node.count;
	}
	
	boolean contains(long kmer){
		KmerLink node=get(kmer);
		return node!=null;
	}
	
	void traversePrefix(ArrayList<KmerLink> list){
		if(next!=null){next.traversePrefix(list);}
		list.add(this);
	}
	
	void traverseInfix(ArrayList<KmerLink> list){
		list.add(this);
		if(next!=null){next.traverseInfix(list);}
	}
	
	KmerLink rebalance(ArrayList<KmerLink> list){
		throw new RuntimeException("Unsupported.");
	}
	
	private static KmerLink rebalance(ArrayList<KmerLink> list, int a, int b){
		throw new RuntimeException("Unsupported.");
	}
	
	long pivot;
	int count;
	KmerLink next;
	
}
