package kmer;

import java.util.ArrayList;

/**
 * @author Brian Bushnell
 * @date Oct 22, 2013
 *
 */
public class KmerNode {

	public KmerNode(long pivot_){
		pivot=pivot_;
	}
	
	public KmerNode(long pivot_, int value_){
		pivot=pivot_;
		count=value_;
	}
	
	int increment(long kmer){
		if(pivot<0){pivot=kmer; return (count=1);} //Allows initializing empty nodes to -1
		if(kmer<pivot){
			if(left==null){left=new KmerNode(kmer, 1); return 1;}
			return left.increment(kmer);
		}else if(kmer>pivot){
			if(right==null){right=new KmerNode(kmer, 1); return 1;}
			return right.increment(kmer);
		}else{
			if(count<Integer.MAX_VALUE){count++;}
			return count;
		}
	}
	
	/** Returns number of nodes added */
	int set(long kmer, int value){
		if(pivot<0){pivot=kmer; count=value; return 1;} //Allows initializing empty nodes to -1
		if(kmer<pivot){
			if(left==null){left=new KmerNode(kmer, value); return 1;}
			return left.set(kmer, value);
		}else if(kmer>pivot){
			if(right==null){right=new KmerNode(kmer, value); return 1;}
			return right.set(kmer, value);
		}else{
			count=value;
		}
		return 0;
	}
	
	/** Returns number of nodes added */
	int setIfNotPresent(long kmer, int value){
		if(pivot<0){pivot=kmer; count=value; return 1;} //Allows initializing empty nodes to -1
		if(kmer<pivot){
			if(left==null){left=new KmerNode(kmer, value); return 1;}
			return left.setIfNotPresent(kmer, value);
		}else if(kmer>pivot){
			if(right==null){right=new KmerNode(kmer, value); return 1;}
			return right.setIfNotPresent(kmer, value);
		}
		return 0;
	}
	
	KmerNode get(long kmer){
		if(kmer<pivot){
			return left==null ? null : left.get(kmer);
		}else if(kmer>pivot){
			return right==null ? null : right.get(kmer);
		}else{
			return this;
		}
	}
	
	KmerNode getNodeOrParent(long kmer){
		if(pivot==kmer || pivot<0){return this;}
		if(kmer<pivot){return left==null ? this : left.getNodeOrParent(kmer);}
		return right==null ? this : right.getNodeOrParent(kmer);
	}
	
	boolean insert(KmerNode n){
		assert(pivot!=-1);
		if(n.pivot<pivot){
			if(left==null){left=n; return true;}
			return left.insert(n);
		}else if(n.pivot>pivot){
			if(right==null){right=n; return true;}
			return right.insert(n);
		}else{
			return false;
		}
	}
	
	int getCount(long kmer){
//		KmerNode node=get(kmer);
//		return node==null ? 0 : node.count;
		
//		if(kmer<pivot){
//			return left==null ? 0 : left.getCount(kmer);
//		}else if(kmer>pivot){
//			return right==null ? 0 : right.getCount(kmer);
//		}else{
//			return count;
//		}
		
//		if(kmer==pivot){
//			return count;
//		}else if(kmer<pivot){
//			return left==null ? 0 : left.getCount(kmer);
//		}else{
//			assert(kmer>pivot);
//			return right==null ? 0 : right.getCount(kmer);
//		}
		
		KmerNode n=this;
		while(n!=null && n.pivot!=kmer){
			n=(kmer<n.pivot ? n.left : n.right);
		}
		return n==null ? 0 : n.count;
		
	}
	
	boolean contains(long kmer){
		KmerNode node=get(kmer);
		return node!=null;
	}
	
	void traversePrefix(ArrayList<KmerNode> list){
		if(left!=null){left.traversePrefix(list);}
		list.add(this);
		if(right!=null){right.traversePrefix(list);}
	}
	
	void traverseInfix(ArrayList<KmerNode> list){
		list.add(this);
		if(left!=null){left.traverseInfix(list);}
		if(right!=null){right.traverseInfix(list);}
	}
	
	KmerNode rebalance(ArrayList<KmerNode> list){
		assert(list.isEmpty());
		traversePrefix(list);
		KmerNode n=this;
		if(list.size()>2){
			n=rebalance(list, 0, list.size()-1);
		}
		list.clear();
		return n;
	}
	
	public StringBuilder dumpKmersAsText(StringBuilder sb, int k){
		if(count<1){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		sb.append(AbstractKmerTable.toText(pivot, count, k).append('\n'));
		if(left!=null){left.dumpKmersAsText(sb, k);}
		if(right!=null){left.dumpKmersAsText(sb, k);}
		return sb;
	}
	
	private static KmerNode rebalance(ArrayList<KmerNode> list, int a, int b){
		final int size=b-a+1;
		final int middle=a+size/2;
		final KmerNode n=list.get(middle);
		if(size<4){
			if(size==1){
				n.left=n.right=null;
			}else if(size==2){
				KmerNode n1=list.get(a);
				n.left=n1;
				n.right=null;
				n1.left=n1.right=null;
			}else{
				assert(size==3);
				KmerNode n1=list.get(a), n2=list.get(b);
				n.left=n1;
				n.right=n2;
				n1.left=n1.right=null;
				n2.left=n2.right=null;
			}
		}else{
			n.left=rebalance(list, a, middle-1);
			n.right=rebalance(list, middle+1, b);
		}
		return n;
	}
	
	long pivot;
	int count;
	KmerNode left, right;
	
}
