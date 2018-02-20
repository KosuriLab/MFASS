package align2;

import java.util.Comparator;

import stream.Read;

/**
 * Sorts longest reads first
 * @author Brian Bushnell
 * @date Jul 19, 2013
 *
 */
public final class ReadLengthComparator implements Comparator<Read> {
	
	private ReadLengthComparator(){}
	
	@Override
	public int compare(Read a, Read b) {
		int x=compare2(a, b);
		if(x==0){x=compare2(a.mate, b.mate);}
		if(x==0){x=a.id.compareTo(b.id);}
		if(x==0){x=a.numericID>b.numericID ? 1 : a.numericID<b.numericID ? -1 : 0;}
		return x;
	}

	private static int compare2(Read a, Read b) {
		if(a==b){return 0;}
		if(a==null){return 1;}
		if(b==null){return -1;}
		int x=compareByLength(a.bases, b.bases);
		return x;
	}
	
	private static int compareByLength(byte[] a, byte[] b){
		if(a==b){return 0;}
		if(a==null){return 1;}
		if(b==null){return -1;}
		return b.length-a.length;
//		if(a.length!=b.length){return b.length-a.length;}
//		for(int i=0; i<a.length; i++){
//			if(a[i]!=b[i]){return a[i]-b[i];}
//		}
//		return 0;
	}
	
	public static final ReadLengthComparator comparator=new ReadLengthComparator();
	
}
