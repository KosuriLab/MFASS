package align2;

import java.util.Comparator;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date May 30, 2013
 *
 */
public final class ReadErrorComparator implements Comparator<Read>{
	
	@Override
	public int compare(Read r1, Read r2) {

		int a=(r1.errors+(r1.mate==null ? 0 : r1.mate.errors));
		int b=(r2.errors+(r2.mate==null ? 0 : r2.mate.errors));
		if(a!=b){return a-b;}
		
		a=(r1.bases.length+(r1.mate==null ? 0 : r1.mate.bases.length));
		b=(r2.bases.length+(r2.mate==null ? 0 : r2.mate.bases.length));
		if(a!=b){return b-a;}
		
		float a2=(r1.expectedErrors()+(r1.mate==null ? 0 : r1.mate.expectedErrors()));
		float b2=(r2.expectedErrors()+(r2.mate==null ? 0 : r2.mate.expectedErrors()));
		if(a2!=b2){return a2>b2 ? 1 : -1;}
		
		if(r1.numericID<r2.numericID){return -1;}
		else if(r1.numericID>r2.numericID){return 1;}
		
		if(!r1.id.equals(r2.id)){return r1.id.compareTo(r2.id);}
		return 0;
	}
	
	public static final ReadErrorComparator comparator=new ReadErrorComparator();
	
}
