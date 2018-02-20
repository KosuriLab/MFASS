package stream;

import java.util.ArrayList;

public abstract class ReadInputStream {
	

	public abstract Read next();
	
//	public abstract Read[] nextBlock();
	
	public Read[] nextBlock(){
		throw new RuntimeException("Deprecated, use nextList() instead.");
	}
	
//	public ArrayList<Read> nextList(){
//		throw new RuntimeException("Not supported.");
//	}
	
	public abstract ArrayList<Read> nextList();
	
	public abstract boolean hasMore();

	public abstract void restart();
	
	/** Returns true if there was an error, false otherwise */
	public abstract boolean close();

	public abstract boolean paired();

	protected static final ArrayList<Read> toList(Read[] array){
		if(array==null || array.length==0){return null;}
		ArrayList<Read> list=new ArrayList<Read>(array.length);
		for(int i=0; i<array.length; i++){list.add(array[i]);}
		return list;
	}
	
	/** Return true if this stream has detected an error */
	public boolean errorState(){return errorState;}
	/** TODO */
	protected boolean errorState=false;

	public abstract boolean preferArrays();
	public final boolean preferLists(){return !preferArrays();}
	public final boolean preferBlocks(){return preferArrays();}

	public abstract void start();
	
}
