package align2;

import java.util.Arrays;



public final class LongList{
	
	public LongList(){this(256);}
	
	public LongList(int initial){
		assert(initial>0);
		array=new long[initial];
	}
	
	public final void set(int loc, long value){
		if(loc>=array.length){
			resize((loc+1L)*2);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void increment(int loc, long value){
		if(loc>=array.length){
			resize((loc+1L)*2);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	public final long get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}
	
	public final void add(long x){
		if(size>=array.length){
			resize((size+1L)*2);
		}
		array[size]=x;
		size++;
	}
	
	public final void resize(long x){
		int size2=(int)min(x, Integer.MAX_VALUE);
		assert(size2>size);
		array=Arrays.copyOf(array, size2);
	}
	
	public final void shrink(){
		if(size==array.length){return;}
		array=Arrays.copyOf(array, size);
	}
	
	public final void shrinkToUnique(){
		//Assumes sorted.
		if(size<=0){
			shrink();
			return;
		}
		
		int unique=1;
		
		for(int i=1; i<size; i++){
			assert(array[i]>=array[i-1]);
			if(array[i]!=array[i-1]){unique++;}
		}
		if(unique==array.length){return;}
		long[] alt=new long[unique];
		
		alt[0]=array[0];
		for(int i=1, j=1; j<unique; i++){
			if(array[i]!=array[i-1]){
				alt[j]=array[i];
				j++;
			}
		}
		
		array=alt;
		size=alt.length;
	}
	
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
			if(array[i]!=0){
				sb.append(comma+"("+i+", "+array[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	private static final long min(long x, long y){return x<y ? x : y;}
	private static final long max(long x, long y){return x>y ? x : y;}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public long[] array;
	public int size=0;
	
}
