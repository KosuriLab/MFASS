package align2;

import java.util.Arrays;



public final class IntList{
	
	public IntList(){this(256);}
	
	public IntList(int initial){
		assert(initial>0);
		array=new int[initial];
	}
	
	public final void set(int loc, int value){
		if(loc>=array.length){
			resize((loc+1)*2);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void increment(int loc, int value){
		if(loc>=array.length){
			resize((loc+1)*2);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	public final int get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}
	
	
	
	public final void add(int x){
		if(size>=array.length){
			resize(max(size*2, 1));
		}
		array[size]=x;
		size++;
	}
	
	public final void resize(int size2){
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
		int[] alt=new int[unique];
		
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
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public int[] array;
	public int size=0;
	
}
