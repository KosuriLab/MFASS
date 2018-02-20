package align2;

import java.util.ArrayList;
import java.util.Random;

import stream.Read;

public final class ListNum<K> {
	
	public ListNum(ArrayList<K> list_, long id_){
		list=list_;
		id=id_;
		if(GEN_RANDOM_NUMBERS && list!=null){
			for(K k : list){
				if(k!=null){
					((Read)k).rand=randy.nextDouble();
				}
			}
		}
	}

	public final ArrayList<K> list;
	public final long id;
	
	public static synchronized void setDeterministicRandom(boolean b){
		GEN_RANDOM_NUMBERS=b;
		if(b){
			randy=new Random(seed);
			seed++;
		}
	}
	public static boolean deterministicRandom(){
		return GEN_RANDOM_NUMBERS;
	}
	
	private static boolean GEN_RANDOM_NUMBERS=false;
	private static Random randy;
	private static long seed=0;
	
}
