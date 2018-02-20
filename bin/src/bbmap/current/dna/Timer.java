package dna;

public class Timer {
	
	public Timer(){}
	
	public long start(){
		time1=time2=System.nanoTime();
		elapsed=0;
		return time1;
	}
	
	public long stop(){
		time2=System.nanoTime();
		elapsed=time2-time1;
		return time2;
	}
	
	public String toString(){
		return String.format("%.3f seconds.", elapsed/1000000000d);
	}

	public long time1;
	public long time2;
	/** in nanos */
	public long elapsed;
	
}
