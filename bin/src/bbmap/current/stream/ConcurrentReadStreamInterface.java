package stream;

import align2.ListNum;

public interface ConcurrentReadStreamInterface extends Runnable{

	public ListNum<Read> nextList(); 

	public void returnList(ListNum<Read> list, boolean poison);

	public void run();
	
	public void shutdown();

	public void restart();
	public void close();
	
	/** Returns true for paired-end stream, false for single-end stream */
	public boolean paired();
	
	public Object[] producers();
	
	/** Return true if this stream has detected an error */
	public boolean errorState();
	
	public void setSampleRate(float rate, long seed);

}
