package jgi;

import java.util.Arrays;
import java.util.Random;

/**
 * @author Brian Bushnell
 * @date Nov 2, 2012
 *
 */
public final class MatrixMult {
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		
		//Grab arguments
		int N=Integer.parseInt(args[0]);
		int threads=(args.length>1 ? Integer.parseInt(args[1]) : 1);
		int iters=(args.length>2 ? Integer.parseInt(args[2]) : 1);
		
		//Initialize arrays
		int[][] A=new int[N][N];
		int[][] B=new int[N][N];
		int[][] C=new int[N][N];
		
		Random randy=new Random(0);
		int max=100;
		
		//Fill random arrays
		for(int i=0; i<N; i++){
			for(int j=0; j<N; j++){
				A[i][j]=randy.nextInt(max);
				B[i][j]=randy.nextInt(max);
			}
		}
		
		//Start timer
		t.start();
		
		//Multiply matrices for some number of iterations
		for(int i=0; i<iters; i++){
			if(threads==1){
				//multiplyST(A, B, C);
				multiplySTT(A, B, C);
			}else{
				multiplyMTT(A, B, C, threads);
			}
		}
		
		t.stop();
		System.out.println("Time: \t"+t);
		
		if(verbose){printMatrices(A, B, C);}
		
	}
	
	
	/*   ~~~~~~~~~~~~~~~~~~~~~~   Singlethreaded Multiply   ~~~~~~~~~~~~~~~~~~~~~~  */
	
	private static void multiplyST(int[][] A, int[][] B, int[][] C) {
		for(int i=0; i<A.length; i++){
			int[] row=A[i];
			for(int j=0; j<B[0].length; j++){
				int sum=0;
				for(int k=0; k<row.length; k++){
					sum=sum+row[k]*B[k][j];
				}
				C[i][j]=sum;
			}
		}
	}
	
	
	/*   ~~~~~~~~~~~~~~~~~~~~~~   Singlethreaded Transposed Multiply   ~~~~~~~~~~~~~~~~~~~~~~  */
	
	private static void multiplySTT(int[][] A, int[][] B, int[][] C) {
		int[][] Bt=makeTranspose(B);
		for(int i=0; i<A.length; i++){
			int[] row=A[i];
			for(int j=0; j<Bt.length; j++){
				int[] col=Bt[j];
				int sum=0;
				for(int k=0; k<row.length; k++){
					sum=sum+row[k]*col[k];
				}
				C[i][j]=sum;
			}
		}
	}
	
	
	/*   ~~~~~~~~~~~~~~~~~~~~~~   Multithreaded Transposed Multiply Setup   ~~~~~~~~~~~~~~~~~~~~~~  */
	
	private static void multiplyMTT(int[][] A, int[][] B, int[][] C, int threads) {
		int[][] Bt=makeTranspose(B);
		
		//Decide work unit size
		int chunk=(A.length+threads-1)/threads;
		assert(chunk*threads>=A.length);
		threads=(A.length+chunk-1)/chunk;
		
		if(!printed){System.out.println("Using "+threads+" threads, chunksize "+chunk);}
		printed=true;
		
		//Create workers
		MultThread[] workers=new MultThread[threads];
		for(int i=0; i<threads; i++){
			workers[i]=new MultThread(A, Bt, C, chunk*i, min(chunk*(i+1), A.length));
			workers[i].start();
		}
		
		//Wait for workers to finish
		for(MultThread mt : workers){
			while(mt.getState()!=Thread.State.TERMINATED){
				try{mt.join();}
				catch (InterruptedException e) {e.printStackTrace();}
			}
		}
	}
	
	
	/*   ~~~~~~~~~~~~~~~~~~~~~~   Multithreaded Transposed Multiply Worker   ~~~~~~~~~~~~~~~~~~~~~~  */
	
	private static class MultThread extends Thread{
		
		/** Create a new worker thread that will multiply A and B to fill rows (min to max-1) in output matrix C */
		MultThread(int[][] A_, int[][] Bt_, int[][] C_, int min_, int max_){
			A=A_;
			Bt=Bt_;
			C=C_;
			min=min_;
			max=max_;
		}
		
		//Thread's run method
		@Override
		public void run(){
			for(int i=min; i<max; i++){
				int[] row=A[i];
				for(int j=0; j<Bt.length; j++){
					int[] col=Bt[j];
					int sum=0;
					for(int k=0; k<row.length; k++){
						sum=sum+row[k]*col[k];
					}
					C[i][j]=sum;
				}
			}
		}
		
		//Local fields for the individual thread
		final int[][] A, Bt, C;
		final int min, max;
	}
	
	

	
	/*   ~~~~~~~~~~~~~~~~~~~~~~   Helper Methods   ~~~~~~~~~~~~~~~~~~~~~~  */

	/** Creates a new matrix that is the transpose of the input matrix */
	private static int[][] makeTranspose(int[][] in){
		int[][] out=new int[in[0].length][in.length];
		for(int i=0; i<in.length; i++){
			for(int j=0; j<in[0].length; j++){
				out[j][i]=in[i][j];
			}
		}
		return out;
	}
	
	private static void printMatrices(int[][] A, int[][] B, int[][] C){
		for(int i=0; i<C.length; i++){
			System.out.println(Arrays.toString(A[i]));
		}
		System.out.println();
		for(int i=0; i<C.length; i++){
			System.out.println(Arrays.toString(B[i]));
		}
		System.out.println();
		for(int i=0; i<C.length; i++){
			System.out.println(Arrays.toString(C[i]));
		}
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	private static class Timer {
		
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
	

	
	/*   ~~~~~~~~~~~~~~~~~~~~~~   Fields   ~~~~~~~~~~~~~~~~~~~~~~  */
	
	private static boolean printed=false;
	private static boolean verbose=false;
	
}
