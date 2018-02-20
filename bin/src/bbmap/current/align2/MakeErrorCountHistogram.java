package align2;

import java.util.ArrayList;

import stream.ConcurrentReadInputStream;
import stream.RTextInputStream;
import stream.Read;

import fileIO.ReadWrite;

public class MakeErrorCountHistogram {
	
	public static void main(String[] args){
		
		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
		
		long maxReads=0;
		RTextInputStream rtis=new RTextInputStream(fname1, fname2, maxReads);
		ConcurrentReadInputStream cris=new ConcurrentReadInputStream(rtis, maxReads);
		
		int[][] errors=process(cris);
		printHistogram(errors);
//		System.out.println("*** main() finished ***");
	}
	
	public static void printHistogram(int[][] errors){
		System.out.println("#Error Count Histogram: Number of Reads with X Mismatches");
		System.out.println("Errors\tRead 1\tRead 2\tPair\tRead 1 %\tRead 2 %\tPair %\t");
		
		long sum1=Tools.sum(errors[0]), sum2=Tools.sum(errors[1]), sum3=Tools.sum(errors[2]);
		
		for(int i=0; i<errors[0].length; i++){
			long e1=errors[0][i];
			long e2=errors[1][i];
			long e3=errors[2][i];
			double p1=e1*100d/(sum1);
			double p2=e2*100d/(sum2);
			double p3=e3*100d/(sum3);
			System.out.println(i+"\t"+e1+"\t"+e2+"\t"+e3+"\t"+
					String.format("%.4f", p1)+"\t"+String.format("%.4f", p2)+"\t"+String.format("%.4f", p3));
		}
	}
	
	public static int[][] process(ConcurrentReadInputStream cris){
		
		new Thread(cris).start();
		
		int[][] counts=new int[3][200];
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> readlist=ln.list;
		while(!readlist.isEmpty()){
			
			processList(readlist, counts);
			
			cris.returnList(ln, readlist.isEmpty());
			//System.err.println("Waiting on a list...");
			ln=cris.nextList();
			readlist=ln.list;
		}
		
		//System.err.println("Returning a list... (final)");
		assert(readlist.isEmpty());
		cris.returnList(ln, readlist.isEmpty());
		ReadWrite.closeStream(cris);
		
		return counts;
	}

	private static void processList(ArrayList<Read> list, int[][] counts) {
		for(Read r : list){
			processRead(r, counts);
		}
	}

	private static void processRead(Read r, int[][] counts) {

//		if(!r.paired()){return;}
		
//		if(r.containsIndels()){return;}
		
//		if(r.countMismatches()>6){return;}
//		if(r.avgQuality()<8){return;}
		
		
		int n1=-1, s1=-1, n2=-1, s2=-1, sum1=-1, sum2=-1;
		Read r2=r.mate;
		
		if(r.mapped() && r.valid() && r.match!=null){
			n1=count(r.match, 'N');
			s1=count(r.match, 'S');
			sum1=n1+s1;
		}
		if(r2!=null && r2.mapped() && r2.valid() && r2.match!=null){
			n2=count(r2.match, 'N');
			s2=count(r2.match, 'S');
			sum2=n2+s2;
		}

		if(sum1>-1){counts[0][sum1]++;}
		if(sum2>-1){counts[1][sum2]++;}
		if(sum1>-1 && sum2>-1){counts[2][sum1+sum2]++;}
		
	}
	


	public static int count(byte[] match, char symbol){
		assert(match!=null);
		int x=0;
		for(byte b : match){
			if(b==symbol){x++;}
		}
		return x;
	}
	
}
