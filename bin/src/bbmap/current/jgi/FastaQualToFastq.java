package jgi;

import dna.Timer;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Oct 25, 2012
 *
 */
public class FastaQualToFastq extends Thread {
	
	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		for(int i=0; i<args.length; i+=3){
			new FastaQualToFastq(args[i], args[i+1], args[i+2]).start();
		}
		t.stop();
//		System.out.println("Time: \t"+t);
	}
	
	final String in1, in2, out;
	
	public FastaQualToFastq(String in1_, String in2_, String out_){
		in1=in1_;
		in2=in2_;
		out=out_;
	}
	
	@Override
	public void run(){
		merge(in1, in2, out);
		System.out.println("Finished writing "+out);
	}

	/**
	 * @param string
	 * @param string2
	 */
	private void merge(String fasta, String qual, String out) {
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, true);
		TextFile tff=new TextFile(fasta, false, false);
		TextFile tfq=new TextFile(qual, false, false);
		tsw.start();

		String header="";
		StringBuilder sbf=new StringBuilder();
		boolean print=false;
		for(String s1=tff.nextLine(); s1!=null; s1=tff.nextLine()){
			if(s1.startsWith(">")){
				if(print){
					tsw.print("@");
					tsw.println(header);
					sbf.append('\n');
					tsw.print(sbf);
					sbf=new StringBuilder(sbf.length());

					String s2=tfq.nextLine();
					s2=tfq.nextLine();
					String[] qs=s2.toString().split(" ");

					StringBuilder sbq=new StringBuilder(qs.length+3);
					sbq.append('+').append('\n');
					for(int i=0; i<qs.length; i++){
						sbq.append((char)(Integer.parseInt(qs[i])+33));
					}
					sbq.append('\n');
					tsw.print(sbq);
				}
				print=true;
				header=s1.substring(1);
			}else{
				sbf.append(s1);
			}
		}

		if(sbf.length()>0){
//			System.out.println(header);
//			System.out.println(sbf);
			tsw.print("@");
			tsw.println(header);
			sbf.append('\n');
			tsw.print(sbf);
			sbf=new StringBuilder(sbf.length());
			
			String s2=tfq.nextLine();
			s2=tfq.nextLine();
//			System.out.println(s2);
			String[] qs=s2.toString().split(" ");
			StringBuilder sbq=new StringBuilder(qs.length+3);
			sbq.append('+').append('\n');
			for(int i=0; i<qs.length; i++){
				sbq.append((char)(Integer.parseInt(qs[i])+33));
			}
			sbq.append('\n');
			tsw.print(sbq);
		}
		
		tsw.poison();
		tff.close();
		tfq.close();
	}
}
