package align2;

import stream.SamLine;
import dna.Timer;
import fileIO.TextFile;
import fileIO.TextStreamWriter;

public class SplitSamFile {
	
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		t.start();
		TextFile tf=new TextFile(args[0], true, false);
		String out1=args[1];
		String out2=args[2];
		assert(!out1.equalsIgnoreCase(out2)) : "Output files are the same.";
		TextStreamWriter tsw1=new TextStreamWriter(out1, true, false, true);
		TextStreamWriter tsw2=new TextStreamWriter(out2, true, false, true);
		
		tsw1.start();
		tsw2.start();
		
		long plus=0;
		long minus=0;
		long other=0;
		
		String s=null;
		for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
			char c=s.charAt(0);
			if(c!='@'/* && c!=' ' && c!='\t'*/){
				SamLine sl=new SamLine(s);
				if(sl.mapped()){
					if(sl.strand()==0){
						tsw1.println(s);
						plus++;
					}else{
						tsw2.println(s);
						minus++;
					}
				}else{
					other++;
				}
			}
		}
		tf.close();
		tsw1.poison();
		tsw2.poison();
		
		System.err.println("Total reads:   \t"+(plus+minus+other));
		System.err.println("Plus reads:    \t"+(plus));
		System.err.println("Minus reads:   \t"+(minus));
		System.err.println("Unmapped reads:\t"+(other));
		
		tsw1.waitForFinish();
		tsw2.waitForFinish();
		t.stop();
		
		System.err.println("Time:          \t"+t);
		
	}
	
	
}
