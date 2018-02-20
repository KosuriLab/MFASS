package jgi;

import java.util.ArrayList;

import stream.SamLine;

import fileIO.TextFile;
import fileIO.TextStreamWriter;

/**
 * Add s
 * @author Brian Bushnell
 * @date Jul 3, 2013
 *
 */
public class AddStopTags {
	
	public static void main(String[] args){
		String in=args[0];
		String out=args[1];
		TextStreamWriter tsw=new TextStreamWriter(out, false, false, true);
		tsw.start();
		
		TextFile tf=new TextFile(in, true, false);
		String line;
		for(line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.charAt(0)=='@'){
				tsw.println(line);
			}else{
				SamLine sl=new SamLine(line);
				ArrayList<String> list=sl.optional;
				if(list==null){sl.optional=list=new ArrayList<String>(1);}
				
				for(int i=0; i<list.size(); i++){
					String s=list.get(i);
					if(s.startsWith("YI:f:PID=")){
						s=s.replace("YI:f:PID=", "YI:f:");
						list.set(i, s);
						break;
					}
				}
				
				if(sl.mapped()){
					list.add(sl.makeStopTag(sl.pos, sl.seq.length, sl.cigar, false));
				}
				tsw.println(sl.toString());
			}
		}
		tf.close();
		tsw.poisonAndWait();
	}
	
}
