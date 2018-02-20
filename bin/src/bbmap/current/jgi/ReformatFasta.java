package jgi;

import align2.Tools;
import dna.Timer;
import fileIO.TextFile;
import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Nov 1, 2012
 *
 */
public class ReformatFasta {
	
	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		
		String in=args[0];
		String out=args[1];
		int minLen=0;
		int maxLen=Integer.MAX_VALUE;
		boolean rename=false;
		if(args.length>2){
			minLen=Integer.parseInt(args[2]);
		}
		if(args.length>3){
			maxLen=Integer.parseInt(args[3]);
		}
		if(args.length>4){
			rename=Tools.parseBoolean(args[4]);
		}
		
		if(in.equalsIgnoreCase(out)){throw new RuntimeException("in == out");}
		
		TextFile tf=new TextFile(in, false, false);
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		tsw.start();
		
		long kept=0;
		long dropped=0;
		
		String header=null;
		StringBuilder sb=new StringBuilder(100);
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.charAt(0)=='>'){
				if(header!=null){
					if(sb.length()>=minLen){
						if(rename){
							tsw.println(">"+(kept+1));
						}else{
							tsw.println(header);
						}
//						tsw.println(sb);
						printAsLines(sb, tsw, maxLen);
						kept++;
					}else{
						dropped++;
					}
				}
				header=s;
				sb=new StringBuilder(100);
			}else{
				sb.append(s);
			}
		}
		if(header!=null){
			if(sb.length()>=minLen){
				if(rename){
					tsw.println(">"+(kept+1));
				}else{
					tsw.println(header);
				}
//				tsw.println(sb);
				printAsLines(sb, tsw, maxLen);
				kept++;
			}else{
				dropped++;
			}
		}
		tsw.poison();
		t.stop();
		System.out.println("Time: \t"+t);
		System.out.println("Kept "+kept+" reads.");
		System.out.println("Dropped "+dropped+" reads.");
		
	}
	
	private static void printAsLines(final CharSequence sb, final TextStreamWriter tsw, final int max){
		
		final int len=sb.length();
		
		if(len<=max){
			tsw.println(sb);
		}else{
			for(int i=0; i<len; i+=max){
				int b=min(len, i+max);
				CharSequence s=sb.subSequence(i, b);
				tsw.println(s);
			}
		}
		
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	
}
