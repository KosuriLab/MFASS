package jgi;

import java.util.ArrayList;

/**
 * @author Brian Bushnell
 * @date Apr 17, 2013
 *
 */
public class AssemblyStatsWrapper {
	
	public static void main(String[] args){
		ArrayList<String> alist=new ArrayList<String>();
		ArrayList<String> ilist=new ArrayList<String>();
		
		alist.add("");
		alist.add("header=t");
		alist.add("showspeed=f");
		alist.add("addname=t");
		alist.add("k=0");
		
		for(String arg : args){
			if(!arg.contains("=")){
				ilist.add("in="+arg);
			}else{
				String[] split=arg.split("=");
				if(split[0].equalsIgnoreCase("in") || split[0].equalsIgnoreCase("ref")){
					ilist.add(arg);
				}else{
					alist.add(arg);
				}
			}
		}
		
		
		String[] args2=alist.toArray(new String[alist.size()]);
		for(int i=0; i<ilist.size(); i++){
			String s=ilist.get(i);
			args2[0]=s;
			if(i>0){
				args2[1]="header=f";
				AssemblyStats2.reset();
				System.gc();
				synchronized(AssemblyStatsWrapper.class){
					try {
						AssemblyStatsWrapper.class.wait(100);
					} catch (InterruptedException e) {}
				}
				Thread.yield();
			}
			AssemblyStats2.main(args2);
		}
		
	}
	
}
