package align2;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;

import dna.ChromosomeArray;
import dna.Data;
import dna.FastaToChromArrays;
import fileIO.ReadWrite;
import fileIO.SummaryFile;

/**
 * @author Brian Bushnell
 * @date Sep 25, 2013
 *
 */
public class RefToIndex {

	public static void makeIndex(String reference, int build, PrintStream sysout, int keylen){
		assert(reference!=null);
		{
			File f=new File(reference);
			if(!f.exists() || !f.isFile() || !f.canRead()){
				if(!reference.startsWith("stdin")){
					throw new RuntimeException("Cannot read file "+f.getAbsolutePath());
				}
			}
		}

		String s=IndexMaker4.fname(1, 1, keylen, 1, false);
		String dir=new File(s).getParent();
		dir=dir.replace('\\', '/');
		final String base=dir.substring(0, dir.length()-7);
		final String args=(Shared.COMMAND_LINE==null ? "null" : Arrays.toString(Shared.COMMAND_LINE));
		final String indexlog=base+"build"+build+"_"+
				(System.nanoTime()&Long.MAX_VALUE)+"."+((args==null ? (reference==null ? "null" : reference) : args).hashCode()&Integer.MAX_VALUE)+".log";
		dir=dir.replace("ref/index/", "ref/genome/");
		String sf=dir+"/summary.txt";
		if(!NODISK && new File(sf).exists() && SummaryFile.compare(sf, reference)){
			//do nothing
			if(LOG && !NODISK){
				if(!new File(base).exists()){new File(base).mkdirs();}
				ReadWrite.writeString(new Date()+"\nFound an already-written genome for build "+build+".\n"+args+"\n", indexlog, true);
			}
			sysout.println("NOTE:\tIgnoring reference file because it already appears to have been processed.");
			sysout.println("NOTE:\tIf you wish to regenerate the index, please manually delete "+dir+"/summary.txt");
		}else{
			if(NODISK){}
			else{//Delete old data if present
				File f=new File(dir);
				if(f.exists()){
					File[] f2=f.listFiles();
					if(f2!=null && f2.length>0){
						if(OVERWRITE || f2[0].getAbsolutePath().equals(new File(reference).getAbsolutePath())){
							sysout.println("NOTE:\tDeleting contents of "+dir+" because reference is specified and overwrite="+OVERWRITE);
							if(LOG && !NODISK){ReadWrite.writeString(new Date()+"\nDeleting genome for build "+build+".\n"+args+"\n", indexlog, true);}
							for(File f3 : f2){
								if(f3.isFile()){
									String f3n=f3.getName();
									if((f3n.contains(".chrom") || f3n.endsWith(".txt") || f3n.endsWith(".txt.gz")) && !f3n.endsWith("list.txt")){
										f3.delete();
									}
								}
							}
						}else{
							sysout.println(Arrays.toString(f2));
							if(LOG && !NODISK){ReadWrite.writeString(new Date()+"\nFailed to overwrite genome for build "+build+".\n"+args+"\n", indexlog, true);}
							throw new RuntimeException("\nThere is already a reference at location '"+f.getAbsolutePath()+"'.  " +
									"Please delete it (and the associated index), or use a different build ID, " +
									"or remove the 'reference=' parameter from the command line, or set overwrite=true.");
						}
					}
				}

				dir=dir.replace("ref/genome/", "ref/index/");
				f=new File(dir);
				if(f.exists()){
					File[] f2=f.listFiles();
					if(f2!=null && f2.length>0){
						if(OVERWRITE){
							sysout.println("NOTE:\tDeleting contents of "+dir+" because reference is specified and overwrite="+OVERWRITE);
							if(LOG && !NODISK){ReadWrite.writeString(new Date()+"\nDeleting index for build "+build+".\n"+args+"\n", indexlog, true);}
							for(File f3 : f2){
								if(f3.isFile()){f3.delete();}
							}
						}else{
							if(LOG && !NODISK){ReadWrite.writeString(new Date()+"\nFailed to overwrite index for build "+build+".\n"+args+"\n", indexlog, true);}
							throw new RuntimeException("\nThere is already an index at location '"+f.getAbsolutePath()+"'.  " +
									"Please delete it, or use a different build ID, or remove the 'reference=' parameter from the command line.");
						}
					}
				}
			}

			if(!NODISK){
				sysout.println("Writing reference.");
				if(LOG && !NODISK){
					if(!new File(base).exists()){new File(base).mkdirs();}
					ReadWrite.writeString(new Date()+"\nWriting genome for build "+build+".\n"+args+"\n", indexlog, true);
				}
			}

			int oldzl=ReadWrite.ZIPLEVEL;
			ReadWrite.ZIPLEVEL=Tools.max(4, ReadWrite.ZIPLEVEL);

			//assert(false) : "minScaf="+minScaf+", midPad="+midPad+", maxChromLen="+maxChromLen+
			//		", startPad="+startPad+", stopPad="+stopPad+", FastaToChromArrays.END_PADDING="+FastaToChromArrays.END_PADDING;
			
			maxChromLen=maxChromLen>0 ? maxChromLen : AUTO_CHROMBITS ? FastaToChromArrays.MAX_LENGTH : ((1L<<(31-(chrombits<0 ? 2 : chrombits)))-200000);
			minScaf=minScaf>-1 ? minScaf : FastaToChromArrays.MIN_SCAFFOLD;
			midPad=midPad>-1 ? midPad : FastaToChromArrays.MID_PADDING;
			startPad=startPad>-1 ? startPad : FastaToChromArrays.START_PADDING;
			stopPad=stopPad>-1 ? stopPad : FastaToChromArrays.END_PADDING;
			
			String[] ftcaArgs=new String[] {reference, ""+build, "writeinthread=false", "genscaffoldinfo="+genScaffoldInfo, "retain", "waitforwriting=false",
					"gzip="+(Data.CHROMGZ), "chromc="+Data.CHROMC, "maxlen="+maxChromLen,
					"writechroms="+(!NODISK), "minscaf="+minScaf, "midpad="+midPad, "startpad="+startPad, "stoppad="+stopPad, "nodisk="+NODISK};

			chromlist=FastaToChromArrays.main2(ftcaArgs);

			ReadWrite.ZIPLEVEL=oldzl;
		}

	}

	public static boolean AUTO_CHROMBITS=true;
	public static boolean LOG=false;
	public static boolean NODISK=false;
	public static boolean OVERWRITE=true;
	public static boolean genScaffoldInfo=true;
	
	public static long maxChromLen=-1;
	
	public static int minScaf=-1, midPad=-1, stopPad=-1, startPad=-1;
	public static int chrombits=-1;
//	public static int minScaf=FastaToChromArrays.MIN_SCAFFOLD;
//	public static int midPad=FastaToChromArrays.MID_PADDING;
//	public static int startPad=FastaToChromArrays.START_PADDING;
//	public static int stopPad=FastaToChromArrays.END_PADDING;
	
	public static ArrayList<ChromosomeArray> chromlist=null;
	
}
