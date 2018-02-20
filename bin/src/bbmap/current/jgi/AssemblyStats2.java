package jgi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import align2.IntList;
import align2.LongList;
import align2.Tools;

import dna.FastaToChromArrays;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;


/**
 * @author Brian Bushnell
 * @date Sep 16, 2013
 *
 */
public final class AssemblyStats2 {
	
	/*
	 
	fasta_stats2.linux -n <number of N between contigs> contigs.fa
	e.g.

	fasta_stats2.linux -n 0 contigs.fa # for aplg assemblies
	fasta_stats2.linux -n 10 contigs.fa # for velvet
	  
	  
	 
	 Main genome scaffold total: 1610
	 Main genome contig total:   7844
	 Main genome scaffold sequence total: 726.6 MB
	 Main genome contig sequence total:   689.4 MB (->  5.1% gap)
	 Main genome scaffold N/L50: 6/62.2 MB
	 Main genome contig N/L50:   331/429.0 KB
	 Number of scaffolds > 50 KB: 122
	 % main genome in scaffolds > 50 KB: 98.9%

	  Minimum    Number    Number     Total        Total     Scaffold
	 Scaffold      of        of      Scaffold      Contig     Contig
	  Length   Scaffolds  Contigs     Length       Length    Coverage
	 --------  ---------  -------  -----------  -----------  --------
	     All     1,610      7,844  726,616,606  689,442,341    94.88%
	    1 kb     1,610      7,844  726,616,606  689,442,341    94.88%
	  2.5 kb     1,468      7,677  726,334,758  689,171,164    94.88%
	    5 kb       537      6,496  723,058,922  685,949,825    94.87%
	   10 kb       321      6,176  721,557,480  684,511,419    94.87%
	   25 kb       138      5,900  718,873,396  681,879,275    94.85%
	   50 kb       122      5,854  718,322,923  681,420,273    94.86%
	  100 kb        83      5,660  715,543,850  679,452,337    94.96%
	  250 kb        47      5,326  709,779,897  675,162,461    95.12%
	  500 kb        32      5,073  704,645,704  671,472,605    95.29%
	    1 mb        19      4,735  695,996,631  664,862,860    95.53%
	  2.5 mb        15      4,587  689,883,367  659,102,480    95.54%
	    5 mb        13      4,463  681,669,379  651,024,951    95.50%
	*/
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		t.start();
		if(args.length==0 || (args.length==1 && 
				(args[0].equalsIgnoreCase("-h") || args[0].equals("-help") || args[0].equals("--help") || args[0].equals("-?") || args[0].equals("?")))){
			System.out.println("\nUsage: java -Xmx120m jgi.AssemblyStats3 <input file>");
			System.out.println("\nOptional flags:");
			System.out.println("in=<file>    \tThe 'in=' flag is only needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.");
			System.out.println("format=1     \tUses variable units like MB and KB, and is designed for compatibility with existing tools.");
			System.out.println("format=2     \tUses only whole numbers of bases, with no commas in numbers, and is designed for machine parsing.");
			System.out.println("format=3     \tOutputs stats in 2 rows of tab-delimited columns: a header row and a data row.");
			System.out.println("format=4     \tLike 3 but with scaffold data only.");
			System.out.println("format=5     \tLike 3 but with contig data only.");
			System.out.println("format=6     \tLike 3 but the header starts with a #.");
			System.out.println("gc=<file>    \tPrint gc statistics per scaffold to a file (or stdout).");
			System.out.println("gcformat=1   \tid start stop A C G T N GC");
			System.out.println("gcformat=2   \tid gc");
			System.out.println("gcformat=3   \tid length A C G T N GC");
			System.out.println("gcformat=4   \tid length gc");
			System.out.println("gchist=<file>\tPrint gc content histogram to this file.");
			System.out.println("gcbins=200   \tNumber of bins in gc histogram.");
			System.out.println("n=10         \tMinimum number of consecutive Ns between contigs.");
			System.out.println("k=13         \tDisplay BBMap's estimated memory usage for this genome with specified kmer length.");
			System.out.println("showspeed=t  \tSet to 'f' to suppress display of processing speed.");
			System.out.println("minscaf=0    \tIgnore scaffolds shorter than this.");
			System.out.println("n_=t         \tThis flag will prefix the terms 'contigs' and 'scaffolds' with 'n_' in formats 3-6.");
//			System.out.println("verbose=t    \tSet to false to remove superfluous info.");
//			System.out.println("Output is always tab-delimited.  AGCT are fractions of defined bases; N is fraction of total bases.");
			System.exit(0);
		}
		
		boolean benchmark=false;
		ReadWrite.USE_UNPIGZ=true;
		
		String in=null, out=null, gc=null, gchist=null, scaffoldHistFile=null;
		int maxNs=-1;
		
		int gchistdecimals=-1;
		
		for(int i=0; i<args.length; i++){

			if(true){
				final String arg=args[i];
				String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if("null".equalsIgnoreCase(b)){b=null;}
				if(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>0)){a=a.substring(1);}

				if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
					//jvm argument; do nothing
				}else if(arg.contains("=") && (a.equals("in") || a.equals("ref"))){
					in=b;
				}else if(a.equals("gc")){
					gc=b;
					if(b==null || "summaryonly".equalsIgnoreCase(b) || "none".equalsIgnoreCase(b)){
						gc=null;
					}
				}else if(a.equals("gchist")){
					gchist=b;
					if(b==null || "none".equalsIgnoreCase(b)){
						gchist=null;
					}
				}else if(a.equals("gchistdecimals")){
					gchistdecimals=Integer.parseInt(b);
				}else if(a.equals("gcbins")){
					int x=Integer.parseInt(b);
					if(x>0){GCBINS=Integer.parseInt(b);}
				}else if(a.equals("shist") || a.equals("scaffoldhist")){
					scaffoldHistFile=b;
					if(b==null || "none".equalsIgnoreCase(b)){
						scaffoldHistFile=null;
					}
				}else if(a.equals("out")){
					out=b;
					if(b==null || "summaryonly".equalsIgnoreCase(b) || "none".equalsIgnoreCase(b)){
						out=null;
					}else if("benchmark".equalsIgnoreCase(b)){
						benchmark=true;
						out=null;
						gc=null;
					}
				}else if(a.equals("benchmark")){
					benchmark=Tools.parseBoolean(b);
					if(benchmark){
						out=null;
						gc=null;
					}
				}else if(a.equals("format")){
					FORMAT=Integer.parseInt(b);
					if(FORMAT<1 || FORMAT>6){
						throw new RuntimeException("\nUnknown format: "+FORMAT+"; valid values are 1 and 2.\n");
					}
				}else if(a.equals("gcformat")){
					GCFORMAT=Integer.parseInt(b);
					if(GCFORMAT!=1 && GCFORMAT!=2 && GCFORMAT!=3 && GCFORMAT!=4){
						throw new RuntimeException("\nUnknown gcformat: "+GCFORMAT+"; valid values are 1 and 2.\n");
					}
				}else if(a.equals("cutoff")){
					cutoff=Long.parseLong(b);
				}else if(a.equals("k") || a.equals("bbmapkmer")){
					bbmapkmer=Integer.parseInt(b);
				}else if(a.equals("overwrite") || a.equals("ow")){
					overwrite=Tools.parseBoolean(b);
				}else if(a.equals("n_")){
					N_UNDERSCORE=Tools.parseBoolean(b);
				}else if(a.equals("bf2")){
					ByteFile.FORCE_MODE_BF1=!(ByteFile.FORCE_MODE_BF2=Tools.parseBoolean(b));
				}else if(a.equals("usegzip") || a.equals("gzip")){
					ReadWrite.USE_GZIP=Tools.parseBoolean(b);
				}else if(a.equals("usepigz") || a.equals("pigz")){
					if(b!=null && Character.isDigit(b.charAt(0))){
						int zt=Integer.parseInt(b);
						if(zt<1){ReadWrite.USE_PIGZ=false;}
						else{
							ReadWrite.USE_PIGZ=true;
							if(zt>1){
								ReadWrite.MAX_ZIP_THREADS=zt;
								ReadWrite.ZIP_THREAD_DIVISOR=1;
							}
						}
					}else{ReadWrite.USE_PIGZ=Tools.parseBoolean(b);}
				}else if(a.equals("usegunzip") || a.equals("gunzip")){
					ReadWrite.USE_GUNZIP=Tools.parseBoolean(b);
				}else if(a.equals("useunpigz") || a.equals("unpigz")){
					ReadWrite.USE_UNPIGZ=Tools.parseBoolean(b);
				}else if(a.equals("ziplevel") || a.equals("zl")){
					ReadWrite.ZIPLEVEL=Integer.parseInt(b);
				}else if(a.equals("header") || a.equals("useheader")){
					useheader=Tools.parseBoolean(b);
				}else if(a.equals("addfilename") || a.equals("addname")){
					addfilename=Tools.parseBoolean(b);
				}else if(a.equals("minscaf")){
					MINSCAF=Integer.parseInt(b);
				}else if(a.equals("showspeed") || a.equals("ss")){
					showspeed=Tools.parseBoolean(b);
				}else if(a.equals("printheadersize") || a.equals("phs")){
					printheadersize=Tools.parseBoolean(b);
				}else if(a.equals("skipduplicatelines") || a.equals("sdl")){
					skipDuplicateLines=Tools.parseBoolean(b);
				}else if(a.equals("printduplicatelines") || a.equals("pdl")){
					skipDuplicateLines=!Tools.parseBoolean(b);
				}else if(a.equals("showbbmap")){
					if(!Tools.parseBoolean(b)){bbmapkmer=0;}
				}else if(a.equals("contigbreak") || (arg.contains("=") && (a.equals("n") || a.equals("-n")))){
					maxNs=Integer.parseInt(b);
				}else if(i>0 && (a.equals("n") || a.equals("-n")) && b!=null){
					maxNs=Integer.parseInt(b);
				}else if(a.equals("-n") && b==null && args.length>i+1){
					maxNs=Integer.parseInt(args[i+1]);
				}else if(in==null && i==0 && !args[i].contains("=")){
					in=args[i];
				}else if(maxNs<0 && i==1 && !args[i].contains("=") && Character.isDigit(args[i].charAt(0))){
					maxNs=Integer.parseInt(args[i]);
				}else if(gc==null && i==2 && !args[i].contains("=")){
					gc=args[i];
				}else{
					throw new RuntimeException("Unknown parameter "+args[i]);
				}
			}
		}
		
		if(gchistdecimals<1){
			if(GCBINS==2 || GCBINS==5 || GCBINS==10){
				GCHIST_DECIMALS=1;
			}else if(GCBINS==20 || GCBINS==25 || GCBINS==50 || GCBINS==100){
				GCHIST_DECIMALS=2;
			}
			if(GCBINS>1000 && GCHIST_DECIMALS<4){GCHIST_DECIMALS=4;}
			if(GCBINS>10000 && GCHIST_DECIMALS<5){GCHIST_DECIMALS=5;}
		}else{
			GCHIST_DECIMALS=gchistdecimals;
		}
		
		if(maxNs<0){maxNs=10;}
		
		long[] counts=null;
		long sum=0;

		if(out==null || out.equalsIgnoreCase("stdout") || out.equalsIgnoreCase("standardout")){out=null;}
		
		InputStream is=null;
		{
			if(in==null){throw new RuntimeException("No input file.");}
			if(in.equalsIgnoreCase("stdin") || in.equalsIgnoreCase("standardin")){
				is=System.in;
			}else{
				File f=new File(in);
				if((!f.exists() || f.isDirectory()) && !in.toLowerCase().startsWith("stdin")){
					throw new RuntimeException("Input file does not appear to be valid: "+in);
				}
			}
		}
		
		
		if(is==null){is=ReadWrite.getInputStream(in, false, true);}
		try {
			if(benchmark){sum=bench(is);}
			else{counts=countFasta(is, gc, maxNs);}
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			if(is!=System.in){is.close();}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if(tlistS!=null && tlistS.size()>0){Collections.sort(tlistS);}
		if(llistS!=null && llistS.size>0){Arrays.sort(llistS.array, 0, llistS.size);}
		
		t.stop();
		
		if(benchmark){
			printBenchResults(t, counts, sum, in);
		}else{
//			System.err.println("\nclistS="+clistS+"\nslistS="+slistS+"\nsclist1S="+sclist1S+"\nsclist2S="+sclist2S+"\nllistS="+llistS+"\ntlistS="+tlistS+"\n"); //***
//			System.err.println("\nclistS.size="+clistS.size+"\nslistS.size="+slistS.size+"\nsclist1S.size="+sclist1S.size+"\nsclist2S.size="+sclist2S.size+"\nllistS.size="+llistS.size+"\ntlistS.size()="+tlistS.size()+"\n"); //***
			printResults(t, counts, sum, in, clistS, slistS, sclist1S, sclist2S, llistS, tlistS, out);
			
			writeHistFile(scaffoldHistFile, slistS, tlistS, false);
			
			if(gchist!=null){printGCHist(gchist);}
		}
		
	}
	
	
	private static void printGCHist(String gchist){
		if(!Tools.canWrite(gchist, overwrite)){
			System.err.println("Can't write gc histogram because file exists and overwrite="+overwrite);
			assert(false);
		}else{
			long gchistsum=Tools.sum(gchistS);
			long gchistsumbb=Tools.sum(gchist_bb_S);
			double invsum=(gchistsum==0 ? 0 : 1.0/gchistsum);
			double invsumbb=(gchistsum==0 ? 0 : 1.0/gchistsumbb);
			double invbins=1.0/(GCBINS==0 ? 1 : GCBINS);
//			assert(false) : Arrays.toString(gchistS);
			StringBuilder sb=new StringBuilder();
			sb.append(String.format("#GC\tscaffolds\tfraction\tlength\tlen_fraction\n"));
			for(int i=0; i<GCBINS; i++){
				sb.append(String.format("%."+GCHIST_DECIMALS+"f\t%d\t%.5f\t%d\t%.5f\n", 
						i*invbins, gchistS[i], gchistS[i]*invsum, gchist_bb_S[i], gchist_bb_S[i]*invsumbb));
			}
			if(gchist.equalsIgnoreCase("stdout")){
				System.out.println(sb);
			}else{
				ReadWrite.writeString(sb, gchist);
			}
		}
	}
	
	
	public static void printBenchResults(Timer t, long[] counts, long sum, String in){
		System.err.println("Time: \t"+t);
		long bytes=new File(in).length();
		if(bytes<1){bytes=LIMSUM;}
		double mbps1=bytes*1000d/t.elapsed;
		double mbps2=sum*1000d/t.elapsed;
		System.err.println(String.format("Raw Speed:         \t%.2f MBytes/s",mbps1));
		System.err.println(String.format("Uncompressed Speed:\t%.2f MBytes/s",mbps2));
	}
	
	
	public static void printResults(Timer t, long[] counts, long sum, String in, LongList clist, LongList slist, LongList sclist1, LongList sclist2, 
			LongList llist, ArrayList<Triple> tlist, String out){
		
		String name=in;
		if(in!=null && !in.toLowerCase().startsWith("stdin")){
			try {
				File f=new File(in);
				name=f.getCanonicalPath();
			} catch (IOException e) {}
		}
		
		long contigs=0;
		long scaffolds=0;
		long contiglen=0;
		long scaflen=0;
		long contigs1;
		long contiglen2;
		long maxScaf=0, maxContig=0;
		long[] carray=clist.array;
		long[] sarray=slist.array;
		long[] scarray1=sclist1.array;
		long[] scarray2=sclist2.array;

		long[] larray=llist.array;
		
		StringBuilder sb=new StringBuilder(), sb2=new StringBuilder();
		
		for(int i=0; i<carray.length; i++){
			long x=carray[i];
			if(x>0){
				contigs+=x;
				contiglen+=(x*i);
				maxContig=i;
			}
		}
		
		for(int i=0; i<sarray.length; i++){
			long x=sarray[i];
			if(x>0){
				scaffolds+=x;
				scaflen+=(x*i);
				maxScaf=i;
			}
		}
		
		contigs+=llist.size;
		for(int i=0; i<llist.size; i++){
			long x=larray[i];
			assert(x>0);
			contiglen+=x;
			maxContig=Tools.max(maxContig, x);
		}
		
		scaffolds+=tlist.size();
		for(Triple tp : tlist){
			scaflen+=tp.length;
			maxScaf=Tools.max(maxScaf, tp.length);
		}
		
		if(FORMAT<3){
			sb.append("Main genome scaffold total:         \t"+scaffolds+"\n");
			sb.append("Main genome contig total:           \t"+contigs+"\n");
		}
		
		if(FORMAT==1){
			sb.append("Main genome scaffold sequence total:\t"+String.format("%.3f MB",scaflen/1000000f)+"\n");
			sb.append("Main genome contig sequence total:  \t"+String.format("%.3f MB  \t%.3f%% gap",contiglen/1000000f,(scaflen-contiglen)*100f/scaflen)+"\n");
		}else if(FORMAT==2){
			sb.append("Main genome scaffold sequence total:\t"+scaflen+"\n");
			sb.append("Main genome contig sequence total:  \t"+String.format("%d  \t%.3f%% gap",contiglen,(scaflen-contiglen)*100f/scaflen)+"\n");
		}else if(FORMAT==3 || FORMAT==6){
			
		}else if(FORMAT==4){
			
		}else if(FORMAT==5){
			
		}else{throw new RuntimeException("Unknown format");}
		
		if(FORMAT<3){
			sb2.append("\nMinimum \tNumber        \tNumber        \tTotal         \tTotal         \tScaffold\n");
			sb2.append("Scaffold\tof            \tof            \tScaffold      \tContig        \tContig  \n");
			sb2.append("Length  \tScaffolds     \tContigs       \tLength        \tLength        \tCoverage\n");
			sb2.append("--------\t--------------\t--------------\t--------------\t--------------\t--------\n");
		}
		
		int[] lims=new int[] {0, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 
				250000, 500000, 1000000, 2500000, 5000000, 10000000, 25000000, 50000000, 100000000, 250000000};
		
		int lidx=0;
		int next=0;
		long csum=contigs;
		long ssum=scaffolds;
		long clen=contiglen;
		long slen=scaflen;
		
		long cn50=-1;
		long ln50=-1;
		long cl50=-1;
		long ll50=-1;
		
		long shalf=slen/2;
		long chalf=clen/2;
		
		long numOver50=0;
		float fractionOver50=0;
		
		
		//Disable printing of 50~500 when not needed
//		{
//			boolean b=true;
//			for(int i=0; i<500 && i<sarray.length && b; i++){
//				b=(sarray[i]==0);
//			}
//			if(b){
//				lidx=Arrays.binarySearch(lims, 1000)-1;
//				lims[lidx]=0;
//			}
//		}
		if(skipDuplicateLines){
			int maxZero=-1;
			for(int i=0; i<sarray.length && sarray[i]==0; i++){
				maxZero=i;
			}
			if(maxZero>=50){
				for(int i=1; i<lims.length; i++){
					if(lims[i]<=maxZero){
						lidx=i-1;
						lims[lidx]=0;
					}else{break;}
				}
			}
		}
		
		final int lim=slist.size;
		assert(lim<=sarray.length);
		
		{
			//These two loops generate the scaffold table
			long prevSsum=-1, prevCsum=-1, prevSlen=-1, prevClen=-1;
			String prevLine=null;
			for(int i=0; i<lim && slen>0; i++){
				//			System.out.println("\n<A1>\ti="+i+", lidx="+lidx+", lims.length="+lims.length+", next="+next+", sarray.length="+sarray.length+", slen="+slen);

				if(i==next){
					prevLine=formatX(next, ssum, csum, slen, clen);
//					if(!skipDuplicateLines || ssum!=prevSsum || csum!=prevCsum || slen!=prevSlen || clen!=prevClen){
						sb2.append(prevLine);
						sb2.append('\n');
						prevLine=null;
//					}
					prevSsum=ssum; prevCsum=csum; prevSlen=slen; prevClen=clen;
					
					lidx++;
					if(lidx<lims.length){next=lims[lidx];}
					else{next=-1;}
					//				System.out.println("<A2>\ti="+i+", lidx="+lidx+", lims.length="+lims.length+", next="+next+", sarray.length="+sarray.length+", slen="+slen);
					//				System.out.println(sb2);
				}
				//			System.out.println("<A3>\ti="+i+", lidx="+lidx+", lims.length="+lims.length+", next="+next+", sarray.length="+sarray.length+", slen="+slen);

				if(i==50001){
					numOver50=ssum;
					fractionOver50=slen*100f/scaflen;
				}

				//			long a=carray[i];
				long b=sarray[i];
				long c=scarray1[i];
				long d=scarray2[i];

				if(b>0){
					csum-=c;
					ssum-=b;
					clen-=d;
					slen-=(b*i);
				}

				if(ln50==-1 && slen<=shalf){
					ln50=i;
					ll50=ssum+b;
				}
				//			System.out.println("<A4>\tb="+b+", c="+c+", d="+d+", csum="+csum+", ssum="+ssum+", clen="+clen+", slen="+slen);

				//			System.out.println("<A5>\ti="+i+", lidx="+lidx+", lims.length="+lims.length+", next="+next+", sarray.length="+sarray.length+", slen="+slen);

			}

			for(Triple tp : tlist){
				//			assert(false) : tlist;
				while(tp.length>=next && lidx<lims.length){
					
					prevLine=formatX(next, ssum, csum, slen, clen);
//					if(!skipDuplicateLines || ssum!=prevSsum || csum!=prevCsum || slen!=prevSlen || clen!=prevClen){
						sb2.append(prevLine);
						sb2.append('\n');
						prevLine=null;
//					}
					prevSsum=ssum; prevCsum=csum; prevSlen=slen; prevClen=clen;
					
					lidx++;
					if(lidx<lims.length){next=lims[lidx];}
					//				System.out.println("<B>\n"+sb2+"\ni="+"?"+", lidx="+lidx+", lims.length="+lims.length+", next="+next);
					//				else{next=-1;}
				}

				if(numOver50==0 && tp.length>50000){
					numOver50=ssum;
					fractionOver50=slen*100f/scaflen;
				}

				//			long a=carray[i];
				long b=tp.length;
				long c=tp.contigs;
				long d=tp.contiglen;

				if(b>0){
					csum-=c;
					ssum-=1;
					clen-=d;
					slen-=b;
				}

				if(ln50==-1 && slen<=shalf){
					ln50=b;
					ll50=ssum+1;
				}

			}
			if(prevLine!=null){
				sb2.append(prevLine);
				prevLine=null;
			}
		}
		
		clen=contiglen;
		csum=contigs;
		for(int i=0; i<carray.length && clen>0; i++){
			long a=carray[i];
			
			csum-=a;
			clen-=a*i;
			
			if(cn50==-1 && clen<=chalf){
				cn50=i;
				cl50=csum+a;
			}
		}

		for(int i=0; i<llist.size && clen>0; i++){
			long a=larray[i];
			
			csum-=1;
			clen-=a;
			
			if(cn50==-1 && clen<=chalf){
				cn50=a;
				cl50=csum+1;
			}
		}
		
		cn50=Tools.max(cn50, 0);
		ln50=Tools.max(ln50, 0);
		cl50=Tools.max(cl50, 0);
		ll50=Tools.max(ll50, 0);
		
		if(FORMAT<3){

			if(addfilename){sb.append("Filename:                           \t"+name+"\n");}
			sb.append("Main genome scaffold N/L50:         \t"+ll50+"/"+formatKB(ln50, 3, 0)+"\n");
			sb.append("Main genome contig N/L50:           \t"+cl50+"/"+formatKB(cn50, 3, 0)+"\n");
			sb.append("Max scaffold length:                \t"+formatKB(maxScaf, 3, 0)+"\n");
			sb.append("Max contig length:                  \t"+formatKB(maxContig, 3, 0)+"\n");
			sb.append("Number of scaffolds > 50 KB:        \t"+numOver50+"\n");
			sb.append("% main genome in scaffolds > 50 KB: \t"+String.format("%.2f%%", fractionOver50)+"\n");
			if(printheadersize){sb.append("Header:\t"+formatKB(HEADERLENSUM, 3, 0)+(HEADERLENSUM<1000 ? " bytes" : ""));}
			
			//		System.out.println();
			//		System.out.println("Scaffolds: "+Tools.sum(slist.array));
			//		for(int i=0; i<slist.size; i++){
			//			if(slist.array[i]>0){System.out.print(i+":"+slist.array[i]+", ");}
			//		}
			//		System.out.println();
			//		System.out.println("Contigs:"+Tools.sum(clist.array));
			//		for(int i=0; i<clist.size; i++){
			//			if(clist.array[i]>0){System.out.print(i+":"+clist.array[i]+", ");}
			//		}
			
			if(GCFORMAT==1 || GCFORMAT==3 || GCFORMAT==4){
				System.out.println("            \tA\tC\tG\tT\tN\tGC\tGC_stdev");
				System.out.println(toString3(new StringBuilder("Base Content"), counts));
			}else{
				System.out.println("            \tGC\tGC_stdev");
				System.out.println(toString3(new StringBuilder("Base Content"), counts));
			}

			System.out.println(sb);
			System.out.println(sb2);
		}else if(FORMAT==3 || FORMAT==6){
			
			if(useheader){
				if(FORMAT==6){sb.append('#');}
				if(N_UNDERSCORE){sb.append("n_");}
				sb.append("scaffolds\t");
				if(N_UNDERSCORE){sb.append("n_");}
				sb.append("contigs\t");
				sb.append("scaf_bp\t");
				sb.append("contig_bp\t");
				sb.append("gap_pct\t");
				sb.append("scaf_N50\t");
				sb.append("scaf_L50\t");
				sb.append("ctg_N50\t");
				sb.append("ctg_L50\t");
				sb.append("scaf_max\t");
				sb.append("ctg_max\t");
				sb.append("scaf_n_gt50K\t");
				sb.append("scaf_pct_gt50K\t");
				sb.append("gc_avg\t");
				sb.append("gc_std");
				if(addfilename){sb.append("\tfilename");}
				
				sb.append("\n");
			}
			
			sb.append(scaffolds+"\t");
			sb.append(contigs+"\t");
			sb.append(scaflen+"\t");
			sb.append(String.format("%d",contiglen)+"\t");
			sb.append(String.format("%.3f",(scaflen-contiglen)*100f/scaflen)+"\t");
			sb.append(ll50+"\t");
			sb.append(formatKB(ln50, 3, 0)+"\t");
			sb.append(cl50+"\t");
			sb.append(formatKB(cn50, 3, 0)+"\t");
			sb.append(formatKB(maxScaf, 3, 0)+"\t");
			sb.append(formatKB(maxContig, 3, 0)+"\t");
			sb.append(numOver50+"\t");
			sb.append(String.format("%.3f", fractionOver50)+"\t");
			sb.append(String.format("%.5f", (counts[1]+counts[2])*1.0/(counts[0]+counts[1]+counts[2]+counts[3]))+"\t");
			sb.append(String.format("%.5f", gc_stdS));
			if(addfilename){sb.append('\t').append(name);}
			
			System.out.println(sb);
		}else if(FORMAT==4){
			
			if(useheader){

			if(N_UNDERSCORE){sb.append("n_");}
			sb.append("scaffolds\t");
//			sb.append("contigs\t");
			sb.append("scaf_bp\t");
//			sb.append("contig_bp\t");
//			sb.append("gap_pct\t");
			sb.append("scaf_N50\t");
			sb.append("scaf_L50\t");
//			sb.append("ctg_N50\t");
//			sb.append("ctg_L50\t");
			sb.append("scaf_max\t");
//			sb.append("ctg_max\t");
			sb.append("scaf_n_gt50K\t");
			sb.append("scaf_pct_gt50K");
//			sb.append("gc_avg");
//			sb.append("gc_std");
			if(addfilename){sb.append("\tfilename");}

			sb.append("\n");
			}
			
			sb.append(scaffolds+"\t");
//			sb.append(contigs+"\t");
			sb.append(scaflen+"\t");
//			sb.append(String.format("%d",contiglen)+"\t");
//			sb.append(String.format("%.3f",(scaflen-contiglen)*100f/scaflen)+"\t");
			sb.append(ll50+"\t");
			sb.append(formatKB(ln50, 3, 0)+"\t");
//			sb.append(cl50+"\t");
//			sb.append(formatKB(cn50, 3, 0)+"\t");
			sb.append(formatKB(maxScaf, 3, 0)+"\t");
//			sb.append(formatKB(maxContig, 3, 0)+"\t");
			sb.append(numOver50+"\t");
			sb.append(String.format("%.3f", fractionOver50));
//			sb.append(String.format("%.5f", (counts[1]+counts[2])*1.0/(counts[0]+counts[1]+counts[2]+counts[3])));
			if(addfilename){sb.append('\t').append(name);}
			System.out.println(sb);
		}else if(FORMAT==5){
			
			if(useheader){
//			sb.append("scaffolds\t");
			if(N_UNDERSCORE){sb.append("n_");}
			sb.append("contigs\t");
//			sb.append("scaf_bp\t");
			sb.append("contig_bp\t");
			sb.append("gap_pct\t");
//			sb.append("scaf_N50\t");
//			sb.append("scaf_L50\t");
			sb.append("ctg_N50\t");
			sb.append("ctg_L50\t");
//			sb.append("scaf_max\t");
			sb.append("ctg_max\t");
//			sb.append("scaf_n_gt50K\t");
//			sb.append("scaf_pct_gt50K\t");
			sb.append("gc_avg\t");
			sb.append("gc_std");
			if(addfilename){sb.append("\tfilename");}

			sb.append("\n");
			}
			
//			sb.append(scaffolds+"\t");
			sb.append(contigs+"\t");
//			sb.append(scaflen+"\t");
			sb.append(String.format("%d",contiglen)+"\t");
			sb.append(String.format("%.3f",(scaflen-contiglen)*100f/scaflen)+"\t");
//			sb.append(ll50+"\t");
//			sb.append(formatKB(ln50, 3, 0)+"\t");
			sb.append(cl50+"\t");
			sb.append(formatKB(cn50, 3, 0)+"\t");
//			sb.append(formatKB(maxScaf, 3, 0)+"\t");
			sb.append(formatKB(maxContig, 3, 0)+"\t");
//			sb.append(numOver50+"\t");
//			sb.append(String.format("%.3f", fractionOver50)+"\t");
			sb.append(String.format("%.5f", (counts[1]+counts[2])*1.0/(counts[0]+counts[1]+counts[2]+counts[3]))+"\t");
			sb.append(String.format("%.5f", gc_stdS));
			if(addfilename){sb.append('\t').append(name);}
			System.out.println(sb);
		}
		
		if(showspeed){
			if(!printheadersize){System.err.println("Header:\t"+formatKB(HEADERLENSUM, 3, 0)+(HEADERLENSUM<1000 ? " bytes" : ""));}
			System.err.println("Time: \t"+t);
			long bytes=new File(in).length();
			if(bytes<1){bytes=LIMSUM;}
			double mbps=bytes*1000d/t.elapsed;
			double mbpps=Tools.sum(counts)*1000d/t.elapsed;
			System.err.println(String.format("Speed:\t%.2f MBytes/s",mbps));
			System.err.println(String.format("      \t%.2f MBases/s",mbpps));
		}
		
		if(bbmapkmer>0){
			System.err.println("BBMap minimum memory estimate at k="+bbmapkmer+":     "+estimateBBMapMemory(counts, scaffolds, HEADERLENSUM, bbmapkmer));
		}
		
		
	}
	
	
	/**
	 * @param counts
	 * @param scaffolds
	 * @param hEADERLENSUM2
	 * @param i
	 * @return
	 */
	private static long bbmapMemoryBytes(long[] acgtn, long scaffolds,
			long headerlen, int k) {

		long keyspace=(1L<<(2*k));
		long defined=acgtn[0]+acgtn[1]+acgtn[2]+acgtn[3];
		long undefined=acgtn[4];
		long midpad=(scaffolds*(FastaToChromArrays.MID_PADDING));
		long total=defined+undefined+midpad;
		int chromlen=FastaToChromArrays.MAX_LENGTH-FastaToChromArrays.END_PADDING-FastaToChromArrays.START_PADDING;
		int chroms=(int)(total/chromlen);
		int chromsperblock=Integer.MAX_VALUE/chromlen;
		int blocks=(chroms+chromsperblock-1)/chromsperblock;
		long memperblock=keyspace*4;
		long memforcounts=keyspace*4;
		
		long mem=0;
		mem+=total; //reference bulk, including inter-scaffold padding
		mem+=(chroms*(FastaToChromArrays.END_PADDING+FastaToChromArrays.START_PADDING)); //reference tip padding
		mem+=headerlen; //Header name byte arrays
		mem+=(scaffolds*(4+4+4+16+8)); //Other structures for scaffold info
		mem+=(blocks*(memperblock)); //start array for each block
		mem+=memforcounts; //count array
		mem+=(defined*4); //key lists
		mem=(long)(mem/0.66); //Expand to compensate for garbage collection
		if(k>13){mem=mem+1000000000;}
		return mem;
	}
	
	
	/**
	 * @param counts
	 * @param scaffolds
	 * @param hEADERLENSUM2
	 * @param i
	 * @return
	 */
	private static CharSequence estimateBBMapMemory(long[] acgtn, long scaffolds,
			long headerlen, int k) {
		long mem=180+bbmapMemoryBytes(acgtn, scaffolds, headerlen, k)/1000000; //in megabytes
		if(mem>4000){
			return "-Xmx"+((mem+1500)/1000)+"g \t(at least "+(long)Math.ceil((((mem+1500)/0.85)/1000))+" GB physical RAM)";
		}else if(mem>2100){
			return "-Xmx"+((mem+999)/1000)+"g \t(at least "+(long)Math.ceil((((mem+999)/0.85)/1000))+" GB physical RAM)";
		}else{
			return "-Xmx"+(((((mem*11)/8+50))/10)*10)+"m \t(at least "+((((long)(((mem*10)/8+50)/0.82))/10)*10)+" MB physical RAM)";
		}
	}


	public static long bench(InputStream is) throws IOException{
		final byte[] buf=new byte[32768];
		long sum=0;
		for(long len=is.read(buf); len>0; len=is.read(buf)){sum+=len;}
		return sum;
	}
	
	public static long[] countFasta(final InputStream is, final String out, final int maxNs) throws IOException{
		
		GCBINS2=(GCBINS>=1000 ? GCBINS : GCBINS*10);
		
		long limsum=0;
		long headerlen=0;
		final byte[] buf=new byte[32768]; 
		final TextStreamWriter tsw=(out==null ? null : new TextStreamWriter(out, overwrite, false, false));
		if(tsw!=null){tsw.start();}
		final int[] counts=new int[6];
		final long[] overall=new long[6];
		final StringBuilder hdr=(out==null ? null : new StringBuilder());
		boolean hdmode=false;
		
		final LongList clist=new LongList((int)Tools.min(1<<15, cutoff+1)); //Number of contigs of length x
		final LongList slist=new LongList((int)Tools.min(1<<15, cutoff+1)); //Number of scaffolds of length x
		final LongList sclist1=new LongList((int)Tools.min(1<<15, cutoff+1)); //Sum of contigs per scaffold of length x
		final LongList sclist2=new LongList((int)Tools.min(1<<15, cutoff+1)); //Sum of contig lengths per scaffold of length x
		
		final LongList llist=new LongList(64); //List of contig lengths for contigs at least cutoff in length
		final ArrayList<Triple> tlist=new ArrayList<Triple>(64); //List of scaf len, contigs, contig sum for scaffolds at least cutoff in length
		
//		final long[] gchist=(calcGChist ? new long[GCBINS] : null);
		final long[] gchist=new long[GCBINS2];
		final long[] gchist_by_base=new long[GCBINS2];
		
		int i=0;
		int lim=is.read(buf);
		limsum+=lim;
		
		int contigs=0;
		int contiglen=0;
//		int contiglensum=0;
		int scaffoldlen=0;
		int ns=0;
		
		final IntList currentContigs=new IntList(10000);
		
		while(lim>0){
			if(hdmode){
				if(hdr==null){
					while(i<lim){
						byte c=buf[i];
						i++;
						if(c<=slashr){
							hdmode=false;
							contiglen=0;
//							contiglensum=0;
							scaffoldlen=0;
							ns=0;
							contigs=0;
							break;
						}
						headerlen++;
					}
				}else{
					while(i<lim){
						byte c=buf[i];
						i++;
						if(c<=slashr){
							hdmode=false;
							contiglen=0;
//							contiglensum=0;
							scaffoldlen=0;
							ns=0;
							contigs=0;
							break;
						}
						hdr.append((char)c);
					}
				}
			}
			
			if(!hdmode){
				while(i<lim){
					byte c=buf[i];
					i++;

					if(c==carrot){
//						assert(false) : scaffoldlen;
						hdmode=true;
						if(scaffoldlen>0){

//							if(contiglen>0 || contigs==0){
//								contigs++;
//								contiglensum+=contiglen;
//								if(contiglen<cutoff){
//									clist.increment(contiglen, 1);
//								}else{
//									llist.add(contiglen);
//								}
//							}
							
							{//NEW
								if(contiglen>0 || contigs==0){
									currentContigs.set(contigs, contiglen);
									contigs++;
//									System.out.println("For header "+hdr+": added contig.  len="+contiglen+", contigs="+contigs);
//									contiglensum+=contiglen;
								}
							}
							
//							assert(false);
							if(scaffoldlen>=MINSCAF){

								int contiglensum=0;
								{//NEW
//									System.out.println("Dumping "+contigs+" contigs.");
									for(int j=0; j<contigs; j++){
										final int cl=currentContigs.get(j);
										if(cl>0 || contigs==0){
											contiglensum+=cl;
											if(cl<cutoff){
												clist.increment(cl, 1);
											}else{
												llist.add(cl);
											}
										}
									}
								}
								
								if(scaffoldlen<cutoff){
									slist.increment(scaffoldlen, 1);
									sclist1.increment(scaffoldlen, contigs);
									sclist2.increment(scaffoldlen, contiglensum);
								}else{
									tlist.add(new Triple(scaffoldlen, contigs, contiglensum));
								}


								if(hdr!=null){
									tsw.print(toString2(hdr, counts));
									headerlen+=hdr.length();
									hdr.setLength(0);
								}
								{
									long gc=counts[1]+counts[2];
									long atgc=gc+counts[0]+counts[3];
									if(atgc>0){
										int index=Tools.min((int)((gc*GCBINS2)/atgc),GCBINS2-1);
										gchist[index]++;
										gchist_by_base[index]+=scaffoldlen;
//										assert(false);
									}
								}
								for(int j=0; j<counts.length; j++){
									overall[j]+=counts[j];
									counts[j]=0;
								}
							}else{
								Arrays.fill(counts, 0);
								if(hdr!=null){hdr.setLength(0);}
							}
						}

						break;
					}
					
					if(c>slashr){
						counts[charToNum[c]]++;
						scaffoldlen++;
						
						if(c!=noref && c!=noref2){
							ns=0;
							contiglen++;
						}else{
							ns++;
							if(ns==maxNs && contiglen>0){
//								if(contiglen<cutoff){
//									clist.increment(contiglen, 1);
//								}else{
//									llist.add(contiglen);
//								}
////								clist.increment(contiglen, 1);
//								contiglensum+=contiglen;
//								contiglen=0;
//								contigs++;

								{//NEW
									currentContigs.set(contigs, contiglen);
									contiglen=0;
									contigs++;
								}
							}
						}
					}
				}
			}
			if(i>=lim){
				i=0;
				lim=is.read(buf);
				limsum+=lim;
			}
		}

		if(scaffoldlen>0){
			
//			if(contiglen>0 || contigs==0){
//				contigs++;
//				contiglensum+=contiglen;
//				if(contiglen<cutoff){
//					clist.increment(contiglen, 1);
//				}else{
//					llist.add(contiglen);
//				}
//			}
			
			{//NEW
				if(contiglen>0 || contigs==0){
					currentContigs.set(contigs, contiglen);
					contigs++;
//					contiglensum+=contiglen;
				}
			}
			
			if(scaffoldlen>=MINSCAF){

				int contiglensum=0;
				{//NEW
//					System.out.println("Dumping "+contigs+" contigs.");
					for(int j=0; j<contigs; j++){
						final int cl=currentContigs.get(j);
						if(cl>0 || contigs==0){
							contiglensum+=cl;
							if(cl<cutoff){
								clist.increment(cl, 1);
							}else{
								llist.add(cl);
							}
						}
					}
				}
				
				if(scaffoldlen<cutoff){
					slist.increment(scaffoldlen, 1);
					sclist1.increment(scaffoldlen, contigs);
					sclist2.increment(scaffoldlen, contiglensum);
				}else{
					tlist.add(new Triple(scaffoldlen, contigs, contiglensum));
				}


//				slist.increment(scaffoldlen, 1);
//				if(contiglen>0 || contigs==0){
//					contigs++;
//					contiglensum+=contiglen;
//					clist.increment(contiglen, 1);
//				}
//				sclist1.increment(scaffoldlen, contigs);
//				sclist2.increment(scaffoldlen, contiglensum);

				if(hdr!=null){
					tsw.print(toString2(hdr, counts));
					hdr.setLength(0);
				}

				{
					long gc=counts[1]+counts[2];
					long atgc=gc+counts[0]+counts[3];
					if(atgc>0){
						int index=Tools.min((int)((gc*GCBINS2)/atgc),GCBINS2-1);
						gchist[index]++;
						gchist_by_base[index]+=scaffoldlen;
					}
				}
				for(int j=0; j<counts.length; j++){
					overall[j]+=counts[j];
					counts[j]=0;
				}
			}
		}
		
//		System.err.println("clist="+clist+"\nslist="+slist+"\nsclist1="+sclist1+"\nsclist2="+sclist2+"\nllist="+llist+"\ntlist="+tlist); //***
		
		
		if(tsw!=null){
			tsw.poison();
			tsw.waitForFinish();
		}
		LIMSUM=limsum;
		HEADERLENSUM=headerlen;
		
		clistS=clist;
		slistS=slist;
		sclist1S=sclist1;
		sclist2S=sclist2;

		tlistS=tlist;
		llistS=llist;

		gc_stdS=Tools.standardDeviationHistogram(gchist)/GCBINS2;
		gchistS=Tools.downsample(gchist, GCBINS);

		gc_bb_stdS=Tools.standardDeviationHistogram(gchist_by_base)/GCBINS2;
		gchist_bb_S=Tools.downsample(gchist_by_base, GCBINS);
		
		return overall;
	}
	
	private static void writeHistFile(String fname, LongList slist, ArrayList<Triple> tlist, boolean ascending){
		if(fname==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#scaffolds\tlength\n");
		long num=0, len=0;
		
		final long[] array=slist.array;
		final int lim=slist.size;
		
		StringBuilder sb=new StringBuilder(32);
		
		if(ascending){
			for(int i=0; i<lim; i++){
				final long a=array[i];
				if(a>0){
					len+=(a*i);
					num+=a;
					sb.append(num);
					sb.append('\t');
					sb.append(len);
					sb.append('\n');
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}

			if(tlist!=null){
				for(Triple t : tlist){
					len+=t.length;
					num++;
					sb.append(num);
					sb.append('\t');
					sb.append(len);
					sb.append('\n');
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}
		}else{

			if(tlist!=null){
				for(int i=tlist.size()-1; i>=0; i--){
					Triple t=tlist.get(i);
					len+=t.length;
					num++;
					sb.append(num);
					sb.append('\t');
					sb.append(len);
					sb.append('\n');
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}
			for(int i=lim-1; i>=0; i--){
				final long a=array[i];
				if(a>0){
					len+=(a*i);
					num+=a;
					sb.append(num);
					sb.append('\t');
					sb.append(len);
					sb.append('\n');
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}
		}
		tsw.poisonAndWait();
	}
	
	private static String toString2(StringBuilder sb, int[] counts){
		long sum=(long)counts[0]+(long)counts[1]+(long)counts[2]+(long)counts[3];
		long sum2=sum+counts[4];
		float inv1=1f/sum;
		if(GCFORMAT==1){
			return sb.append(String.format("\t0\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n", 
					sum, counts[0]*inv1, counts[1]*inv1, counts[2]*inv1, counts[3]*inv1, counts[4]*1f/sum2, (counts[1]+counts[2])*inv1)).toString();
		}else if(GCFORMAT==2){
			return sb.append(String.format("\t%.5f\n", (counts[1]+counts[2])*inv1)).toString();
		}else if(GCFORMAT==3){
			return sb.append(String.format("\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n", 
					sum, counts[0]*inv1, counts[1]*inv1, counts[2]*inv1, counts[3]*inv1, counts[4]*1f/sum2, (counts[1]+counts[2])*inv1)).toString();
		}else if(GCFORMAT==4){
			return sb.append(String.format("\t%d\t%.5f\n", sum, (counts[1]+counts[2])*inv1)).toString();
		}else{
			throw new RuntimeException("Unknown format.");
		}
	}
	
	private static String toString2(StringBuilder sb, long[] counts){
		long sum=(long)counts[0]+(long)counts[1]+(long)counts[2]+(long)counts[3];
		long sum2=sum+counts[4];
		float inv1=1f/sum;
		if(GCFORMAT==1){
			return sb.append(String.format("\t0\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n", 
					sum, counts[0]*inv1, counts[1]*inv1, counts[2]*inv1, counts[3]*inv1, counts[4]*1f/sum2, (counts[1]+counts[2])*inv1)).toString();
		}else if(GCFORMAT==2){
			return sb.append(String.format("\t%.5f\n", (counts[1]+counts[2])*inv1)).toString();
		}else if(GCFORMAT==3){
			return sb.append(String.format("\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n", 
					sum, counts[0]*inv1, counts[1]*inv1, counts[2]*inv1, counts[3]*inv1, counts[4]*1f/sum2, (counts[1]+counts[2])*inv1)).toString();
		}else if(GCFORMAT==4){
			return sb.append(String.format("\t%d\t%.5f\n", sum, (counts[1]+counts[2])*inv1)).toString();
		}else{
			throw new RuntimeException("Unknown format.");
		}
	}
	
	private static String toString3(StringBuilder sb, long[] counts){
		long sum=(long)counts[0]+(long)counts[1]+(long)counts[2]+(long)counts[3];
		long sum2=sum+counts[4];
		float inv1=1f/sum;
		if(GCFORMAT==1 || GCFORMAT==3 || GCFORMAT==4){
			return sb.append(String.format("\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", 
					counts[0]*inv1, counts[1]*inv1, counts[2]*inv1, counts[3]*inv1, counts[4]*1f/sum2, (counts[1]+counts[2])*inv1, gc_stdS)).toString();
		}else if(GCFORMAT==2){
			return sb.append(String.format("\t%.5f\t%.5f\n", (counts[1]+counts[2])*inv1, gc_stdS)).toString();
		}else{
			throw new RuntimeException("Unknown format.");
		}
		
	}
	
	
	private static final String formatX(int next, long ssum, long csum, long slen, long clen){
		float cov=clen*100f/slen;
		
		final String s;
		if(FORMAT==1){
			s=formatKB_all(next, 1, 7)+" \t"+formatComma(ssum, 14)+"\t"+formatComma(csum, 14)+"\t"+formatComma(slen, 14)+"\t"+formatComma(clen, 14)+"\t"+formatPercent(cov);
		}else if(FORMAT>=2){
			s=formatKB_all(next, 1, 7)+" \t"+formatComma(ssum, 14)+"\t"+formatComma(csum, 14)+"\t"+formatComma(slen, 14)+"\t"+formatComma(clen, 14)+"\t"+formatPercent(cov);
		}else{
			throw new RuntimeException("Unknown format: "+FORMAT);
		}
		return s;
	}
	
	private static final String formatKB(long x, int precision, int width){
		String s;
		if(FORMAT>=2 || x<1000){
			s=Long.toString(x);
		}else if(x<1000000){
			s=String.format("%."+precision+"f",x/1000f);
			while(s.contains(".") && (s.endsWith("0") || s.endsWith("."))){s=s.substring(0, s.length()-1);}
			s=s+" KB";
		}else{
			s=String.format("%."+precision+"f",x/1000000f);
			while(s.contains(".") && (s.endsWith("0") || s.endsWith("."))){s=s.substring(0, s.length()-1);}
			s=s+" MB";
		}
		while(s.length()<width){s=" "+s;}
		return s;
	}
	
	private static final String formatKB_all(long x, int precision, int width){
		String s;
			
		if(x==0){s="All";}
		else if(FORMAT>=2 || x<1000){
			s=Long.toString(x);
		}else if(x<1000000){
			s=String.format("%."+precision+"f",x/1000f);
			while(s.contains(".") && (s.endsWith("0") || s.endsWith("."))){s=s.substring(0, s.length()-1);}
			s=s+" KB";
		}else{
			s=String.format("%."+precision+"f",x/1000000f);
			while(s.contains(".") && (s.endsWith("0") || s.endsWith("."))){s=s.substring(0, s.length()-1);}
			s=s+" MB";
		}
		
		while(s.length()<width){s=" "+s;}
		return s;
	}
	
	private static final StringBuilder formatComma(long x, int width){
		StringBuilder sb=new StringBuilder(width);
		if(FORMAT==1){
			sb.append(x%1000);
			x/=1000;
			int len=3;
			while(x>0){
				while(sb.length()<len){sb.insert(0, '0');}
				sb.insert(0, (x%1000)+",");
				x/=1000;
				len+=4;
			}
		}else if(FORMAT>=2){
			sb.append(x);
		}else{
			throw new RuntimeException("Unknown format: "+FORMAT);
		}
		while(sb.length()<width){
			sb.insert(0, ' ');
		}
		return sb;
	}
	
	private static final String formatPercent(float x){
		String s=String.format("%.2f%%", x);
		while(s.length()<8){s=" "+s;}
		return s;
	}
	
	protected static void reset(){
		clistS=null;
		slistS=null;
		sclist1S=null;
		sclist2S=null;

		gchistS=null;
		gc_stdS=0;

		llistS=null;
		tlistS=null;
	}
	
	/**
	 * @return
	 */
	private static byte[] makeCharToNum() {
		byte[] r=new byte[256];
		Arrays.fill(r, (byte)4);
		r['a']=r['A']=0;
		r['c']=r['C']=1;
		r['g']=r['G']=2;
		r['t']=r['T']=3;
		r['\n']=r['\r']=r['>']=r['@']=r['+']=5;
		return r;
	}
	
	private static final byte[] charToNum=makeCharToNum();
	public static int GCFORMAT=1;
	public static int FORMAT=1;
	private static long cutoff=1000000;
	
	private static long LIMSUM=0;
	private static long HEADERLENSUM=0;
	private static int GCBINS=200;
	private static int GCBINS2;
	private static int MINSCAF=0;
	private static int GCHIST_DECIMALS=3;
	
	private static int bbmapkmer=0;//13;
	public static boolean overwrite=false;
	public static boolean useheader=true;
	public static boolean addfilename=false;
	public static boolean showspeed=false;//true;
	public static boolean printheadersize=false;
	public static boolean skipDuplicateLines=true;
	public static boolean N_UNDERSCORE=true;

	private final static byte slashr='\r', slashn='\n', carrot='>', at='@', noref='N', noref2='n';
	
	/** Number of contigs of length x */
	private static LongList clistS=null;
	
	/** Number of scaffolds of length x */
	private static LongList slistS=null;
	
	/** Sum of contigs per scaffold of length x */
	private static LongList sclist1S=null;
	
	/** Sum of contig lengths per scaffold of length x */
	private static LongList sclist2S=null;
	
	/** List of contig lengths for contigs at least cutoff in length */
	private static LongList llistS=null;
	
	/** List of scaf len, contigs, contig sum for scaffolds at least cutoff in length */
	private static ArrayList<Triple> tlistS=null;

	/** Downsampled gc histogram */
	private static long[] gchistS;
	
	/** gc standard deviation */
	private static double gc_stdS;
	
	/** Downsampled gc histogram, using base counts rather than scaffold counts */
	private static long[] gchist_bb_S;
	
	/** gc standard deviation, using base counts rather than scaffold counts */
	private static double gc_bb_stdS;
	
	private static class Triple implements Comparable<Triple>{
		
		Triple(long len_, long contigs_, long contiglen_){
			length=len_;
			contigs=contigs_;
			contiglen=contiglen_;
		}
		
		@Override
		public int compareTo(Triple o) {
			if(length>o.length){return 1;}
			if(length<o.length){return -1;}
			return (int)(contiglen-o.contiglen);
		}
		
		public boolean equals(Object o){return equals((Triple)o);}
		
		public boolean equals(Triple o){
			return length==o.length && contiglen==o.contiglen;
		}
		
		public String toString(){return length+","+contigs+","+contiglen;}
		
		public final long length;
		public final long contigs;
		public final long contiglen;
		
	}
	
}
