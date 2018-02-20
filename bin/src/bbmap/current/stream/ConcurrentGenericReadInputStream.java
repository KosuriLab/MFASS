package stream;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

import align2.ListNum;
import align2.Shared;
import align2.Tools;

import dna.Data;
import dna.Timer;

import fileIO.ReadWrite;
import fileIO.FileFormat;

public class ConcurrentGenericReadInputStream implements ConcurrentReadStreamInterface {
	
	public static void main(String[] args){
		String in1=args[0];
		String in2=(args.length<2 || args[1].equalsIgnoreCase("null") || args[1].contains("=") ? null : args[1]);
		if(in2!=null){
			assert(!in1.equalsIgnoreCase(in2));
			FASTQ.TEST_INTERLEAVED=false;
		}else{
			FASTQ.TEST_INTERLEAVED=true;
			FASTQ.FORCE_INTERLEAVED=true;
		}
		
		long maxReads=-1;
		for(int i=1; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=(split.length>1 ? split[1] : "true");
			
			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("null") || (split.length==1 && i==1)){
				// do nothing
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.startsWith("fastareadlen")){
				FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
				FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
			}else if(a.startsWith("fastaminread") || a.startsWith("fastaminlen")){
				FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		assert(FastaReadInputStream.settingsOK());
		Timer t=new Timer();
		t.start();
		
		ConcurrentReadStreamInterface cris=getReadInputStream(maxReads, false, false, true, in1, in2);
		System.out.println("Fetched "+cris.getClass().getName());
		{
			Object[] p=cris.producers();
//			while(p[0]==null){
//				p=cris.producers();
//			}
			System.out.print("Producers: ");
			String comma="";
			for(Object o : p){
				System.out.print(comma+(o==null ? "null" : o.getClass().getName()));
				comma=", ";
			}
			System.out.println();
		}
		boolean paired=cris.paired();
		System.out.println("paired="+paired);
		Thread cristhread=new Thread(cris);
		cristhread.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((r.mate!=null)==paired);
		}
		
		long readCount=0;
		long baseCount=0;
		
		while(reads!=null && reads.size()>0){
			
			for(Read r : reads){
				Read r2=r.mate;
				if(r!=null){
					readCount++;
					if(r.bases!=null){
						baseCount+=r.bases.length;
					}
				}
				if(r2!=null){
					readCount++;
					if(r2.bases!=null){
						baseCount+=r2.bases.length;
					}
				}	
			}
			cris.returnList(ln, ln.list.isEmpty());
//			System.err.println("fetching list");
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
//			System.out.println("reads: "+(reads==null ? "null" : reads.size()));
		}
		System.err.println("Finished reading");
		cris.returnList(ln, ln.list.isEmpty());
		
		cris.close();
		t.stop();

		System.out.println("Reads:      \t"+readCount);
		System.out.println("Bases:      \t"+baseCount);
		System.out.println("Avg Length: \t"+String.format("%.2f",baseCount*1.0/readCount));
		System.out.println("Time:      \t"+t);
	}
	
	public ConcurrentGenericReadInputStream(ReadInputStream source1, ReadInputStream source2, long maxReadsToGenerate){
		assert(source1!=source2);
		producer1=source1;
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		producer2=source2;
		assert(source2==null || !FASTQ.FORCE_INTERLEAVED) : "Please do not set 'interleaved=true' with dual input files.";
		maxReads=maxReadsToGenerate>=0 ? maxReadsToGenerate : Long.MAX_VALUE;
		if(maxReads==0){
			System.err.println("Warning - created a read stream for 0 reads.");
			assert(false);
		}
//		if(maxReads<Long.MAX_VALUE){System.err.println("maxReads="+maxReads);}

		if(producer1!=null){p1q=new ArrayBlockingQueue<ArrayList<Read>>(4);}
		if(producer2!=null){p2q=new ArrayBlockingQueue<ArrayList<Read>>(4);}
	}
	
	public synchronized ListNum<Read> nextList() {
		ArrayList<Read> list=null;
		if(verbose){System.err.println("**************** nextList() was called; shutdown="+shutdown+", depot.full="+depot.full.size());}
		while(list==null){
			if(shutdown){
				if(verbose){System.err.println("**************** nextList() returning null; shutdown="+shutdown+", depot.full="+depot.full.size());}
				return null;
			}
			try {
				list=depot.full.take();
				assert(list!=null);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		if(verbose){System.err.println("**************** nextList() returning list of size "+list.size()+"; shutdown="+shutdown+", depot.full="+depot.full.size());}
		ListNum<Read> ln=new ListNum<Read>(list, listnum);
		listnum++;
		return ln;
	}
	
	public void returnList(ListNum<Read> ln, boolean poison){
		if(ln!=null){
			ln.list.clear();
		}else{
			System.err.println("Warning, null list returned: ");  //System.err.println("Warning from class "+getClass().getName()+", null list returned: ");
			new Exception().printStackTrace();
		}
		if(poison){
			if(verbose){System.err.println("A: Adding empty list to full.");}
			depot.full.add(ln==null ? new ArrayList<Read>(0) : ln.list);
		}else{
			if(ln!=null){depot.empty.add(ln.list);}
//			depot.empty.add(ln==null ? new ArrayList<Read>(0) : ln.list);
		}
	}
	
	@Override
	public void run() {
//		producer.start();
		synchronized(running){
			assert(!running[0]) : "This cris was started by multiple threads.";
			running[0]=true;
		}

		ReadThread rt1=null;
		ReadThread rt2=null;
		if(producer1.preferLists() || producer1.preferBlocks()){
			rt1=new ReadThread(producer1, p1q);
			rt2=(producer2==null ? null : new ReadThread(producer2, p2q));
			rt1.start();
			if(rt2!=null){rt2.start();}
		}
		
		threads=(rt1==null ? new Thread[] {Thread.currentThread()} : 
			rt2==null ? new Thread[] {Thread.currentThread(), rt1} : 
				new Thread[] {Thread.currentThread(), rt1, rt2});

		if(producer1.preferLists() || producer1.preferBlocks()){
			readLists();
			//System.err.println("Done reading lists.");
		}else if(producer1.preferBlocks()){
			assert(false);
//			readBlocks();
		}else{
			readSingles();
		}
		
		addPoison();
		
		//End thread

		if(verbose){System.err.println("cris finished addPoison.");}
		while(!depot.empty.isEmpty() && !shutdown){
//			System.out.println("Ending");
			if(verbose){System.err.println("B: Adding empty lists to full.");}
			depot.full.add(depot.empty.poll());
		}
		if(verbose){System.err.println("cris thread syncing before shutdown.");}
		
		synchronized(running){//TODO Note: for some reason syncing on 'this' instead of 'running' causes a hang.  Something else must be syncing improperly on this.
			assert(running[0]);
			running[0]=false;
		}
		if(verbose){System.err.println("cris thread terminated. Final depot size: "+depot.full.size()+", "+depot.empty.size());}
	}
	
	private final void addPoison(){
		//System.err.println("Adding poison.");
		//Add poison pills
		if(verbose){System.err.println("C: Adding poison to full.");}
		depot.full.add(new ArrayList<Read>());
		for(int i=1; i<depot.bufferCount; i++){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=depot.empty.poll(1000, TimeUnit.MILLISECONDS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
//					System.err.println("Do not be alarmed by the following error message:");
//					e.printStackTrace();
					if(shutdown){
						i=depot.bufferCount;
						break;
					}
				}
			}
			if(list!=null){
				if(verbose){System.err.println("D: Adding list("+list.size()+") to full.");}
				depot.full.add(list);
			}
		}
		if(verbose){System.err.println("Added poison.");}
	}
	
	private final void readSingles(){

		while(!shutdown && producer1.hasMore() && generated<maxReads){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=depot.empty.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					if(shutdown){break;}
				}
			}
			if(shutdown || list==null){break;}
			
			long bases=0;
			while(list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
				Read a=producer1.next();
				Read b=(producer2==null ? null : producer2.next());
				if(a==null){break;}
				if(randy==null || randy.nextFloat()<samplerate){
					list.add(a);
					if(b!=null){
						assert(a.numericID==b.numericID) : "\n"+a.numericID+", "+b.numericID+"\n"+a.toText(false)+"\n"+b.toText(false)+"\n";
						assert(a.mate==null) : "Please set interleaved=false when using dual input files.\n"+a.id+"\n"+a.mate.id+"\n"+producer1+"\n"+producer2;
						assert(b.mate==null) : "Please set interleaved=false when using dual input files.";
						a.mate=b;
						b.mate=a;

						assert(a.pairnum()==0);
						b.setPairnum(1);
						bases+=(b.bases==null ? 0 : b.bases.length);
					}
					bases+=(a.bases==null ? 0 : a.bases.length);
				}
				incrementGenerated(1);
			}

			if(verbose){System.err.println("E: Adding list("+list.size()+") to full.");}
			depot.full.add(list);
		}
	}
	
	private final void readLists(){
		ArrayList<Read> buffer1=null;
		ArrayList<Read> buffer2=null;
		ArrayList<Read> list=null;
		int next=0;
		
//		System.out.println("a");
		if(verbose){System.err.println(getClass().getName()+" entering read lists loop.");}
		while(buffer1!=poison && (buffer1!=null || (!shutdown && generated<maxReads))){
//			System.out.println("b");
			if(verbose){System.err.println("looping: buffer1==null "+(buffer1==null)+", buffer1==poison "+(buffer1==poison)
					+", shutdown="+shutdown+", generated<maxReads="+(generated<maxReads));}
			while(list==null){
				if(verbose){System.err.println("Fetching an empty list: generated="+generated+"/"+maxReads);}
				try {
					list=depot.empty.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					if(shutdown){break;}
				}
				if(verbose){System.err.println("Fetched "+(list==null ? "null" : ""+list.size()));}
			}
//			System.out.println("c");
			if(verbose){System.err.println("Left empty fetch loop.");}
			if(shutdown || list==null){
				//System.err.println("Shutdown triggered; breaking.");
				break;
			}
//			System.out.println("d");
			
			if(verbose){System.err.println("Entering full fetch loop.");}
			long bases=0;
			while(list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
				if(verbose){System.err.println("list.size()="+list.size()+", depot.bufferSize="+depot.bufferSize+", generated="+generated);}
				if(buffer1==null || next>=buffer1.size()){
					buffer1=null;
					while(!shutdown && buffer1==null){
						try {
							buffer1=p1q.take();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
//					System.out.println("e");
					
					if(buffer1!=null && p2q!=null){
						buffer2=null;
						while(!shutdown && buffer2==null){
							try {
								buffer2=p2q.take();
							} catch (InterruptedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
						if(buffer2!=null){pair(buffer1, buffer2);}
						if(REMOVE_DISCARDED_READS){removeDiscarded(buffer1, buffer2);}
					}
//					System.out.println("f");
					next=0;
				}
//				System.out.println("g");
				if(buffer1==null || buffer1==poison || shutdown){
//					if(list!=null && list.size()>0){
//						if(verbose){System.err.println("G: Adding list("+list.size()+") to full.");}
//						depot.full.add(list);
//						list=null;
//					}
					if(verbose){System.err.println("Breaking because buffer1==null: "+(buffer1==null)+" || buffer1==poison: "+(buffer1==poison)+" || shutdown: "+shutdown);}
					break;
				}
				assert(buffer1.size()<=BUF_LEN); //Although this is not really necessary.
				
//				assert(!set2.contains(buffer1)) : buffer1.hashCode();
//				set2.add(buffer1);
//				System.out.println(buffer1.hashCode());
				
				if(buffer2!=null){
//					System.out.println("h");
					
					if(buffer2!=null && (buffer1==null || buffer2.size()!=buffer1.size())){
						System.err.println("Error: Misaligned read streams.");
						errorState=true;
						return;
					}
					assert(buffer2==null || buffer2.size()==buffer1.size());
				}
//				System.out.println("i");
				if(buffer1.size()<=(BUF_LEN-list.size()) && (buffer1.size()+generated)<maxReads && randy==null){
					//System.out.println("j");
					//Then do a quicker bulk operation
					
					list.addAll(buffer1);
					incrementGenerated(buffer1.size());
					for(Read a : buffer1){
						bases+=(a.bases==null ? 0 : a.bases.length);
//						bases+=(a.mate==null || a.mate.bases==null ? 0 : a.mate.bases.length);
						if(a.mate!=null){
							bases+=a.mate.bases.length;
							assert(a.pairnum()==0 && a.mate.pairnum()==1);
						}
					}
					
					next=0;
					buffer1=null;
					buffer2=null;
				}else{

					while(next<buffer1.size() && list.size()<depot.bufferSize && generated<maxReads && bases<MAX_DATA){
						Read a=buffer1.get(next);
						
						if(p2q!=null){
							Read b=buffer2.get(next);
//							assert(a.numericID==b.numericID) : "\n"+a.numericID+", "+b.numericID+"\n"+a.toText(false)+"\n"+b.toText(false)+"\n";
//							a.mate=b;
//							b.mate=a;
//							
//							assert(a.pairnum()==0);
//							b.setPairnum(1);
							assert(a.pairnum()==0 && b.pairnum()==1 && a.mate==b && b.mate==a && a.numericID==b.numericID);
						}
						
						if(randy==null || randy.nextFloat()<samplerate){
							list.add(a);
							bases+=a.bases.length;
							if(a.mate!=null){
								bases+=a.mate.bases.length;
//								assert(a.pairnum()==0 && a.mate.pairnum()==1);
							}
						}
						incrementGenerated(1);
						next++;
					}
//					System.out.println("l");

					
					if(next>=buffer1.size()){
						buffer1=null;
						buffer2=null;
						next=0;
//						System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
					}else{
//						System.out.println("------------------------------------------------");
					}
//					System.out.println("m");
				}
				if(verbose){System.err.println("Loop end: list.size()="+(list.size()+", depot.bufferSize="+depot.bufferSize+", generated="+generated));}
//				System.out.println("n");
				if(verbose){System.err.println(Thread.currentThread().getName());}
			}
			
//			System.out.println("p");
//			System.err.println("Adding list to full depot.  Shutdown="+shutdown);
			if(verbose){System.err.println("F: Adding list("+list.size()+") to full.");}
			depot.full.add(list);
//			System.err.println("Added.");
			
//			System.out.println("o");
			if(buffer1==poison){
				if(verbose){System.err.println("Detected poison from buffer1.");}
				break;
			}
			list=null;
			if(verbose){System.err.println("Finished loop iteration.\n");}
			if(verbose){System.err.println("loop end: buffer1==null "+(buffer1==null)+", buffer1==poison "+(buffer1==poison)
					+", shutdown="+shutdown+", generated<maxReads="+(generated<maxReads));}
//			System.out.println("q");
		}
//		System.out.println("r");
		
		
		p1q.clear();
		if(p2q!=null){p2q.clear();}
	}
	
	private final void pair(ArrayList<Read> buffer1, ArrayList<Read> buffer2){
		for(int i=0; i<buffer1.size(); i++){
			Read a=buffer1.get(i);
			Read b=buffer2.get(i);

			assert(a.numericID==b.numericID) : "\n"+a.numericID+", "+b.numericID+"\n"+a.toText(false)+"\n"+b.toText(false)+"\n";
			assert(a.mate==null) : "Please set interleaved=false when using dual input files.\n"+a.id+"\n"+a.mate.id+"\n"+b.id+"\n"+producer1+"\n"+producer2;
			assert(b.mate==null) : "Please set interleaved=false when using dual input files.";
			a.mate=b;
			b.mate=a;

			assert(a.pairnum()==0);
			b.setPairnum(1);
			//		assert(a.pairnum()!=b.pairnum());
		}
	}
	
	private final int removeDiscarded(ArrayList<Read> buffer1, ArrayList<Read> buffer2){
		int removed=0;
		if(buffer2==null){
			for(int i=0; i<buffer1.size(); i++){
				Read a=buffer1.get(i);
				if(a.discarded()){
					buffer1.set(i, null);
					removed++;
				}
			}
		}else{
			for(int i=0; i<buffer1.size(); i++){
				Read a=buffer1.get(i);
				Read b=buffer2.get(i);
				if(a.discarded() || b.discarded()){
					buffer1.set(i, null);
					buffer2.set(i, null);
					removed++;
				}
			}
		}
		if(removed>0){
			Tools.condenseStrict(buffer1);
			if(buffer2!=null){Tools.condenseStrict(buffer2);}
		}
		return removed;
	}
	
	private boolean shutdown=false;
	
	@Override
	public void shutdown(){
//		System.err.println("Called shutdown.");
		shutdown=true;
		if(!shutdown){
			for(Thread t : threads){
				if(t!=null && t.isAlive()){
					t.interrupt();
				}
			}
		}
	}
	
	@Override
	public synchronized void restart(){
		shutdown=false;
		p1q.clear();
		if(p2q!=null){p2q.clear();}
		producer1.restart();
		if(producer2!=null){producer2.restart();}
		depot=new ConcurrentDepot<Read>(BUF_LEN, NUM_BUFFS);
		generated=0;
		nextProgress=PROGRESS_INCR;
	}
	
	@Override
	public synchronized void close(){
		if(verbose){System.err.println("Called shutdown for "+producer1+"; "+threads[0].getState());}
//		if(verbose){System.err.println(((FastqReadInputStream)producer1).tf.isOpen());}
		shutdown();
		errorState|=producer1.close();
		if(producer2!=null){errorState|=producer2.close();}
		if(threads!=null && threads[0]!=null && threads[0].isAlive()){
			
			while(threads[0].isAlive()){
//				System.out.println("B");
				ArrayList<Read> list=null;
				for(int i=0; i<1000 && list==null && threads[0].isAlive(); i++){
					try {
						list=depot.full.poll(200, TimeUnit.MILLISECONDS);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						System.err.println("Do not be alarmed by the following error message:");
						e.printStackTrace();
						break;
					}
				}
				
				if(list!=null){
					list.clear();
					depot.empty.add(list);
				}
			}
			
		}
		
		if(threads!=null){
			for(int i=1; i<threads.length; i++){
				while(threads[i]!=null && threads[i].getState()!=Thread.State.TERMINATED){
					try {
						threads[i].join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		assert(threads==null || threads.length<2 || threads[1]==null || !threads[1].isAlive()) : ((ReadThread)threads[1]).generatedLocal;
//		threads=null;
//		System.out.println("C");

		if(verbose){System.err.println("shutdown exited; errorState="+errorState);}
	}

	@Override
	public boolean paired() {
		return producer1.paired() || producer2!=null;
	}
	
//	public static ConcurrentReadStreamInterface getReadInputStream(long maxReads, boolean colorspace, String...args){
//		return getReadInputStream(maxReads, colorspace, false, args);
//	}

	private static ConcurrentReadStreamInterface getReadInputStream(long maxReads, boolean colorspace, boolean keepSamHeader, boolean allowSubprocess, String...args){
		assert(args.length>0) : Arrays.toString(args);
		for(int i=0; i<args.length; i++){
			if("null".equalsIgnoreCase(args[i])){args[i]=null;}
		}
		assert(args[0]!=null) : Arrays.toString(args);
		
		assert(args.length<2 || !args[0].equalsIgnoreCase(args[1]));
		String in1=args[0], in2=null, qf1=null, qf2=null;
		if(args.length>1){in2=args[1];}
		if(args.length>2){qf1=args[2];}
		if(args.length>3){qf2=args[3];}

		final FileFormat ff1=FileFormat.testInput(in1, null, allowSubprocess);
		final FileFormat ff2=FileFormat.testInput(in2, null, allowSubprocess);
		
		if(verbose){
			System.err.println("getReadInputStream("+maxReads+", "+colorspace+", "+keepSamHeader+", "+allowSubprocess+", "+in1+", "+in2+", "+qf1+", "+qf2+")");
		}
		
		return getReadInputStream(maxReads, colorspace, keepSamHeader, ff1, ff2, qf1, qf2);
	}
	
	public static ConcurrentReadStreamInterface getReadInputStream(long maxReads, boolean colorspace, boolean keepSamHeader, FileFormat ff1, FileFormat ff2){
		return getReadInputStream(maxReads, colorspace, keepSamHeader, ff1, ff2, (String)null, (String)null);
	}
	
	public static ConcurrentReadStreamInterface getReadInputStream(long maxReads, boolean colorspace, boolean keepSamHeader, FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		
		if(verbose){
			System.err.println("getReadInputStream("+maxReads+", "+colorspace+", "+keepSamHeader+", "+ff1+", "+ff2+", "+qf1+", "+qf2+")");
		}
		
		assert(ff1!=null);
		assert(ff2==null || ff1.name()==null || !ff1.name().equalsIgnoreCase(ff2.name()));
		assert(qf1==null || ff1.name()==null || !ff1.name().equalsIgnoreCase(qf2));
		assert(qf1==null || qf2==null || qf1.equalsIgnoreCase(qf2));
		
		final ConcurrentReadStreamInterface cris;
		
		if(ff1.fastq()){
			
			ReadInputStream ris1=new FastqReadInputStream(ff1, colorspace);
			ReadInputStream ris2=(ff2==null ? null : new FastqReadInputStream(ff2, colorspace));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);
			
		}else if(ff1.fasta()){
			
			ReadInputStream ris1=(qf1==null ? new FastaReadInputStream(ff1, colorspace, (FASTQ.FORCE_INTERLEAVED && ff2==null), ff2==null ? Shared.READ_BUFFER_MAX_DATA : -1)
				: new FastaQualReadInputStream(ff1, qf1, colorspace));
			ReadInputStream ris2=(ff2==null ? null : qf2==null ? new FastaReadInputStream(ff2, colorspace, false, -1) : new FastaQualReadInputStream(ff2, qf2, colorspace));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);
			
		}else if(ff1.scarf()){
			
			ReadInputStream ris1=new ScarfReadInputStream(ff1, colorspace);
			ReadInputStream ris2=(ff2==null ? null : new ScarfReadInputStream(ff2, colorspace));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);
			
		}else if(ff1.samOrBam()){
			
			ReadInputStream ris1=new SamReadInputStream(ff1, colorspace, keepSamHeader, FASTQ.FORCE_INTERLEAVED);
			ReadInputStream ris2=(ff2==null ? null : new SamReadInputStream(ff2, colorspace, false, false));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);
			
		}else if(ff1.bread()){
			
			RTextInputStream rtis=new RTextInputStream(ff1, ff2, maxReads);
			cris=new ConcurrentReadInputStream(rtis, maxReads); //TODO: Change to generic
			
			
		}else if(ff1.sequential()){
			SequentialReadInputStream ris=new SequentialReadInputStream(maxReads, 200, 50, 0, false);
			cris=new ConcurrentReadInputStream(ris, maxReads);
		}else if(ff1.csfasta()){
			
			if(ff2!=null){
				cris=new ConcurrentSolidInputStream(ff1, qf1, ff2, qf2, maxReads);
			}else{
				cris=new ConcurrentSolidInputStream(ff1, qf1, maxReads, null);
			}
		}else{
			cris=null;
			throw new RuntimeException(""+ff1);
		}
		
		return cris;
	}
	
	
	private class ReadThread extends Thread{
		ReadThread(ReadInputStream producer_, ArrayBlockingQueue<ArrayList<Read>> pq_){
			producer=producer_;
			pq=pq_;
		}
		
		@Override
		public void run(){
			readLists();
		}
		
		private final void readLists(){
			
			ArrayList<Read> list=null;
			
			if(verbose){System.err.println(getClass().getName()+" entering read lists loop.");}
			while(list!=null || (!shutdown && producer.hasMore() && generatedLocal<maxReads)){

				if(verbose){System.err.println(getClass().getName()+" looping: buffer1==null "+(list==null)+", shutdown="+shutdown+
						", producer.hasMore()="+producer.hasMore()+", generated<maxReads="+(generatedLocal<maxReads));}

				
				
				if(verbose){System.err.println(getClass().getName()+" Entering full fetch loop.");}
				while(generatedLocal<maxReads){
//					System.out.println("E");
					if(verbose){System.err.println(getClass().getName()+" depot.bufferSize="+depot.bufferSize+", generated="+generatedLocal);}
//					System.out.println("F");
					try {
						list=producer.nextList();
					} catch (Throwable e1) {
						// TODO
//						System.err.print('*');
						e1.printStackTrace();
						list=null;
						shutdown=true;
						try {
							pq.put(new ArrayList<Read>(1));
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						errorState=true;
					}
					if(verbose){System.err.println(getClass().getName()+" grabbed a list of size "+(list==null ? "null" : list.size()+""));}
//					System.out.println("G");
					if(list==null){
//						System.out.println("H");
						if(verbose){System.err.println(getClass().getName()+" broke loop on null list.");}
						break;
					}
					assert(list.size()>0);
					assert(list.size()<=BUF_LEN); //Although this is not really necessary.
//					System.out.println("I");
					if(list.size()+generatedLocal>maxReads){
//						System.out.println("J");
						if(verbose){System.err.println("Removing extra reads.");}
						while(list.size()+generatedLocal>maxReads){list.remove(list.size()-1);}
//						System.out.println("K");
					}
//					System.out.println("A");
					while(list!=null && !shutdown){
//						System.out.println("B");
						try {
							if(verbose){System.err.println("Trying to add list");}
							pq.put(list);
							generatedLocal+=list.size();
							list=null;
							if(verbose){
								System.out.println("Added list; pq.size() = "+pq.size());
							}
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
//						System.out.println("C");
					}
//					System.out.println("D");
					if(verbose){System.err.println("looping");}
				}

				if(verbose){System.err.println(getClass().getName()+" Finished inner loop iteration.\n");}
			}
			

			if(verbose){System.err.println(getClass().getName()+" attempting to poison output queue.");}
			boolean b=true;
			while(b){
				//TODO Note that this could cause a deadlock if there was a premature shutdown, so the consumer died while the queue was full.
				try {
//					pq.offer(poison, 10000, TimeUnit.SECONDS);
					pq.put(poison);
					b=false;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			

			if(verbose){System.err.println(getClass().getName()+" exited read lists loop: "+(list==null)+", "+shutdown+", "+producer.hasMore()+", "+generatedLocal+", "+maxReads);}

		}
		
		private final ArrayBlockingQueue<ArrayList<Read>> pq;
		private final ReadInputStream producer;
		private long generatedLocal=0;
	}
	
	private void incrementGenerated(long amt){
		generated+=amt;
		if(SHOW_PROGRESS && generated>=nextProgress){
			Data.sysout.print('.');
			nextProgress+=PROGRESS_INCR;
		}
	}
	
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		if(rate>=1f){
			randy=null;
		}else if(seed>-1){
			randy=new java.util.Random(seed);
		}else{
			randy=new java.util.Random();
		}
	}
	
	@Override
	public boolean errorState(){return errorState || 
			(producer1==null ? false : producer1.errorState()) || (producer2==null ? false : producer2.errorState());}
	/** TODO */
	private boolean errorState=false;
	
	private boolean[] running=new boolean[] {false};
	
	private float samplerate=1f;
	private java.util.Random randy=null;
	
	private ArrayBlockingQueue<ArrayList<Read>> p1q;
	private ArrayBlockingQueue<ArrayList<Read>> p2q;
	
	
	public Object[] producers(){return producer2==null ? new Object[] {producer1} : new Object[] {producer1, producer2};}

	private Thread[] threads;
	
	public final ReadInputStream producer1;
	public final ReadInputStream producer2;
	private ConcurrentDepot<Read> depot;
	
	private long maxReads;
	private long generated=0;
	private long listnum=0;
	private long nextProgress=PROGRESS_INCR;

	private final int BUF_LEN=Shared.READ_BUFFER_LENGTH;
	private final int NUM_BUFFS=Shared.READ_BUFFER_NUM_BUFFERS;
	private final long MAX_DATA=Shared.READ_BUFFER_MAX_DATA;
	
	public static boolean verbose=false;
	
	private static final ArrayList<Read> poison=new ArrayList<Read>(0);

	public static boolean SHOW_PROGRESS=false;
	public static long PROGRESS_INCR=1000000;
	public static boolean REMOVE_DISCARDED_READS=false;
	
}
