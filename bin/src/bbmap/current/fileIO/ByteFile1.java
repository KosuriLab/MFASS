package fileIO;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import dna.Data;
import dna.Timer;


/**
 * @author Brian Bushnell
 *
 */
public class ByteFile1 extends ByteFile {
	
	
	public static void main(String[] args) throws IOException{
		ByteFile1 tf=new ByteFile1(args.length>0 ? args[0] : "stdin", false, true);
		long first=0, last=100;
		boolean speedtest=false;
		if(args.length>1){
			if(args[1].equalsIgnoreCase("speedtest")){
				speedtest=true;
				first=0;
				last=Long.MAX_VALUE;
			}else{
				first=Integer.parseInt(args[1]);
				last=first+100;
			}
		}
		if(args.length>2){
			last=Integer.parseInt(args[2]);
		}
		speedtest(tf, first, last, !speedtest);
		
		tf.close();
		tf.reset();
		tf.close();
	}
	
	private static void speedtest(ByteFile1 tf, long first, long last, boolean reprint){
		Timer t=new Timer();
		t.start();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){tf.readLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=tf.readLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length;
				System.out.println(new String(s));
			}
			
			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
			for(long i=first; i<last; i++){
				byte[] s=tf.readLine();
				if(s==null){break;}
				lines++;
				bytes+=s.length;
			}
		}
		t.stop();
		
		if(!reprint){
			double rpnano=lines/(double)(t.elapsed);
			double bpnano=bytes/(double)(t.elapsed);

			String rpstring=(lines<100000 ? ""+lines : lines<100000000 ? (lines/1000)+"k" : (lines/1000000)+"m");
			String bpstring=(bytes<100000 ? ""+bytes : bytes<100000000 ? (bytes/1000)+"k" : (bytes/1000000)+"m");

			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}

			System.err.println("Time:                         \t"+t);
			System.err.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk lines/sec", rpnano*1000000));
			System.err.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bytes/sec", bpnano*1000));
		}
	}
	
//	public ByteFile1(String name){this(name, false);}
	
	public ByteFile1(String fname, boolean tryAllExtensions, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false), tryAllExtensions);
	}
	
	public ByteFile1(FileFormat ff, boolean tryAllExtensions){
		super(ff, tryAllExtensions);
		if(verbose){System.err.println("ByteFile1("+ff+", "+tryAllExtensions+")");}
		is=open();
	}
	
	public final void reset(){
		close();
		is=open();
	}
	
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		if(!open){return errorState;}
		open=false;
		assert(is!=null);
		
		errorState|=ReadWrite.finishReading(is, name(), allowSubprocess());
		
		is=null;
		lineNum=-1;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		return errorState;
	}
	
	public byte[] nextLine(){
//		throw new RuntimeException("Please implement "+getClass()+".nextLine()");
		return readLine();
	}
	
	public final byte[] readLine(){
		
		if(!open || is==null){
			if(Data.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

//		System.out.println("\nCalled readLine() for line "+lineNum);
//		System.out.println("A: bstart="+bstart+", bstop="+bstop);
		
		if(bstart<bstop && lasteol==slashr && buffer[bstart]==slashn){bstart++;}
		assert(bstart>=bstop || (buffer[bstart]!=slashr || buffer[bstart]!=slashn)/*buffer[bstart]>slashr || buffer[bstart]==slashn*/);
		int nlpos=bstart;
		
//		System.out.println("B: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
//		while(nlpos<bstop && (buffer[nlpos]>slashr || buffer[nlpos]==tab)){nlpos++;}
		while(nlpos<bstop && (buffer[nlpos]!=slashr && buffer[nlpos]!=slashn)){nlpos++;}
//		System.out.println("C: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
		if(nlpos>=bstop){
			nlpos=fillBuffer();
//			System.out.println("Filled buffer.");
		}else{
			lasteol=buffer[nlpos];
		}
//		System.out.println("D: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
		
		if(nlpos<0 || bstop<1){
			close();
			return null;
		}

		lineNum++;
		if(bstart==nlpos){//Empty line.
			bstart=nlpos+1;
			return null;
		}
		byte[] line=Arrays.copyOfRange(buffer, bstart, nlpos);
		assert(line.length>0) : bstart+", "+nlpos;
		bstart=nlpos+1;
//		System.out.println("E: bstart="+bstart+", bstop="+bstop+", nlpos="+nlpos);
		return line;
	}
	
	private int fillBuffer(){
		if(bstart<bstop){ //Shift end bytes to beginning
			assert(bstart>0);
//			assert(bstop==buffer.length);
			int extra=bstop-bstart;
			for(int i=0; i<extra; i++, bstart++){
				buffer[i]=buffer[bstart];
//				assert(buffer[i]>=slashr || buffer[i]==tab);
				assert(buffer[i]!=slashr && buffer[i]!=slashn);
			}
			bstop=extra;
		}else{
			bstop=0;
		}

		bstart=0;
		int len=bstop;
		int r=-1;
		while(len==bstop){//hit end of input without encountering a newline
			if(bstop==buffer.length){
				buffer=Arrays.copyOf(buffer, buffer.length*2);
			}
			try {
				r=is.read(buffer, bstop, buffer.length-bstop);
			} catch (IOException e) {
				e.printStackTrace();
				System.err.println("open="+open);
			}
			if(r>0){
				bstop=bstop+r;
//				while(len<bstop && (buffer[len]>slashr || buffer[len]==tab)){len++;}
				while(len<bstop && (buffer[len]!=slashr && buffer[len]!=slashn)){len++;}
			}else{
				len=bstop;
				break;
			}
		}
		
//		System.out.println("Filled buffer; r="+r+", returning "+len);
		assert(r==-1 || buffer[len]<=slashr);
		if(len>0){lasteol=buffer[len];}
		return len;
	}
	
	private final synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened TextFile "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), false, allowSubprocess());
		bstart=-1;
		bstop=-1;
		lasteol=-1;
		return is;
	}
	
	public boolean isOpen(){return open;}
	
	public final InputStream is(){return is;}
	
	public final long lineNum(){return lineNum;}
	
	private boolean open=false;
	private byte[] buffer=new byte[16384];
	private int bstart=0, bstop=0;
	public InputStream is;
	public long lineNum=-1;
	
	private byte lasteol=-1;
	
	public static boolean verbose=false;

	private boolean errorState=false;
	
}
