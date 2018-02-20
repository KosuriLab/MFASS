package stream;

import java.util.ArrayList;

import align2.Shared;

import dna.Data;
import fileIO.FileFormat;
import fileIO.TextFile;

public class FastaQualReadInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		
		FastaQualReadInputStream fris=new FastaQualReadInputStream(args[0], args[1], false);
		
		Read r=fris.next();
		int i=0;
		while(r!=null){
			System.out.println(r.toText(false));
			r=fris.next();
			if(i++>3){break;}
		}
		
	}
	
	public FastaQualReadInputStream(FileFormat ff, String qfname, boolean colorspace_){
//		assert(false) : "In progress";
		colorspace=colorspace_;

		tf=new TextFile(ff, false);
		qtf=new TextFile(qfname, false, false);
		interleaved=false;

	}
	
	public FastaQualReadInputStream(String fname, String qfname, boolean colorspace_){
//		assert(false) : "In progress";
		colorspace=colorspace_;
		
		if(!fileIO.FileFormat.hasFastaExtension(fname)){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+fname);
		}

		tf=new TextFile(fname, false, false);
		qtf=new TextFile(qfname, false, false);
		interleaved=false;

	}

	@Override
	public void start() {
//		if(cris!=null){new Thread(cris).start();}
	}
	
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.length){
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.length);
	}

	@Override
	public Read next() {
		if(!hasMore()){
			if(verbose){System.err.println("hasMore() returned false;  buffer="+(buffer==null ? null : buffer.length)+", next="+next+", consumed="+consumed);}
			return null;
		}
		Read r=buffer[next];
		buffer[next]=null;
		next++;
		consumed++;
		return r;
	}
	
	@Override
	public synchronized Read[] nextBlock() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.length){fillBuffer();}
		Read[] r=buffer;
		buffer=null;
		if(r!=null && r.length==0){r=null;}
		consumed+=(r==null ? 0 : r.length);
//		System.err.println(hashCode()+" produced "+r[0].numericID);
		return r;
	}
	
	@Override
	public synchronized ArrayList<Read> nextList() {
		return toList(nextBlock());
	}
	public final boolean preferArrays(){return true;}
	
	private synchronized void fillBuffer(){
		if(verbose){System.err.println("Filling buffer.  buffer="+(buffer==null ? null : buffer.length));}
		assert(buffer==null || next>=buffer.length);
		
		buffer=null;
		next=0;
		
		buffer=toReads(tf, BUF_LEN, nextReadID, interleaved, headerA);

		if(verbose){System.err.println("Filled buffer.  buffer="+(buffer==null ? null : buffer.length));}
		
		nextReadID+=buffer.length;
		if(buffer.length<BUF_LEN){
			if(verbose){System.err.println("Closing tf");}
			tf.close();
			qtf.close();
		}
		
		generated+=buffer.length;
		if(verbose){System.err.println("generated="+generated);}
		
		if(buffer.length>0 && readlen<0){readlen=buffer[0].bases.length;}
	}
	
	private Read[] toReads(TextFile tf, int maxReadsToReturn, long numericID, boolean interleaved, String[] headerA){
		ArrayList<Read> list=toReadList(tf, maxReadsToReturn, numericID, interleaved, headerA);
		assert(list.size()<=maxReadsToReturn);
		return list.toArray(new Read[list.size()]);
	}
	
	private ArrayList<Read> toReadList(TextFile tf, int maxReadsToReturn, long numericID, boolean interleaved, String[] headerA){
		if(finished){return null;}
		if(verbose){System.err.println("FastaRIS fetching a list.");}
		
		ArrayList<Read> list=new ArrayList<Read>(Data.min(16384, maxReadsToReturn));
		
		int added=0;
		
		Read prev=null;
		
		while(added<maxReadsToReturn){
			Read r=makeNextRead(tf, maxReadsToReturn, numericID, headerA);
			if(verbose){System.err.println("Made "+r);}
			if(r==null){
				if(verbose){System.err.println("makeNextRead returned null.");}
				break;
			}
			if(interleaved){
				if(prev==null){prev=r;}
				else{
					prev.mate=r;
					r.mate=prev;
					list.add(prev);
					added++;
					numericID++;
					prev=null;
				}
			}else{
				list.add(r);
				added++;
				numericID++;
			}
		}
		
		assert(list.size()<=maxReadsToReturn);
		if(verbose){System.err.println("FastaRIS returning a list.  Size="+list.size());}
		return list;
	}
	
	private Read makeNextRead(TextFile tf, int maxReadsToReturn, long numericID, String[] headerA){
		if(finished){
			if(verbose){System.err.println("Returning null because finished.");}
			return null;
		}
		assert(currentLine==null || currentLine.charAt(0)!='>');
		
		while(currentLine==null){
			currentLine=tf.nextLine();
			currentQLine=nextQtfLine(qtf);
			if(currentLine==null){
				if(verbose){System.err.println("Returning null because tf.nextLine()==null: A");}
				return null;
			}
			assert(currentLine.length()==currentQLine.length());
			if(currentLine.charAt(0)=='>'){
				assert(currentLine.equals(currentQLine));
				headerA[0]=currentLine;
				currentLoc=0;
				currentSection=0;
				currentLine=null;
			}
		}

		final boolean SPLIT_READS=FastaReadInputStream.SPLIT_READS;
		final int TARGET_READ_LEN=FastaReadInputStream.TARGET_READ_LEN;
		final int MIN_READ_LEN=FastaReadInputStream.MIN_READ_LEN;
		
		assert(currentLine==null || currentLine.charAt(0)!='>');

		StringBuilder sb=new StringBuilder();
		StringBuilder sbq=new StringBuilder();
		Read r=null;
		while(r==null){
			if(!SPLIT_READS || (currentLoc==0 && (currentLine.length()<=(TARGET_READ_LEN-sb.length())))){
				sb.append(currentLine);
				sbq.append(currentQLine);
				currentLoc=currentLine.length();
			}else{
				while(sb.length()<TARGET_READ_LEN && currentLoc<currentLine.length()){
					sb.append(currentLine.charAt(currentLoc));
					sbq.append(currentQLine.charAt(currentLoc));
					currentLoc++;
				}
			}
			assert(currentLine==null || currentLine.charAt(0)!='>');
			assert(sb.length()<=TARGET_READ_LEN) : sb.length()+"\n"+sb;
			
			if(sb.length()==TARGET_READ_LEN){
				assert(currentLine==null || currentLine.charAt(0)!='>');
				r=makeRead(sb, sbq, numericID);
				currentSection++;
				return r;
			}else{
				assert(currentLine==null || currentLine.charAt(0)!='>');
				assert(currentLoc>=currentLine.length()) : currentLoc+", "+currentLine.length()+", "+
						TARGET_READ_LEN+", "+sb.length()+"\n"+currentLine+"\n"+sb;
				currentLine=null;
				currentQLine=null;
				currentLoc=0;
				while(currentLine==null){
					currentLine=tf.nextLine();
					currentQLine=nextQtfLine(qtf);
					assert(currentLine==null || currentLine.length()==currentQLine.length());
					assert(currentLine==null || currentLine.charAt(0)!='>' || currentLine.equals(currentQLine));
					if(currentLine==null || currentLine.charAt(0)=='>'){
						if(sb.length()>=MIN_READ_LEN){
							r=makeRead(sb, sbq, numericID);
						}else{
							sb.setLength(0);
							sbq.setLength(0);
						}
						headerA[0]=currentLine;
						currentLoc=0;
						currentSection=0;
						currentLine=null;
						currentQLine=null;
						if(r!=null){return r;}
						if(headerA[0]==null){
							if(verbose){System.err.println("Returning null because tf.nextLine()==null: B");}
							return null;
						}
					}
					assert(currentLine==null || currentLine.charAt(0)!='>');
				}
				assert(currentLine==null || currentLine.charAt(0)!='>');
			}
			assert(currentLine==null || currentLine.charAt(0)!='>');
		}
		assert(currentLine==null || currentLine.charAt(0)!='>');
		if(verbose){System.err.println("Returning null because loop exited (should be unreachable).");}
		return null;
	}
	
	private final CharSequence nextQtfLine(TextFile qtf){
		String s=qtf.nextLine();
		if(!NUMERIC_QUAL || s==null || s.length()==0 || s.charAt(0)=='>'){return s;}
		s=s.trim();
		final StringBuilder sb=new StringBuilder((s.length()+1)/2);
		
		int x=0;
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c==' '){
				assert(i>0);
				sb.append((char)(x+FASTQ.ASCII_OFFSET));
				x=0;
			}else{
				x=10*x+(c-'0');
			}
		}
		sb.append((char)(x+FASTQ.ASCII_OFFSET));
		return sb;
	}
	
	private Read makeRead(StringBuilder sb, StringBuilder sbq, long numericID){
//		assert(!sb.equals(sbq)) : sb+"\n"+sbq;
		byte[] quals=new byte[sbq.length()];
		byte[] bases=new byte[sb.length()];
//		if(FAKE_QUALITY){
//			quals=new byte[sb.length()];
//			Arrays.fill(quals, (byte)(30));
//		}
		for(int i=0; i<bases.length; i++){
			bases[i]=(byte)Character.toUpperCase(sb.charAt(i));
		}
		for(int i=0; i<quals.length; i++){
			quals[i]=(byte)(sbq.charAt(i)-FASTQ.ASCII_OFFSET);
		}
		assert(bases[0]!='>') : new String(bases)+"\n"+numericID+"\n"+headerA[0];
		String hd=(currentSection>0 ? headerA[0].substring(1)+"_"+currentSection : new String(headerA[0].substring(1)));
		Read r=new Read(bases, (byte)0, (byte)0, 0, 0, hd, quals, colorspace, numericID);
		return r;
	}
	
	public boolean close(){
		boolean a=tf.close();
		boolean b=qtf.close();
		return a | b;
	}

	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		
		currentLine=null;
		currentLoc=0;
		currentSection=0;
		finished=false;

		tf.reset();
		qtf.reset();
	}

	@Override
	public boolean paired() {
		return interleaved;
	}

	private Read[] buffer=null;
	private int next=0;

	private final TextFile tf;
	private final TextFile qtf;
	private final boolean interleaved;

	public static final int BUF_LEN=Shared.READ_BUFFER_LENGTH;

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public int readlen=-1;
	
	public final boolean colorspace;
	private final String[] headerA=new String[1];
	
	public static boolean NUMERIC_QUAL=true;
	
	public static boolean verbose=false;
	public static boolean FAKE_QUALITY=false;

	private String currentLine=null;
	private CharSequence currentQLine=null;
	private int currentLoc=0;
	private int currentSection=0;
	private boolean finished=false;
	
}
