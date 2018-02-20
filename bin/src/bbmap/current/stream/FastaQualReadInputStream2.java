package stream;

import java.util.ArrayList;

import align2.Shared;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;

public class FastaQualReadInputStream2 extends ReadInputStream {
	
	public static void main(String[] args){
		
		FastaQualReadInputStream2 fris=new FastaQualReadInputStream2(args[0], args[1], false, true);
		
		Read r=fris.next();
		int i=0;
		while(r!=null){
			System.out.println(r.toText(false));
			r=fris.next();
			if(i++>3){break;}
		}
		
	}
	
	public FastaQualReadInputStream2(String fname, String qfname, boolean colorspace_, boolean allowSubprocess_){
//		assert(false) : "In progress";
		colorspace=colorspace_;
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false);
		if(!ff.fasta() && !ff.stdio()){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+fname);
		}

		tf=ByteFile.makeByteFile(ff, false);
		qtf=ByteFile.makeByteFile(FileFormat.testInput(qfname, FileFormat.QUAL, null, allowSubprocess_, false), false);
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
	}
	
	private Read[] toReads(ByteFile tf, int maxReadsToReturn, long numericID, boolean interleaved, String[] headerA){
		ArrayList<Read> list=toReadList(tf, maxReadsToReturn, numericID, interleaved, headerA);
		assert(list.size()<=maxReadsToReturn);
		return list.toArray(new Read[list.size()]);
	}
	
	private ArrayList<Read> toReadList(ByteFile tf, int maxReadsToReturn, long numericID, boolean interleaved, String[] headerA){
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
	
	private Read makeNextRead(ByteFile tf, int maxReadsToReturn, long numericID, String[] headerA){
		if(finished){
			if(verbose){System.err.println("Returning null because finished.");}
			return null;
		}
		assert(currentLine==null || currentLine[0]!=carrot);
		
		while(currentLine==null){
			currentLine=tf.nextLine();
			currentQLine=nextQtfLine(qtf);
			if(currentLine==null){
				if(verbose){System.err.println("Returning null because tf.nextLine()==null: A");}
				return null;
			}
			assert(currentLine.length==currentQLine.length);
			if(currentLine[0]==carrot){
				assert(currentLine.equals(currentQLine));
				headerA[0]=new String(currentLine);
				currentLoc=0;
				currentSection=0;
				currentLine=null;
			}
		}

		final boolean SPLIT_READS=FastaReadInputStream.SPLIT_READS;
		final int TARGET_READ_LEN=FastaReadInputStream.TARGET_READ_LEN;
		final int MIN_READ_LEN=FastaReadInputStream.MIN_READ_LEN;
		
		assert(currentLine==null || currentLine[0]!=carrot);

		StringBuilder sb=new StringBuilder();
		StringBuilder sbq=new StringBuilder();
		Read r=null;
		while(r==null){
			if(!SPLIT_READS || (currentLoc==0 && (currentLine.length<=(TARGET_READ_LEN-sb.length())))){
//				sb.append(currentLine);
//				sbq.append(currentQLine);
				for(byte b : currentLine){sb.append((char)b);}
				for(byte b : currentQLine){sbq.append((char)b);}
				currentLoc=currentLine.length;
			}else{
				while(sb.length()<TARGET_READ_LEN && currentLoc<currentLine.length){
					sb.append((char)currentLine[currentLoc]);
					sbq.append((char)currentQLine[currentLoc]);
					currentLoc++;
				}
			}
			assert(currentLine==null || currentLine[0]!=carrot);
			assert(sb.length()<=TARGET_READ_LEN) : sb.length()+"\n"+sb;
			
			if(sb.length()==TARGET_READ_LEN){
				assert(currentLine==null || currentLine[0]!=carrot);
				r=makeRead(sb, sbq, numericID);
				currentSection++;
				return r;
			}else{
				assert(currentLine==null || currentLine[0]!=carrot);
				assert(currentLoc>=currentLine.length) : currentLoc+", "+currentLine.length+", "+
						TARGET_READ_LEN+", "+sb.length()+"\n"+currentLine+"\n"+sb;
				currentLine=null;
				currentQLine=null;
				currentLoc=0;
				while(currentLine==null){
					currentLine=tf.nextLine();
					currentQLine=nextQtfLine(qtf);
					assert(currentLine==null || currentLine.length==currentQLine.length);
					assert(currentLine==null || currentLine[0]!=carrot || currentLine.equals(currentQLine));
					if(currentLine==null || currentLine[0]==carrot){
						if(sb.length()>=MIN_READ_LEN){
							r=makeRead(sb, sbq, numericID);
						}else{
							sb.setLength(0);
							sbq.setLength(0);
						}
						headerA[0]=currentLine==null ? null : new String(currentLine);
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
					assert(currentLine==null || currentLine[0]!=carrot);
				}
				assert(currentLine==null || currentLine[0]!=carrot);
			}
			assert(currentLine==null || currentLine[0]!=carrot);
		}
		assert(currentLine==null || currentLine[0]!=carrot);
		if(verbose){System.err.println("Returning null because loop exited (should be unreachable).");}
		return null;
	}
	
	private final byte[] nextQtfLine(ByteFile qtf){
		byte[] s=qtf.nextLine();
		if(!NUMERIC_QUAL || s==null || s.length==0 || s[0]==carrot){return s;}
		
		int last=s.length-1;
		while(last>=0 && Character.isWhitespace(s[last])){
			last--;
		}
		if(last<0){return new byte[0];}
		final int lim=last+1;
		int spaces=0;
		for(int i=0; i<lim; i++){
			if(s[i]==space){spaces++;}
		}
		
		final byte[] sb=new byte[s.length==0 ? 0 : spaces+1];
		
		int x=0, j=0;
		for(int i=0; i<lim; i++){
			byte b=s[i];
			if(b==space){
				assert(i>0);
				sb[j]=(byte)(x+FASTQ.ASCII_OFFSET);
				x=0;
				j++;
			}else{
				x=10*x+(b-zero);
			}
		}
		sb[j]=(byte)(x+FASTQ.ASCII_OFFSET);
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
		assert(bases[0]!=carrot) : new String(bases)+"\n"+numericID+"\n"+headerA[0];
		String hd=(currentSection>0 ? headerA[0].substring(1)+"_"+currentSection : new String(headerA[0].substring(1)));
		Read r=new Read(bases, (byte)0, (byte)0, 0, 0, hd, quals, colorspace, numericID);
		return r;
	}
	
	public boolean close(){
		boolean a=tf.close();
		boolean b=qtf.close();
		return a|b;
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

	private final ByteFile tf;
	private final ByteFile qtf;
	private final boolean interleaved;

	public static final int BUF_LEN=Shared.READ_BUFFER_LENGTH;

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public final boolean colorspace;
	private final String[] headerA=new String[1];
	
	public static boolean NUMERIC_QUAL=true;
	
	public static boolean verbose=false;
	public static boolean FAKE_QUALITY=false;

	private byte[] currentLine=null;
	private byte[] currentQLine=null;
//	private CharSequence currentQLine=null;
	private int currentLoc=0;
	private int currentSection=0;
	private boolean finished=false;
	private final byte carrot='>', space=' ', zero='0';
	
}
