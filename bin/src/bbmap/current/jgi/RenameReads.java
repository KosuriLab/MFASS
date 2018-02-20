package jgi;

import java.util.ArrayList;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.Read;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;

import align2.ListNum;

/**
 * @author Brian Bushnell
 * @date Aug 23, 2013
 *
 */
public class RenameReads {
	
	public static void main(String[] args){
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=8;
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		
		long maxReads=-1;
		final ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(args[0], FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, null);
			Thread th=new Thread(cris);
			th.start();
		}
		
		TextStreamWriter tsw=new TextStreamWriter(args[2], false, false, true);
		tsw.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		long x=0;
		while(reads!=null && reads.size()>0){

			for(Read r : reads){
				r.id=args[1]+"_"+x;
				if(r.mate!=null){
					r.id=r.id+" /1";
					r.mate.id=args[1]+"_"+x+" /2";
				}
				tsw.println(r);
				if(r.mate!=null){tsw.println(r.mate);}
				x++;
			}

			cris.returnList(ln, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln, ln.list.isEmpty());
		ReadWrite.closeStream(cris);
		
		tsw.poisonAndWait();
		
	}
	
}
