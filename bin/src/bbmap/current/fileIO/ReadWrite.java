package fileIO;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.Reader;
import java.lang.ProcessBuilder.Redirect;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import stream.ConcurrentReadStreamInterface;
import stream.RTextOutputStream3;

import align2.Shared;
import align2.Tools;

import dna.Data;

public class ReadWrite {
	
	
	public static void main(String[] args){
		File f=new File(args[1]);
		assert(!f.exists()) : "Destination file already exists.";
		copyFile(args[0], args[1]);
	}
	
	public static void writeStringInThread(CharSequence x, String fname){
		writeStringInThread(x, fname, false);
	}
	
	public static void writeStringInThread(CharSequence x, String fname, boolean append){
		addThread(1);
		new Thread(new WriteStringThread(x, fname, append)).start();
	}
	
	public static void writeObjectInThread(Object x, String fname, boolean allowSubprocess){
		addThread(1);
		new Thread(new WriteObjectThread(x, fname, allowSubprocess)).start();
	}
	
	private static class WriteStringThread implements Runnable{
		
		private final CharSequence x;
		private final String fname;
		private final boolean append;
		WriteStringThread(CharSequence x_, String fname_, boolean append_){
			x=x_;
			fname=fname_;
			append=append_;
		}
		
		@Override
		public void run() {
			if(verbose){System.err.println("WriteStringThread.run() started for fname "+fname);}
			addRunningThread(1);
			writeStringAsync(x, fname, append);
			addThread(-1);
			if(verbose){System.err.println("WriteStringThread.run() finished for fname "+fname);}
		}
		
	}
	
	private static class WriteObjectThread implements Runnable{
		
		private final Object x;
		private final String fname;
		private final boolean allowSubprocess;
		WriteObjectThread(Object x_, String fname_, boolean allowSubprocess_){
			x=x_;
			fname=fname_;
			allowSubprocess=allowSubprocess_;
		}
		
		@Override
		public void run() {
			if(verbose){System.err.println("WriteObjectThread.run() started for fname "+fname);}
			addRunningThread(1);
//			System.out.println(fname+" began writing.");
			writeAsync(x, fname, allowSubprocess);
//			System.out.println(fname+" finished writing.");
			addThread(-1);
//			System.out.println(fname+" reports "+countActiveThreads()+" active threads.");
			if(verbose){System.err.println("WriteObjectThread.run() finished for fname "+fname);}
		}
		
	}
	
	public static boolean setPermissions(String fname, boolean read, boolean write, boolean execute, boolean ownerOnly){
		File f=new File(fname);
		if(!f.exists()){return false;}
		try {
			f.setReadable(read, ownerOnly);
			f.setWritable(write, ownerOnly);
			f.setExecutable(execute, ownerOnly);
		} catch (Exception e) {
			return false;
		}
		return true;
	}

	public static void writeString(CharSequence x, String fname){writeString(x, fname, false);}
	public static void writeString(CharSequence x, String fname, boolean append){
		if(verbose){System.err.println("writeString(x, "+fname+", "+append+")");}
		OutputStream os=getOutputStream(fname, append, true, false);
		
		try {

			synchronized(diskSync){
				PrintWriter out=new PrintWriter(os);
				out.print(x);
				out.flush(); 

				if(os.getClass()==ZipOutputStream.class){
					ZipOutputStream zos=(ZipOutputStream)os;
					zos.closeEntry();
					zos.finish();
				}
//				else if(PROCESS_XZ && os.getClass()==org.tukaani.xz.XZOutputStream.class){
//					org.tukaani.xz.XZOutputStream zos=(org.tukaani.xz.XZOutputStream)os;
//					zos.finish();
//				}
				out.close();
			}
//			System.out.println("Wrote to "+fname);
			
//			String read=readString(fname);
//			assert(x.equals(read)) : x.length()+", "+read.length();
			
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public static void writeStringAsync(CharSequence x, String fname){writeStringAsync(x, fname, false);}
	public static void writeStringAsync(CharSequence x, String fname, boolean append){
		if(verbose){System.err.println("writeStringAsync(x, "+fname+", "+append+")");}
		
		OutputStream os=getOutputStream(fname, append, true, false);
		
		try {

			synchronized(diskSync){
				PrintWriter out=new PrintWriter(os);
				out.print(x);
				out.flush(); 

				if(os.getClass()==ZipOutputStream.class){
					ZipOutputStream zos=(ZipOutputStream)os;
					zos.closeEntry();
					zos.finish();
				}
//				else if(PROCESS_XZ && os.getClass()==org.tukaani.xz.XZOutputStream.class){
//					org.tukaani.xz.XZOutputStream zos=(org.tukaani.xz.XZOutputStream)os;
//					zos.finish();
//				}
				out.close();
			}
//			System.out.println("Wrote to "+fname);
			
//			String read=readString(fname);
//			assert(x.equals(read)) : x.length()+", "+read.length();
			
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	public static <X> void write(X x, String fname, boolean allowSubprocess){
		if(verbose){System.err.println("write(x, "+fname+", "+allowSubprocess+")");}
		
		OutputStream os=getOutputStream(fname, false, true, allowSubprocess);
		
		try {

			synchronized(diskSync){
				ObjectOutputStream out=new ObjectOutputStream(os);
				out.writeObject(x);
				close(out);
			}
			
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	public static <X> void writeAsync(X x, String fname, boolean allowSubprocess){
		if(verbose){System.err.println("writeAsync(x, "+fname+", "+allowSubprocess+")");}
		
		OutputStream os=getOutputStream(fname, false, true, allowSubprocess);
		
		try {

			ObjectOutputStream out=new ObjectOutputStream(os);
			out.writeObject(x);
			close(out);

		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	public static final boolean finishReading(InputStream is, String fname, boolean killProcess, Reader...ra){
		if(verbose){System.err.println("finishReading("+is+", "+fname+", "+killProcess+", "+ra.length+")");}
		boolean error=false;
		if(ra!=null){
			for(Reader r : ra){
				try {
					r.close();
				} catch (IOException e) {
					error=true;
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		error|=finishReading(is, fname, killProcess);
		if(verbose){System.err.println("finishReading("+is+", "+fname+", "+killProcess+", "+ra.length+") returned "+error);}
		return error;
	}
	
	public static final boolean finishReading(InputStream is, String fname, boolean killProcess){
		if(verbose){System.err.println("finishReading("+is+", "+fname+", "+killProcess+")");}
		boolean error=false;
		if(is!=System.in){
			try {
				is.close();
			} catch (IOException e) {
				error=true;
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(killProcess && fname!=null && is!=System.in){error|=ReadWrite.killProcess(fname);}
		if(verbose){System.err.println("finishReading("+is+", "+fname+", "+killProcess+") returned "+error);}
		return error;
	}
	
//	public static final boolean finishWriting(PrintWriter writer, OutputStream outStream, String fname){
//		return finishWriting(writer, outStream, fname, fname!=null);
//	}
	
	public static final boolean finishWriting(PrintWriter writer, OutputStream outStream, String fname, boolean killProcess){
		if(verbose){System.err.println("finishWriting("+writer+", "+outStream+" , "+fname+", "+killProcess+")");}
		boolean error=false;
		if(writer!=null){writer.flush();}
		if(outStream!=System.out && outStream!=System.err){
			close(outStream);
		}
		if(writer!=null){writer.close();}
		if(killProcess && fname!=null && outStream!=System.err && outStream!=System.out){error|=ReadWrite.killProcess(fname);}
		if(verbose){System.err.println("finishWriting("+writer+", "+outStream+" , "+fname+", "+killProcess+") returned "+error);}
		return error;
	}
	
	public static final boolean close(OutputStream os, String fname){
		if(verbose){System.err.println("close("+os+", "+fname+")");}
		boolean error=false;
		if(os!=null){error|=close(os);}
		if(fname!=null && os!=System.err && os!=System.out){error|=killProcess(fname);}
		if(verbose){System.err.println("close("+os+", "+fname+") returned "+error);}
		return error;
	}
	
	public static final boolean close(OutputStream os){
		if(verbose){System.err.println("close("+os+")");}
		boolean error=false;
		try {
			os.flush();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			error=true;
		}
		if(os.getClass()==ZipOutputStream.class){
			ZipOutputStream zos=(ZipOutputStream)os;
			try {
				zos.closeEntry();
				zos.finish();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				error=true;
			}
		}
//		else if(PROCESS_XZ && os.getClass()==org.tukaani.xz.XZOutputStream.class){
//			org.tukaani.xz.XZOutputStream zos=(org.tukaani.xz.XZOutputStream)os;
//			try {
//				zos.finish();
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
		try {
			os.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			error=true;
		}
		if(verbose){System.err.println("close("+os+") returned "+error);}
		return error;
	}
	
	
	
	@Deprecated
	public static OutputStream getOutputStream_old(String fname, boolean append){

//		fname=fname.replaceAll("\\\\", "/");
		fname=fname.replace('\\', '/');
		assert(fname.indexOf('\\')<0);
//		assert(!fname.contains("//"));

		boolean xz=fname.endsWith(".xz");
		boolean gzipped=fname.endsWith(".gz");
		boolean zipped=fname.endsWith(".zip");
		boolean bzipped=PROCESS_BZ2 && fname.endsWith(".bz2");
		final String basename=basename(fname);
		
		OutputStream out=null;
		
		try {
			if(fname.equals("stdout") || fname.startsWith("stdout.")){
				out=System.out;
			}else{
				FileOutputStream fos=new FileOutputStream(fname, append);
				out=new BufferedOutputStream(fos);
			}
			if(RAWMODE){return out;}
			
			if(zipped){
				ZipOutputStream zos=new ZipOutputStream(out);
				zos.setLevel(ZIPLEVEL);
				zos.putNextEntry(new ZipEntry(basename));
				out=zos;
			}else if(gzipped){
				GZIPOutputStream gos=new GZIPOutputStream(out, 8192){
					{
						//					        def.setLevel(Deflater.DEFAULT_COMPRESSION);
						def.setLevel(ZIPLEVEL);
					}
				};

				out=gos;
			}else if(bzipped){
				throw new RuntimeException("bz2 compression not supported in public version.");
//				out.write('B');
//				out.write('Z');
//				CBZip2OutputStream zos=new CBZip2OutputStream(out, 8192);
//				out=zos;
			}
			//				else if(PROCESS_XZ && xz){
			//					org.tukaani.xz.LZMA2Options options = new org.tukaani.xz.LZMA2Options();
			//					options.setPreset(ZIPLEVEL);
			//					org.tukaani.xz.XZOutputStream zos=new org.tukaani.xz.XZOutputStream(out, options);
			//					out=zos;
			//				}

			
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		assert(out!=null);
		return out;
	}
	

	public static OutputStream getOutputStream(String fname, boolean append, boolean buffered, boolean allowSubprocess){
		
		if(verbose){
			System.err.println("getOutputStream("+fname+", "+append+", "+buffered+", "+allowSubprocess+")");
			new Exception().printStackTrace(System.err);
		}
		
//		assert(false) : fname; //TODO: for testing
//		fname=fname.replaceAll("\\\\", "/");
		fname=fname.replace('\\', '/');
		assert(fname.indexOf('\\')<0);
//		assert(!fname.contains("//"));
		
		boolean gzipped=fname.endsWith(".gz") || fname.endsWith(".gzip");
		boolean zipped=fname.endsWith(".zip");
		boolean bzipped=PROCESS_BZ2 && fname.endsWith(".bz2");
		boolean xz=PROCESS_XZ && fname.endsWith(".xz");
		
//		assert(false) : fname;
		
		allowSubprocess=(allowSubprocess && Shared.THREADS>1);
				
		if(gzipped){
			assert(!append);
			return getGZipOutputStream(fname, allowSubprocess);
		}else if(zipped){
			assert(!append);
			return getZipOutputStream(fname, buffered, allowSubprocess);
		}else if(bzipped){
			assert(!append);
			return getBZipOutputStream(fname, buffered, allowSubprocess);
		}else if(xz){
			assert(!append);
			return getXZOutputStream(fname, buffered, allowSubprocess);
		}
		return getRawOutputStream(fname, append, buffered);
	}
	
	public static OutputStream getRawOutputStream(String fname, boolean append, boolean buffered){
		
		if(verbose){System.err.println("getRawOutputStream("+fname+", "+append+", "+buffered+")");}
		
		if(fname.equals("stdout") || fname.startsWith("stdout.")){
			return System.out;
		}else if(fname.equals("stderr") || fname.startsWith("stderr.")){
			return System.err;
		}
		FileOutputStream fos=null;
		try {
			fos = new FileOutputStream(fname, append);
		} catch (FileNotFoundException e) {
			synchronized(ReadWrite.class){
				try {
					File f=new File(fname);
					String parent=f.getParent();
					f=new File(parent);
					if(!f.exists()){f.mkdirs();}
					fos = new FileOutputStream(fname, append);
				} catch (Exception e2) {
					throw new RuntimeException(e2);
				}
			}
		}
		assert(fos!=null);
		if(buffered){return new BufferedOutputStream(fos);}
		return fos;
	}
	
	public static OutputStream getXZOutputStream(String fname, boolean buffered, boolean allowSubprocess){
		final OutputStream raw=getRawOutputStream(fname, false, buffered);
		if(RAWMODE){return raw;}
		throw new RuntimeException("Unsupported format: XZ");
//		try {
//			org.tukaani.xz.LZMA2Options options = new org.tukaani.xz.LZMA2Options();
//			options.setPreset(ZIPLEVEL);
//			org.tukaani.xz.XZOutputStream out=new org.tukaani.xz.XZOutputStream(raw, options);
//			return out;
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		assert(false);
//		return null;
	}
	
	public static OutputStream getBZipOutputStream(String fname, boolean buffered, boolean allowSubprocess){
		if(verbose){System.err.println("getBZipOutputStream("+fname+", "+buffered+", "+allowSubprocess+")");}
		final OutputStream raw=getRawOutputStream(fname, false, RAWMODE);
		if(RAWMODE){return raw;}
		throw new RuntimeException("bz2 compression not supported in public version.");
//		try {
//			raw.write('B');
//			raw.write('Z');
//			CBZip2OutputStream out=new CBZip2OutputStream(raw, 8192);
//			return out;
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		assert(false);
//		return null;
	}
	
	public static OutputStream getZipOutputStream(String fname, boolean buffered, boolean allowSubprocess){
		if(verbose){System.err.println("getZipOutputStream("+fname+", "+buffered+", "+allowSubprocess+")");}
		final OutputStream raw=getRawOutputStream(fname, false, buffered);
		if(RAWMODE){return raw;}
		try {
			ZipOutputStream out=new ZipOutputStream(raw);
			out.setLevel(ZIPLEVEL);
			final String basename=basename(fname);
			out.putNextEntry(new ZipEntry(basename));
			return out;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert(false);
		return null;
	}
	
	public static OutputStream getGZipOutputStream(String fname, boolean allowSubprocess){
		if(verbose){System.err.println("getGZipOutputStream("+fname+", "+allowSubprocess+")");}
		if(allowSubprocess && Shared.THREADS>2){
			if(USE_PIGZ && Data.PIGZ() && (Data.SH() /*|| fname.equals("stdout") || fname.startsWith("stdout.")*/)){return getPigzStream(fname);}
			if(USE_GZIP && Data.GZIP() && (Data.SH() /*|| fname.equals("stdout") || fname.startsWith("stdout.")*/)){return getGzipStream(fname);}
		}
		
		final OutputStream raw=getRawOutputStream(fname, false, false);
		if(RAWMODE){return raw;}
		try {
			final GZIPOutputStream out=new GZIPOutputStream(raw, 8192){
				{
					//					        def.setLevel(Deflater.DEFAULT_COMPRESSION);
					def.setLevel(ZIPLEVEL);
				}
			};
			return out;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert(false);
		return null;
	}
	
	public static OutputStream getPigzStream(String fname){
		if(verbose){System.err.println("getPigzStream("+fname+")");}
		int threads=Tools.min(MAX_ZIP_THREADS, Tools.max(Shared.THREADS/Tools.max(ZIP_THREAD_DIVISOR, 1), 1));
		threads=Tools.max(1, threads);
		int zl=ZIPLEVEL;
		if(threads>=4 && zl>0 && zl<4){zl=4;}
		OutputStream out=getOutputStreamFromProcess(fname, "pigz -c -p "+threads+" -"+zl, true);
		return out;
	}
	
	public static OutputStream getGzipStream(String fname){
		if(verbose){System.err.println("getGzipStream("+fname+")");}
		OutputStream out=getOutputStreamFromProcess(fname, "gzip -c -"+ZIPLEVEL, true);
		return out;
	}
	
	public static OutputStream getOutputStreamFromProcess(String fname, String command, boolean sh){
		if(verbose){System.err.println("getOutputStreamFromProcess("+fname+", "+command+", "+sh+")");}
		
		OutputStream out=null;
		Process p=null;
		
		boolean useProcessBuilder=false;
		
		if(useProcessBuilder){
			ProcessBuilder pb=new ProcessBuilder();
			pb.redirectError(Redirect.INHERIT);
			
			if(fname.equals("stdout") || fname.startsWith("stdout.")){
				pb.redirectOutput(Redirect.INHERIT);
				pb.command(command.split(" "));
			}else{
				
				if(fname!=null){
					pb.redirectOutput(new File(fname));
				}
				
				pb.command(command.split(" "));
				
//				if(sh){
//					pb.command(("sh -c "+command+" 1>"+fname).split(" "));
//				}else{
//					pb.command(command.split(" "));
//				}
			}
			try {
				p=pb.start();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			addProcess(fname, p);
			out=p.getOutputStream();
			return out;
		}
		
		if(fname.equals("stdout") || fname.startsWith("stdout.")){
			try {
				p = Runtime.getRuntime().exec(command);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			InputStream is=p.getInputStream();
			PipeThread it=new PipeThread(is, System.out);
			addPipeThread(fname, it);
			it.start();
//		}else if(fname.equals("stderr") || fname.startsWith("stderr.")){
//			try {
//				p = Runtime.getRuntime().exec(command);
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//			InputStream is=p.getErrorStream();
//			PipeThread it=new PipeThread(is, System.err);
//			it.start();
		}else{
			try {
				if(sh){
					String[] cmd = {
							"sh",
							"-c",
							command+" 1>"+fname
					};
					p=Runtime.getRuntime().exec(cmd);
				}else{
					p=Runtime.getRuntime().exec(command);
				}
//				p = Runtime.getRuntime().exec("gzip -c -"+ZIPLEVEL+" 1>"+fname);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		addProcess(fname, p);
		out=p.getOutputStream();
		InputStream es=p.getErrorStream();
		assert(es!=null);
		PipeThread et=new PipeThread(es, System.err);
		addPipeThread(fname, et);
		et.start();

		return out;
	}
	
	public static String readString(String fname){
		if(verbose){System.err.println("readString("+fname+")");}
		String x=null;
		InputStream is=getInputStream(fname, false, false);
		
		try {
			
			StringBuilder sb=new StringBuilder();
			
//			synchronized(diskSync){
				BufferedReader in=new BufferedReader(new InputStreamReader(is), INBUF);
				String temp=in.readLine();
				while(temp!=null){
					sb.append(temp).append('\n');
					temp=in.readLine();
				}
				in.close();
//			}
			
			x=sb.toString();
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		
		return x;
	}
	
	public static Object readObject(String fname){
		if(verbose){System.err.println("readObject("+fname+")");}
		Object x=null;
		InputStream is=getInputStream(fname, true, false);
		
		try {
//			synchronized(diskSync){
				ObjectInputStream in=new ObjectInputStream(is);
				x=in.readObject();
				in.close();
//			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		}
		
		return x;
	}
	
//	public static InputStream getInputStream(String fname){
//		
//		boolean gzipped=fname.endsWith(".gz");
//		boolean zipped=fname.endsWith(".zip");
//		boolean bzipped=false;//fname.endsWith(".bz2");
//		final String basename=basename(fname);
//		
//		InputStream in=null;
//		try {
//
//			FileInputStream fis=new FileInputStream(fname);
//			BufferedInputStream bis=new BufferedInputStream(fis, INBUF);
//			
//			if(zipped && !RAWMODE){
//				ZipInputStream zis=new ZipInputStream(bis);
//				ZipEntry ze=zis.getNextEntry();
//				assert(ze!=null);
//				assert(basename.equals(ze.getName())) : basename+" != "+ze.getName();
//				in=zis;
//			}else if(gzipped && !RAWMODE){
//				in=new GZIPInputStream(bis, 4096);
//			}else if(bzipped && !RAWMODE){
//				
//				
//				in=new CBZip2InputStream(bis);
//				
//				/*
//				From http://www.kohsuke.org/bzip2/:
//				
//				Note
//				
//				Jacek Bilski told me that he had to read two bytes from the stream 
//				before he uses CBZip2InputStream. Those two bytes ('B' and 'Z')
//				 are used by the command line bzip program to mark the stream.
//				*/ 
//				
//				
//			}else{
//				in=bis;
//			}
//		} catch (FileNotFoundException e) {
//			throw new RuntimeException(e);
//		} catch (IOException e) {
//			throw new RuntimeException(e);
//		}
//		
//		return in;
//	}
	
//	public static InputStream getInputStream(String fname){
//		return getInputStream(fname, true);
//	}
//	
//	public static InputStream getInputStream(String fname, boolean buffer){
//		return getInputStream(fname, buffer, true);
//	}
	
	public static InputStream getInputStream(String fname, boolean buffer, boolean allowSubprocess){
		if(verbose){System.err.println("getInputStream("+fname+", "+buffer+", "+allowSubprocess+")");}
		boolean xz=fname.endsWith(".xz");
		boolean gzipped=fname.endsWith(".gz") || fname.endsWith(".gzip");
		boolean zipped=fname.endsWith(".zip");
		boolean bzipped=PROCESS_BZ2 && fname.endsWith(".bz2");
		boolean bam=fname.endsWith(".bam") && Data.SAMTOOLS();
		
		allowSubprocess=(allowSubprocess && Shared.THREADS>1);
		
		if(!RAWMODE){
			if(zipped){return getZipInputStream(fname);}
			if(gzipped){return getGZipInputStream(fname, allowSubprocess);}
			if(bzipped){return getBZipInputStream(fname);}
			if(bam){return getInputStreamFromProcess(fname, "samtools view -h", false);}
		}

		return getRawInputStream(fname, buffer);
	}
	
	public static InputStream getRawInputStream(String fname, boolean buffer){
		if(verbose){System.err.println("getRawInputStream("+fname+", "+buffer+")");}
		
		assert(fname!=null);
		fname=fname.replace('\\', '/');
		assert(fname.indexOf('\\')<0);
		assert(!fname.contains("\\\\"));
//		assert(!fname.contains("//")) : fname;
		
		final boolean jar=fname.startsWith("jar:");
		
		if(!jar){
			File f=new File(fname);
			if(!f.exists()){
				String f2=fname.toLowerCase();
				if(f2.equals("stdin") || f2.startsWith("stdin.")){
					//				System.err.println("Returning stdin: A");
					return System.in;
				}
				throw new RuntimeException("Can't find file "+fname);
			}
		}
		
//		System.err.println("Getting input stream for "+fname);
//		assert(!fname.contains("\\"));
//		assert(!loadedFiles.contains(fname)) : "Already loaded "+fname;
//		loadedFiles.add(fname);
//		assert(!fname.contains("custom_summary_unionGene_build36.txt"));
		
		
		InputStream in=null;
		if(jar){
			try {
				
				URL url=new URL(fname);
				
				InputStream is=url.openStream();

				if(buffer){
					BufferedInputStream bis=new BufferedInputStream(is, INBUF);
					in=bis;
				}else{
					in=is;
				}

			} catch (FileNotFoundException e) {
				System.err.println("Error when attempting to read "+fname);
				throw new RuntimeException(e);
			} catch (MalformedURLException e) {
				System.err.println("Error when attempting to read "+fname);
				throw new RuntimeException(e);
			} catch (IOException e) {
				System.err.println("Error when attempting to read "+fname);
				throw new RuntimeException(e);
			}
		}else{
			try {

				FileInputStream fis=new FileInputStream(fname);

				if(buffer){
					BufferedInputStream bis=new BufferedInputStream(fis, INBUF);
					in=bis;
				}else{
					in=fis;
				}

			} catch (FileNotFoundException e) {
				throw new RuntimeException(e);
			}
		}
		
		return in;
	}
	
	public static InputStream getZipInputStream(String fname){return getZipInputStream(fname, true);}
	public static InputStream getZipInputStream(String fname, boolean buffer){
		if(verbose){System.err.println("getZipInputStream("+fname+", "+buffer+")");}
		InputStream raw=getRawInputStream(fname, buffer);
		InputStream in=null;

		final String basename=basename(fname);

		try {

			ZipInputStream zis=new ZipInputStream(raw);
			ZipEntry ze=zis.getNextEntry();
			assert(ze!=null);
			assert(basename.equals(ze.getName())) : basename+" != "+ze.getName();
			in=zis;

		} catch (FileNotFoundException e) {
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		} catch (IOException e) {
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		}

		return in;
	}
	
	public static InputStream getGZipInputStream(String fname, boolean allowSubprocess){
		if(verbose){System.err.println("getGZipInputStream("+fname+", "+allowSubprocess+")");}
		if(allowSubprocess && Shared.THREADS>2){
			if(!fname.startsWith("jar:")){
				if(verbose){System.err.println("Fetching gzip input stream: "+fname+", "+allowSubprocess+", "+USE_UNPIGZ+", "+Data.UNPIGZ());}
				if(USE_UNPIGZ && Data.UNPIGZ()){return getUnpigzStream(fname);}
				if(USE_GUNZIP && Data.GUNZIP()){return getGunzipStream(fname);}
			}
		}
		
		InputStream raw=getRawInputStream(fname, false);
		InputStream in=null;
		
		try {
			in=new GZIPInputStream(raw, INBUF);
		} catch (FileNotFoundException e) {
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		} catch (IOException e) {
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		}

		return in;
	}
	
//	public static InputStream getGunzipStream(String fname){
//
//		//InputStream raw=getRawInputStream(fname, false);
//		InputStream in=null;
//
//		Process p=null;
//		if(fname.equals("stdin") || fname.startsWith("stdin.")){
//			try {
//				p = Runtime.getRuntime().exec("gzip -c -d");
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//			OutputStream os=p.getOutputStream();
//			PipeThread it=new PipeThread(System.in, os);
//			it.start();
//		}else{
//			try {
//				p = Runtime.getRuntime().exec("gzip -c -d "+fname);
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
//		
//		in=p.getInputStream();
//		InputStream es=p.getErrorStream();
//		assert(es!=null);
//		PipeThread et=new PipeThread(es, System.err);
//		et.start();
//
//		return in;
//	}
	
	public static InputStream getGunzipStream(String fname){
		if(verbose){System.err.println("getGunzipStream("+fname+")");}
		return getInputStreamFromProcess(fname, "gzip -c -d", false);
	}
	
	public static InputStream getUnpigzStream(String fname){
		if(verbose){System.err.println("getUnpigzStream("+fname+")");}
		return getInputStreamFromProcess(fname, "pigz -c -d", false);
	}
	
	public static InputStream getInputStreamFromProcess(String fname, String command, boolean cat){
		if(verbose){System.err.println("getInputStreamFromProcess("+fname+", "+command+", "+cat+")");}

		//InputStream raw=getRawInputStream(fname, false);
		InputStream in=null;

		Process p=null;
		if(fname.equals("stdin") || fname.startsWith("stdin.")){
			try {
				if(cat){
					throw new RuntimeException();
				}else{
					p=Runtime.getRuntime().exec(command);
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			OutputStream os=p.getOutputStream();
			PipeThread it=new PipeThread(System.in, os);
			addPipeThread(fname, it);
			it.start();
		}else{
			try {
				if(cat){
					assert(false) : "This mode is untested.";
					String[] cmd = {
							"sh","cat "+fname,
							" | "+command
					};
					p=Runtime.getRuntime().exec(cmd);
				}else{
					p = Runtime.getRuntime().exec(command+" "+fname);
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		addProcess(fname, p);
		in=p.getInputStream();
		InputStream es=p.getErrorStream();
		assert(es!=null);
		PipeThread et=new PipeThread(es, System.err);
		addPipeThread(fname, et);
		et.start();

		return in;
	}
	
	
	public static InputStream getBZipInputStream(String fname){
		if(verbose){System.err.println("getBZipInputStream("+fname+")");}
		InputStream in=null;
		
		try {in=getBZipInputStream(fname, true);} 
		catch (IOException e) {}
		
		if(in==null){
			try {in=getBZipInputStream(fname, false);} 
			catch (IOException e) {
				System.err.println("Error when attempting to read "+fname);
				throw new RuntimeException(e);
			}catch (NullPointerException e) {
				System.err.println("Error when attempting to read "+fname);
				throw new RuntimeException(e);
			}
		}
		
		assert(in!=null);
		return in;
	}
	
	
	private static InputStream getBZipInputStream(String fname, boolean stripBZ) throws IOException{
		if(verbose){System.err.println("getBZipInputStream("+fname+", "+stripBZ+")");}
		throw new RuntimeException("bz2 compression not supported in public version.");
//		InputStream raw=getRawInputStream(fname, true);
//		InputStream in=null;
//		
//		if(stripBZ){
////			System.err.println("Attempting to strip BZ");
//			byte[] header=new byte[2];
//			try {raw.read(header);} 
//			catch (IOException e) {throw new RuntimeException(e);}
//			
//			if(header[0]!='B' || header[1]!='Z'){
//				throw new IOException("Found BZ2 file that does not start with BZ: "+"("+header[0]+", "+header[1]+")");
//			}
//		}
//		
//		try {
//			
//			in=new CBZip2InputStream(raw);
//			
//			/*
//				From http://www.kohsuke.org/bzip2/:
//
//				Note
//
//				Jacek Bilski told me that he had to read two bytes from the stream 
//				before he uses CBZip2InputStream. Those two bytes ('B' and 'Z')
//				 are used by the command line bzip program to mark the stream.
//			 */ 
//			
//		} catch (Exception e) {
//			
//			try {
//				in=null;
//				raw.close();
//				raw=null;
//			} catch (IOException e1) {
//				System.err.println("Error when attempting to read "+fname);
//				e.printStackTrace();
//				e1.printStackTrace();
//				throw new RuntimeException(e1);
//			}
//		}
//
//		return in;
	}
	
	public static InputStream getXZInputStream(String fname){
		
		InputStream in=null;
		
//		if(PROCESS_XZ){
//			InputStream raw=getRawInputStream(fname, true);
//			try {
//				in=new org.tukaani.xz.XZInputStream(raw);
//			} catch (FileNotFoundException e) {
//				throw new RuntimeException(e);
//			} catch (IOException e) {
//				throw new RuntimeException(e);
//			}
//		}

		return in;
	}
	
	
	public static <X> X read(Class<X> cx, String fname){
		X x=(X)readObject(fname);
		return x;
	}
	
	public static <X> X[] readArray(Class<X> cx, String fname){
		X[] x=(X[])readObject(fname);
		return x;
	}
	
	public static <X> X[][] readArray2(Class<X> cx, String fname){
		X[][] x=(X[][])readObject(fname);
		return x;
	}
	
	public static <X> X[][][] readArray3(Class<X> cx, String fname){
		X[][][] x=(X[][][])readObject(fname);
		return x;
	}
	
	
	private static String basename(String fname){
		fname=fname.replace('\\', '/');
		boolean xz=fname.endsWith(".xz");
		boolean gzipped=fname.endsWith(".gz");
		boolean zipped=fname.endsWith(".zip");
		boolean bzipped=PROCESS_BZ2 && fname.endsWith(".bz2");
		String basename=fname;
//		if(basename.contains("\\")){basename=basename.substring(basename.lastIndexOf("\\")+1);}
		if(basename.contains("/")){basename=basename.substring(basename.lastIndexOf('/')+1);}
		if(zipped || bzipped){basename=basename.substring(0, basename.length()-4);}
		else if(gzipped){basename=basename.substring(0, basename.length()-3);}
		return basename;
	}
	
	public static String rawName(String fname){
		for(String s : extensions){
			while(fname.endsWith(s)){fname=fname.substring(0, fname.length()-s.length());}
		}
		return fname;
	}
	
	public static String compressionType(String fname){
		fname=fname.toLowerCase();
		for(String s : extensions){
			if(fname.endsWith(s)){return s.substring(1);}
		}
		return null;
	}
	
	public static boolean isCompressed(String fname){
		return compressionType(fname)!=null;
	}
	
	public static boolean isSam(String fname){
		if(fname.endsWith(".sam")){return true;}
		String s=compressionType(fname);
		if(s==null){return false;}
		return fname.substring(0, fname.lastIndexOf('.')).endsWith(".sam");
	}
	
	public static String rawExtension(String fname){
		fname=rawName(fname);
		int x=fname.lastIndexOf('.');
		if(x<0){return "";}
		return fname.substring(x+1);
	}
	
	public static String parseRoot(String path){
		File f=new File(path);
		if(f.isDirectory()){
			if(!path.endsWith(FILESEP)){
				path=path+FILESEP;
			}
			return path;
		}else if(f.isFile()){
			int slash=path.lastIndexOf(FILESEP);
			if(slash<0){
				return "";
			}else{
				return path.substring(0, slash+1);
			}
		}else{
			throw new RuntimeException("Can't find "+path); //Try using parseRoot2 instead.
		}
	}
	
	/** This one does not throw an exception for non-existing paths */
	public static String parseRoot2(String path){
		File f=new File(path);
		
		if(!f.exists()){
			if(path.endsWith(FILESEP)){return path;}
			int slash=path.lastIndexOf(FILESEP);
			if(slash<0){
				return "";
			}else{
				return path.substring(0, slash+1);
			}
		}
		
		if(f.isDirectory()){
			if(!path.endsWith(FILESEP)){
				path=path+FILESEP;
			}
			return path;
		}else if(f.isFile()){
			int slash=path.lastIndexOf(FILESEP);
			if(slash<0){
				return "";
			}else{
				return path.substring(0, slash+1);
			}
		}else{
			throw new RuntimeException("Can't find "+path);
		}
	}
	
	public static String findFileExtension(final String fname){

		File file=new File(fname);
		if(file.exists()){return fname;}

		String basename=fname, temp;
		if(fname.endsWith(".zip") || fname.endsWith(".gz") || (PROCESS_BZ2 && fname.endsWith(".bz2")) || (PROCESS_XZ && fname.endsWith(".xz"))){
			basename=fname.substring(0, fname.lastIndexOf('.'));
		}
		temp=basename;
		file=new File(temp);
		if(!file.exists()){
			temp=basename+".gz";
			file=new File(temp);
		}
//		System.err.println(temp+" "+(file.exists() ? " exists" : " does not exist"));
		if(!file.exists()){
			temp=basename+".zip";
			file=new File(temp);
		}
//		System.err.println(temp+" "+(file.exists() ? " exists" : " does not exist"));
		if(!file.exists() && PROCESS_BZ2){
			temp=basename+".bz2";
			file=new File(temp);
		}
//		System.err.println(temp+" "+(file.exists() ? " exists" : " does not exist"));
		if(!file.exists() && PROCESS_XZ){
			temp=basename+".xz";
			file=new File(temp);
		}
//		System.err.println(temp+" "+(file.exists() ? " exists" : " does not exist"));
		if(!file.exists()){temp=fname;}
		
		return temp;
	}
	
	public static synchronized void copyFile(String source, String dest){copyFile(source, dest, false);}
	public static synchronized void copyFile(String source, String dest, boolean createPathIfNeeded){
		
		assert(!new File(dest).exists()) : "Destination file already exists: "+dest;
		if(createPathIfNeeded){
			File parent=new File(dest).getParentFile();
			if(parent!=null && !parent.exists()){
				parent.mkdirs();
			}
		}
		
		final boolean oldRawmode=RAWMODE;
		if((source.endsWith(".zip") && dest.endsWith(".zip"))
				 || (source.endsWith(".gz") && dest.endsWith(".gz")
						 || (source.endsWith(".bz2") && dest.endsWith(".bz2"))
						 || (source.endsWith(".xz") && dest.endsWith(".xz")))){
			RAWMODE=true;
		}
		
		try{
			InputStream in=getInputStream(source, false, false);
			OutputStream out=getOutputStream(dest, false, false, true);

			byte[] buffer=new byte[INBUF];
			int len;
			
			while((len = in.read(buffer)) > 0){
				out.write(buffer, 0, len);
			}
			
			in.close();
			out.flush();
			if(out.getClass()==ZipOutputStream.class){
				ZipOutputStream zos=(ZipOutputStream)out;
				zos.closeEntry();
				zos.finish();
			}
//			else if(PROCESS_XZ && out.getClass()==org.tukaani.xz.XZOutputStream.class){
//				org.tukaani.xz.XZOutputStream zos=(org.tukaani.xz.XZOutputStream)out;
//				zos.finish();
//			}
			out.close();
			
		}catch(FileNotFoundException e){
			RAWMODE=oldRawmode;
			throw new RuntimeException(e);
		}catch(IOException e){
			RAWMODE=oldRawmode;
			throw new RuntimeException(e);   
		}
		
		RAWMODE=oldRawmode;
	}
	
	public static void copyDirectoryContents(String from, String to){
		assert(!from.equalsIgnoreCase(to));
		
		if(to.indexOf('\\')>0){to=to.replace('\\', '/');}
		
		File d1=new File(from);
		assert(d1.exists());
		assert(d1.isDirectory());
		
		File d2=new File(to);
		assert(!d1.equals(d2));
		if(d2.exists()){
			assert(d2.isDirectory());
		}else{
			d2.mkdirs();
		}
		if(!to.endsWith("/")){to=to+"/";}
		
		File[] array=d1.listFiles();
		
		for(File f : array){
			String name=f.getName();
			String dest=to+name;
			if(f.isFile()){
				copyFile(f.getAbsolutePath(), dest);
			}else{
				assert(f.isDirectory());
				File f2=new File(dest);
				if(!f2.exists()){
					f2.mkdir();
				}else{
					assert(f2.isDirectory());
				}
				copyDirectoryContents(f.getAbsolutePath(), f2.getAbsolutePath());
			}
		}
		
	}
	
	
	private static final int addThread(int x){
		if(verbose){System.err.println("addThread("+x+")");}
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				activeThreads[0]+=x;
				activeThreads[1]+=x;
			}else{
				addRunningThread(x);
			}
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 && 
					activeThreads[2]>=0 && activeThreads[2]<=maxWriteThreads) : Arrays.toString(activeThreads);
					
			return activeThreads[0];
		}
	}
	
	private static final int addRunningThread(int x){
		if(verbose){System.err.println("addRunningThread("+x+")");}
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				assert(activeThreads[1]>=x);
				while(activeThreads[2]>=maxWriteThreads){
					try {
						activeThreads.wait();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				activeThreads[1]-=x; //Remove from waiting
			}else{
				activeThreads[0]+=x; //Remove from active
			}
			activeThreads[2]+=x; //Change number running
			
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 && 
					activeThreads[2]>=0 && activeThreads[2]<=maxWriteThreads) : Arrays.toString(activeThreads);
			
			if(activeThreads[2]==0 || (activeThreads[2]<maxWriteThreads && activeThreads[1]>0)){activeThreads.notify();}
			return activeThreads[2];
		}
	}
	
	public static final int countActiveThreads(){
		if(verbose){System.err.println("countActiveThreads()");}
		synchronized(activeThreads){
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 && 
					activeThreads[2]>=0 && activeThreads[2]<=maxWriteThreads) : Arrays.toString(activeThreads);
			return activeThreads[0];
		}
	}
	
	public static final void waitForWritingToFinish(){
		if(verbose){System.err.println("waitForWritingToFinish()");}
		synchronized(activeThreads){
			while(activeThreads[0]>0){
				assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 && 
						activeThreads[2]>=0 && activeThreads[2]<=maxWriteThreads) : Arrays.toString(activeThreads);
				try {
					activeThreads.wait(8000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(activeThreads[2]==0 || (activeThreads[2]<maxWriteThreads && activeThreads[1]>0)){activeThreads.notify();}
			}
		}
	}

	
	public static final boolean closeStream(ConcurrentReadStreamInterface cris){return closeStreams(cris, (RTextOutputStream3[])null);}
	public static final boolean closeStream(RTextOutputStream3 ross){return closeStreams(null, ross);}
	
	public static final boolean closeStreams(ConcurrentReadStreamInterface cris, RTextOutputStream3...ross){
		if(verbose){
			System.err.println("closeStreams("+cris+", "+(ross==null ? "null" : ross.length)+")");
			new Exception().printStackTrace(System.err);
		}
		boolean errorState=false;
		if(cris!=null){
			if(verbose){System.err.println("Closing cris; error="+errorState);}
			cris.close();
			errorState|=cris.errorState();
//			Object[] prods=cris.producers();
//			for(Object o : prods){
//				if(o!=null && o.getClass()==ReadInputStream.class){
//					ReadInputStream ris=(ReadInputStream)o;
//					ris.
//				}
//			}
			if(verbose){System.err.println("Closed cris; error="+errorState);}
		}
		if(ross!=null){
			for(RTextOutputStream3 ros : ross){
				if(verbose){System.err.println("Closing ros "+ros+"; error="+errorState);}
				if(ros!=null){
					ros.close();
					ros.join();
					errorState|=(ros.errorState() || !ros.finishedSuccessfully());
				}
				if(verbose){System.err.println("Closed ros; error="+errorState);}
			}
		}
		return errorState;
	}
	
	public static boolean killProcess(String fname){
		if(verbose){
			System.err.println("killProcess("+fname+")");
			new Exception().printStackTrace(System.err);
		}
		if(fname==null || (!isCompressed(fname) && !fname.endsWith(".bam"))){return false;}
		
		boolean error=false;
		synchronized(processMap){
			Process p=processMap.remove(fname);
			if(p!=null){
				if(verbose){System.err.println("Found Process for "+fname);}
				int x=-1, tries=0; 
				for(; tries<20; tries++){
					if(verbose){System.err.println("Trying p.waitFor()");}
					try {
						x=p.waitFor();
						if(verbose){System.err.println("success; return="+x);}
						break;
					} catch (InterruptedException e) {
						if(verbose){System.err.println("Failed.");}
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				error|=(tries>=20 || x!=0);
				if(tries>=20){
					if(verbose){System.err.println("Calling p.destroy because tries=="+tries+"; error="+error);}
					p.destroy();
					if(verbose){System.err.println("destroyed");}
				}
			}
		}
		synchronized(pipeThreadMap){
			ArrayList<PipeThread> atp=pipeThreadMap.remove(fname);
			if(atp!=null){
				for(PipeThread p : atp){
					if(p!=null){
						if(verbose){System.err.println("Found PipeThread for "+fname);}
						p.terminate();
						if(verbose){System.err.println("Terminated PipeThread");}
					}
				}
			}
		}
		if(verbose){System.err.println("killProcess("+fname+") returned "+error);}
		return error;
	}
	
	private static void addProcess(String fname, Process p){
		if(verbose){System.err.println("addProcess("+fname+", "+p+")");}
		synchronized(processMap){
//			System.err.println("Adding Process for "+fname);
			Process old=processMap.put(fname, p);
			if(old!=null){
				old.destroy();
				throw new RuntimeException("Duplicate process for file "+fname);
			}
		}
	}
	
	private static void addPipeThread(String fname, PipeThread pt){
		if(verbose){System.err.println("addPipeThread("+fname+", "+pt+")");}
		synchronized(pipeThreadMap){
//			System.err.println("Adding PipeThread for "+fname);
			ArrayList<PipeThread> atp=pipeThreadMap.get(fname);
			if(atp==null){
				atp=new ArrayList<PipeThread>(2);
				pipeThreadMap.put(fname, atp);
			}
			atp.add(pt);
		}
	}
	
	/** {active, waiting, running} <br>
	 * Active means running or waiting.
	 */
	public static int[] activeThreads={0, 0, 0};
	public static int maxWriteThreads=Shared.THREADS;
	
	public static boolean verbose=false;
	
	public static boolean RAWMODE=false; //Does not automatically compress and decompress when true

	public static boolean USE_GZIP=false;
	public static boolean USE_PIGZ=false;
	public static boolean USE_GUNZIP=false;
	public static boolean USE_UNPIGZ=false;

	public static boolean PROCESS_BZ2=true;
	public static final boolean PROCESS_XZ=false;
	
	public static final int INBUF=16384;
	public static final int OUTBUF=16384;

	public static int ZIPLEVEL=4;
	public static int MAX_ZIP_THREADS=8;
	public static int ZIP_THREAD_DIVISOR=4;
	
	public static final String FILESEP=System.getProperty("file.separator");

	private static final String diskSync=new String("DISKSYNC");
	
	public static final HashSet<String> loadedFiles=new HashSet<String>();
	
	private static final String[] extensions=new String[] {".gz", ".gzip", ".zip", ".bz2", ".xz"};

//	private static HashMap<String, Process> inputProcesses=new HashMap<String, Process>(8);
//	private static HashMap<String, Process> outputProcesses=new HashMap<String, Process>(8);
	private static HashMap<String, Process> processMap=new HashMap<String, Process>(8);
	private static HashMap<String, ArrayList<PipeThread>> pipeThreadMap=new HashMap<String, ArrayList<PipeThread>>(8);
	
}
