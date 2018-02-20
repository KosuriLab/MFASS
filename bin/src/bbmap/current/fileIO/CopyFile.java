package fileIO;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.zip.ZipOutputStream;

import align2.Tools;

import dna.Timer;


/**
 * Unlike ReadWrite's version, this one forces compression and decompression even with same extensions.
 * Mainly for benchmarking.
 * @author Brian Bushnell
 * @date Jan 23, 2013
 *
 */
public class CopyFile {
	
	public static void main(String[] args){

		String in=null, out=null;
		boolean overwrite=true;

		for(int i=0; i<args.length; i++){

			if(true){
				final String arg=args[i];
				final String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if("null".equalsIgnoreCase(b)){b=null;}

				if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
					//jvm argument; do nothing
				}else if(a.equals("in")){
					in=b;
				}else if(a.equals("out")){
					out=b;
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
					ReadWrite.USE_UNPIGZ=Tools.parseBoolean(b);;
				}else if(a.equals("overwrite") || a.equals("ow")){
					overwrite=Tools.parseBoolean(b);
				}else if(a.equals("ziplevel") || a.equals("zl")){
					ReadWrite.ZIPLEVEL=Integer.parseInt(b);
				}else if(in==null && i==0 && !args[i].contains("=")){
					in=args[i];
				}else if(out==null && i==1 && !args[i].contains("=")){
					out=args[i];
				}
			}
		}
		assert(in!=null && out!=null);
		long bytes=new File(in).length();
		Timer t=new Timer();
		t.start();
		copyFile(in, out, false, overwrite);
		t.stop();
		double mbps1=bytes*1000d/t.elapsed;
		System.err.println("Time:  \t"+t);
		System.err.println(String.format("Speed: \t%.2f MB/s", mbps1));
	}
	
	
	public static synchronized void copyFile(String source, String dest, boolean createPathIfNeeded, boolean overwrite){

		assert(overwrite || !new File(dest).exists()) : "Destination file already exists: "+dest;
		if(createPathIfNeeded){
			File parent=new File(dest).getParentFile();
			if(parent!=null && !parent.exists()){
				parent.mkdirs();
			}
		}

		try{
			InputStream in=ReadWrite.getInputStream(source, false, true);
			OutputStream out=ReadWrite.getOutputStream(dest, false, false, true);

			final byte[] buffer=new byte[16384];
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
			throw new RuntimeException(e);
		}catch(IOException e){
			throw new RuntimeException(e);   
		}
	}
	
}
