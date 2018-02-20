package fileIO;

import java.io.File;

import dna.Data;


public class RenameFiles {
	
	
	public static void main(String[] args){
		for(String s : args){
			renameFiles(s);
		}
	}
	
	
	public static void renameFiles(String path){
		File f=new File(path);
		renameFiles(f);
	}
	
	public static void renameFiles(File path){
		
		if(path.isDirectory()){
			File[] array=path.listFiles();
			for(File f : array){renameFiles(f);}
		}else{
			rename(path);
		}
		
	}
	
	public static void rename(File in){
		assert(in.exists()) : in.toString();
		assert(in.isFile()) : in.toString();
		String abs=in.getAbsolutePath();
		
		
		int dot=abs.lastIndexOf('.');
		int slash=abs.lastIndexOf('/');
		
//		String[] split=Person.parsePath(abs.substring(0, slash));
//		String name=split[0];
//		String out=abs.substring(0, dot)+"_"+name+".txt";
		
		
		
		String fname=abs.substring(slash+1);
		
//		System.out.println(fname);
		

//		if(fname.startsWith("chr") && fname.endsWith(".txt")){
//			
//			String out=abs.replace(".txt", ".flow");
//			assert(!out.equals(abs)) : out+", "+abs;
//			
//			System.out.println("Renaming "+abs+" to "+out);
//			in.renameTo(new File(out));
//		}
		
		int build=36;
		if(abs.contains("FL5-") || abs.contains("630-") || abs.contains("618-")){
			build=37;
		}
		
//		if(fname.startsWith("var") && fname.endsWith(".vla") && !fname.contains("build")){
//			
//			String out=abs.replace(".vla", "-build"+build+".vla");
//			assert(!out.equals(abs)) : out+", "+abs;
//			
//			System.out.println("Renaming "+abs+" to "+out);
//			in.renameTo(new File(out));
//		}
//		
//		if(fname.startsWith("gene") && fname.endsWith(".gvla") && !fname.contains("build")){
//			
//			String out=abs.replace(".gvla", "-build"+build+".gvla");
//			assert(!out.equals(abs)) : out+", "+abs;
//			
//			System.out.println("Renaming "+abs+" to "+out);
//			in.renameTo(new File(out));
//		}
//		
//		if(fname.endsWith(".tsv.zip") && !fname.contains("build")){
//			
//			String out=abs.replace(".tsv.zip", "-build"+build+".tsv.zip");
//			assert(!out.equals(abs)) : out+", "+abs;
//			
//			System.out.println("Renaming "+abs+" to "+out);
//			in.renameTo(new File(out));
//		}
//		
//		if(fname.endsWith(".tsv.gz") && !fname.contains("build")){
//			
//			String out=abs.replace(".tsv.gz", "-build"+build+".tsv.gz");
//			assert(!out.equals(abs)) : out+", "+abs;
//			
//			System.out.println("Renaming "+abs+" to "+out);
//			in.renameTo(new File(out));
//		}
//		
//		if(fname.endsWith(".tsv") && !fname.contains("build")){
//			
//			String out=abs.replace(".tsv", "-build"+build+".tsv");
//			assert(!out.equals(abs)) : out+", "+abs;
//			
//			System.out.println("Renaming "+abs+" to "+out);
//			in.renameTo(new File(out));
//		}
//		
//		if(fname.endsWith(".ca") && !fname.contains("build")){
//			
//			String out=abs.replace(".ca", "-build"+build+".ca");
//			assert(!out.equals(abs)) : out+", "+abs;
//			
//			System.out.println("Renaming "+abs+" to "+out);
//			in.renameTo(new File(out));
//		}
//		
//		if(fname.endsWith(".ca.zip") && !fname.contains("build")){
//			
//			String out=abs.replace(".ca.zip", "-build"+build+".ca.zip");
//			assert(!out.equals(abs)) : out+", "+abs;
//			
//			System.out.println("Renaming "+abs+" to "+out);
//			in.renameTo(new File(out));
//		}
//		
//		if(fname.contains("-ASM-") && fname.contains("build36")){
//			
//			String out=abs.replace("-build36", "");
//			assert(!out.equals(abs)) : out+", "+abs;
//			
//			System.out.println("Renaming "+abs+" to "+out);
//			in.renameTo(new File(out));
//		}
		
		if(fname.contains("READMEtxt")){
			String out=abs.replace("READMEtxt", "README.txt");
			assert(!out.equals(abs)) : out+", "+abs;
			
			System.out.println("Renaming "+abs+" to "+out);
			in.renameTo(new File(out));
		}
		
		if(fname.contains("-1.8.0.")){
			
			String out=abs.replace("-1.8.0.", ".");
			assert(!out.equals(abs)) : out+", "+abs;
			
			System.out.println("Renaming "+abs+" to "+out);
			in.renameTo(new File(out));
		}
	}
	
}
