package dna;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URL;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.HashMap;

import kmer.Primes;

import var.Variation;

import align2.AbstractIndex;
import align2.ChromLoadThread;
import align2.Tools;


import driver.Search;
import fileIO.ChainBlock;
import fileIO.ChainLine;
import fileIO.ReadWrite;
import fileIO.TextFile;

public class Data {
	
	
	public static void main(String[] args){}
	
	
	//TODO IMPORTANT!  Ensure that this unloads everything big, AND that reloading subsequently works OK.
	public static void unloadAll(){
		chromosomePlusMatrix=null;
		chromosomeMinusMatrix=null;
		AbstractIndex.clear();
	}
	
	
	//TODO IMPORTANT!  Ensure that this unloads everything big, AND that reloading subsequently works OK.
	public static void unload(int chrom, boolean unloadSoft){

//		unloadGenes(chrom);
		
		chromosomePlusMatrix[chrom]=null;
	}
	
	public static void unloadGenes(int chrom){
		geneMatrix[chrom]=null;
		geneSetMatrix[chrom]=null;
		geneTxRangeMatrix[chrom]=null;
		geneSetRangeMatrix[chrom]=null;
		geneCodeRangeMatrix[chrom]=null;
		geneCodeAndExonRangeMatrix[chrom]=null;
		geneNearbyRangeMatrix[chrom]=null;
		exonRangeMatrix[chrom]=null;
	}
	
	public static byte find(int x, byte[] array){
		for(byte i=0; i<array.length; i++){
			if(array[i]==x){return i;}
		}
		return -1;
	}
	
	public static void reverse(byte[] array){
		int mid=array.length/2;
		for(int i=0; i<mid; i++){
			byte temp=array[i];
			array[i]=array[array.length-i-1];
			array[array.length-i-1]=temp;
		}
	}
	
	
	public static Gene[] getGenes(int chrom){
		if(geneMatrix[chrom]==null){
			loadGenes(chrom);
		}
		return geneMatrix[chrom];
	}
	
	
	public static Gene[] getGenes(int chrom, byte strand){
		ArrayList<Gene> genes=new ArrayList<Gene>();
		for(Gene g : getGenes(chrom)){
			if(g.strand==strand){
				genes.add(g);
			}
		}
		return genes.toArray(new Gene[genes.size()]);
	}
	
	
	public static GeneSet[] getGeneSets(int chrom){
		if(geneSetMatrix[chrom]==null){
			loadGenes(chrom);
		}
		return geneSetMatrix[chrom];
	}
	
	
	public static HashMap<Integer, ArrayList<GeneSet>> getGeneIDTable(){
		if(geneIDTable==null){

//			System.err.println("WAITING FOR CS");
			synchronized(GENEIDLOCK){
//				System.err.println("ENTER CS");
				if(geneIDTable==null){
//					System.err.println("ENTER CS2");
					HashMap<Integer, ArrayList<GeneSet>> temp=new HashMap<Integer, ArrayList<GeneSet>>(2048);
					for(byte chrom=1; chrom<=25; chrom++){
						GeneSet[] set=getGeneSets(chrom);
						for(GeneSet gs : set){
							int id=-1;
							for(Gene g : gs.genes){
								if(id==-1){id=g.id;}
								else{assert(id==g.id) : gs+"\n"+gs.genes+"\n";}
							}
							assert(id>-1);
							Integer key=new Integer(id);
							ArrayList<GeneSet> value=temp.get(key);
							//					assert(old==null || chrom>22) : "\nCollision!\n\n"+gs+"\n\nis overwriting\n\n"+old;
							if(value==null){
								value=new ArrayList<GeneSet>(2);
								temp.put(key, value);
							}
							value.add(gs);
						}
					}
//					System.err.println("EXIT CS2");
					geneIDTable=temp;
				}
//				System.err.println("EXIT CS");
			}

		}
//		System.err.println("GeneIDTable contains "+geneIDTable.size()+" entries.");
		return geneIDTable;
	}
	
//	public static ChromosomeArray getChromosome(int chrom, boolean colorspace){
//		assert(!colorspace) : "Colorspace no longer supported";
//		if(colorspace){throw new RuntimeException("Colorspace no longer supported");}
//		if(true){throw new RuntimeException("Colorspace no longer supported");}
//		return getChromosome(chrom);
//	}
	
	public static ChromosomeArray getChromosome(int chrom, byte strand){
		if(strand==Gene.PLUS){return getChromosome(chrom);}
		
		if(chromosomeMinusMatrix[chrom]==null){
			synchronized(CHROMLOCKS[chrom%CHROMLOCKS.length]){
				if(chromosomeMinusMatrix[chrom]==null){
					ChromosomeArray p=getChromosome(chrom);
					chromosomeMinusMatrix[chrom]=p.complement();
				}
			}
		}
		return chromosomeMinusMatrix[chrom];
	}
	
	public static ChromosomeArray getChromosome(int chrom){
		assert(chromosomePlusMatrix!=null); 
		assert(chromosomePlusMatrix.length>chrom) : chromosomePlusMatrix.length+", "+chrom; 
		if(chromosomePlusMatrix[chrom]==null){
			synchronized(CHROMLOCKS[chrom%CHROMLOCKS.length]){
				if(chromosomePlusMatrix[chrom]==null){loadChromosome(chrom);}
			}
		}
		assert(chromosomePlusMatrix[chrom].array[0]=='N') : (char)chromosomePlusMatrix[chrom].array[0]+
			"\nIf you see this message, please regenerate your index.";//startpad was too low or for some reason invalid.
		return chromosomePlusMatrix[chrom];
	}
	
	private static void loadGenes(int chrom){
		
		if(geneMatrix[chrom]!=null){return;} //In case another thread already loaded the chromosome
		synchronized(CHROMLOCKS[chrom%CHROMLOCKS.length]){
			if(geneMatrix[chrom]==null){

				//		Gene[] genes=FindExons.readGenes(ROOT_GENE+"ref/chr"+Gene.chromCodes[chrom]+".Ref.Table", Gene.FORMAT_NM);
				//		Gene[] genes=FindExons.readGenes(ROOT_GENE+"ref2/ccds-chr"+Gene.chromCodes[chrom]+"-genes.txt", Gene.FORMAT_CCDS);
				//		Gene[] genes=FindExons.readGenes(ROOT_GENE+"ref3/ccds-chr"+Gene.chromCodes[chrom]+"-genes.txt", Gene.FORMAT_CCDS);

				//		Gene[] genes=ReadWrite.readArray(Gene.class, ROOT_GENE+"seqGene/chr"+Gene.chromCodes[chrom]+".ga");
				
				String fname=ROOT_GENE+"Build"+GENOME_BUILD+"/"+GENE_MAP+"/chr"+Gene.chromCodes[chrom]+".ga";
				
				Gene[] genes=ReadWrite.readArray(Gene.class, fname);

				Arrays.sort(genes);
//				geneMatrix[chrom]=genes;

				geneTxRangeMatrix[chrom]=findGeneRanges(genes, TX_RANGE);
				geneCodeRangeMatrix[chrom]=findGeneRanges(genes, CODE_RANGE);
				geneCodeAndExonRangeMatrix[chrom]=findCodeAndExonRanges(genes, false, true);
				exonRangeMatrix[chrom]=findCodeAndExonRanges(genes, false, false);
				geneNearbyRangeMatrix[chrom]=findCodeAndExonRanges(genes, true, true);

				HashMap<String, ArrayList<Gene>> temp=new HashMap<String, ArrayList<Gene>>();
				HashMap<String, GeneSet> gntable=new HashMap<String, GeneSet>();
				HashMap<String, Gene> tntable=new HashMap<String, Gene>();

				for(Gene g : genes){

					String trkey=g.mrnaAcc;
					if(trkey==null){trkey=g.chromosome+"_"+g.id;}
					if(trkey!=null){
						Gene old=tntable.get(trkey);
						if(old!=null){
							//					stdout.println("For transcript '"+g.nameTranscript+"': Overwrote \n"+old+"\nwith\n"+g+"\n");
						}
						tntable.put(trkey, g);
					}

					String key=g.symbol;
					if(key==null){key=g.mrnaAcc;}
					ArrayList<Gene> list=temp.get(key);
					if(list==null){
						list=new ArrayList<Gene>(8);
						temp.put(key, list);
					}
					list.add(g);
				}

				GeneSet[] gsm=new GeneSet[temp.size()];
				String[] keys=temp.keySet().toArray(new String[temp.size()]);
				for(int i=0; i<keys.length; i++){
					String key=keys[i];
					ArrayList<Gene> list=temp.get(key);
					GeneSet gs=new GeneSet(key, list);
					gsm[i]=gs;
					gntable.put(key, gs);
				}

				geneNameTable[chrom]=gntable;
				transcriptNameTable[chrom]=tntable;
				geneSetMatrix[chrom]=gsm;
				Arrays.sort(geneSetMatrix[chrom]);

				geneSetRangeMatrix[chrom]=findGeneSetRanges(geneSetMatrix[chrom]);
				
				{
					assert(geneMatrix[chrom]==null) : "Need to sync.";
					geneMatrix[chrom]=genes;
				}
			}
		}
	}
	
	public static void loadChromosomes(int min, int max){
//		assert(false);
		synchronized(CHROMLOCKS){
			String pattern=chromFname(GENOME_BUILD);
			ChromLoadThread.loadAll(pattern, min, max, chromosomePlusMatrix);
		}
	}
	
	private static void loadChromosome(int chrom){
//		assert(false);
		assert(chromosomePlusMatrix[chrom]==null);
		assert(chrom>0) : chrom;
		
		String fname=chromFname(chrom, GENOME_BUILD);
		sysout.println("Loading "+fname);
		if(CHROMC){
			chromosomePlusMatrix[chrom]=ReadWrite.read(ChromosomeArrayCompressed.class, fname).toChromosomeArray();
		}else{
			chromosomePlusMatrix[chrom]=ReadWrite.read(ChromosomeArray.class, fname);
		}
		assert(chromosomePlusMatrix[chrom].chromosome==chrom);
	}
	
	public static final String chromExtension(){
		return ".chrom"+(CHROMC ? "C" : "") + (CHROMGZ ? ".gz" : "");
	}
	
	public static final String chromFname(int chrom, int genome){
		return ROOT_GENOME+genome+"/chr"+chrom+chromExtension();
	}
	
	public static final String chromFname(int genome){
		return ROOT_GENOME+genome+"/chr#"+chromExtension();
	}
	
	public static Range[] findGeneRanges(Gene[] genes, final int mode){
		
		ArrayList<Range> list=new ArrayList<Range>(8192);
		ArrayList<Gene> glist=new ArrayList<Gene>(64);
		
		Range current=null;
		for(int i=0; i<genes.length; i++){
			Gene g=genes[i];
			Range r;
			
			int a, b;
			
			switch(mode){
				case TX_RANGE: {a=g.txStart; b=g.txStop;}
				break;
				
				case CODE_RANGE: {a=g.codeStart; b=g.codeStop;}
				break;
				
				default: {throw new RuntimeException();}
			}
			
			if(b>=a){
				r=new Range(a, b);

				if(current==null){
					current=r;
					glist.add(g);
				}else if(current.touches(r)){
					current=current.merge(r);
					glist.add(g);
				}else{
					current.obj1=glist.toArray(new Gene[glist.size()]);
					glist.clear();
					glist.add(g);
					list.add(current);
					current=r;
				}
			}
		}
		if(current!=null){ //i.e., if there were any genes
			current.obj1=glist.toArray(new Gene[glist.size()]);
			list.add(current);
		}
		
		return list.toArray(new Range[list.size()]);
	}
	
	public static Range[] findGeneSetRanges(GeneSet[] genes){
		
		ArrayList<Range> list=new ArrayList<Range>(8192);
		ArrayList<GeneSet> glist=new ArrayList<GeneSet>(64);
		
		Range current=null;
		for(int i=0; i<genes.length; i++){
			GeneSet g=genes[i];
			Range r;
			
			int a=g.minStart-NEAR, b=g.maxEnd+NEAR;
			
			if(b>=a){
				r=new Range(a, b);

				if(current==null){
					current=r;
					glist.add(g);
				}else if(current.touches(r)){
					current=current.merge(r);
					glist.add(g);
				}else{
					current.obj1=glist.toArray(new GeneSet[glist.size()]);
					glist.clear();
					glist.add(g);
					list.add(current);
					current=r;
				}
			}
		}
		if(current!=null){ //i.e., if there were any genes
			current.obj1=glist.toArray(new GeneSet[glist.size()]);
			list.add(current);
		}
		
		return list.toArray(new Range[list.size()]);
	}
	
	
	public static Range[] findCodeAndExonRanges(Gene[] genes, boolean nearby, boolean codingOnly){
		
		
		ArrayList<Range> list=new ArrayList<Range>(32768);
		
		for(int i=0; i<genes.length; i++){
			Gene g=genes[i];
			Range r;
			
			for(Exon ex : g.exons){

				int a=ex.a, b=ex.b;
				
				if(codingOnly){
					a=max(ex.a, g.codeStart);
					b=min(ex.b, g.codeStop);
				}
				
				assert(ex.a<=ex.b);
				assert(g.codeStart<=g.codeStop+1) : g;
				
				if(nearby){
					a=a-NEAR;
					b=b+NEAR;
				}
				
				if(a<=b){
					r=new Range(a, b);
					r.obj1=g;
//					r.obj2=ex;
					list.add(r);
					
				}
				
			}
		}
		
		ArrayList<Range> list2=new ArrayList<Range>(list.size());
		Collections.sort(list);
		
		
		HashSet<Gene> gset=new HashSet<Gene>(64);
		Range current=null;
		for(Range r : list){
			if(current==null){
				gset.add((Gene)r.obj1);
				current=r;
			}else if(current.touches(r)){
				gset.add((Gene)r.obj1);
				current=current.merge(r);
			}else{
				current.obj1=gset.toArray(new Gene[gset.size()]);
				list2.add(current);
				gset.clear();
				gset.add((Gene)r.obj1);
				current=r;
			}
		}
		
		if(current!=null){
			current.obj1=gset.toArray(new Gene[gset.size()]);
			list2.add(current);
			Collections.sort(list2);
		}
		
		return list2.toArray(new Range[list2.size()]);
	}
	
	public static Range[] geneSetRangeMatrix(int chrom){
		if(geneSetRangeMatrix[chrom]==null){
			loadGenes(chrom);
		}
		assert(geneSetRangeMatrix[chrom]!=null);
		return geneSetRangeMatrix[chrom];
	}
	
	public static Range[] exonRangeMatrix(int chrom){
		if(exonRangeMatrix[chrom]==null){
			loadGenes(chrom);
		}
		assert(exonRangeMatrix[chrom]!=null);
		return exonRangeMatrix[chrom];
	}
	
	public static Range[] geneCodeAndExonRangeMatrix(int chrom){
		if(geneCodeAndExonRangeMatrix[chrom]==null){
			loadGenes(chrom);
		}
		assert(geneCodeAndExonRangeMatrix[chrom]!=null);
		return geneCodeAndExonRangeMatrix[chrom];
	}
	
	public static Range[] geneNearbyRangeMatrix(int chrom){
		if(geneNearbyRangeMatrix[chrom]==null){
			loadGenes(chrom);
		}
		assert(geneNearbyRangeMatrix[chrom]!=null);
		return geneNearbyRangeMatrix[chrom];
	}
	
	public static HashMap<String, GeneSet> geneNameTable(int chrom){
		if(geneNameTable[chrom]==null){
			loadGenes(chrom);
		}
		assert(geneNameTable[chrom]!=null);
		return geneNameTable[chrom];
	}
	
	public static HashMap<String, Gene> transcriptNameTable(int chrom){
		if(transcriptNameTable[chrom]==null){
			loadGenes(chrom);
		}
		assert(transcriptNameTable[chrom]!=null);
		return transcriptNameTable[chrom];
	}
	
	
	public static GeneSet[] getNearestGeneSets(int chrom, int loc){
		Range[] r=geneSetRangeMatrix(chrom);
		int index=driver.Search.findPointBinary(loc, r);
		GeneSet[] sets=(GeneSet[]) r[index].obj1;
		if(sets==null || sets.length==0){
			assert(false);
			return null;
		}
		return sets;
	}
	
	/** Returns genesets overlapping the range */
	public static GeneSet[] getNearestGeneSets(int chrom, int loc1, int loc2){
		assert(loc2>=loc1);
		
//		boolean flag=(chrom==21 && loc1<38540895 && loc2>38540895);//TODO UNDO
//		
//		if(flag){
//			stdout.println(loc1+", "+loc2+", "+((loc1+loc2)/2));
//			for(GeneSet gs : Data.geneNameTable[chrom].values()){
//				if(gs.intersects(loc1, loc2)){
//					stdout.println("%%% "+gs);
//				}
//			}
//		}
		
		Range[] ranges=geneSetRangeMatrix(chrom);
		if(ranges==null || ranges.length==0){return null;}
		int index=driver.Search.findPointBinary(loc1, ranges);
		
		
//		if(flag){
//			Range r0=ranges[index-1];
//			Range r1=ranges[index];
//			Range r2=ranges[index+1];
//
//			stdout.println("r0: "+r0+"\n"+Arrays.toString((GeneSet[])r0.obj1)+"\n");
//			stdout.println("r1: "+r1+"\n"+Arrays.toString((GeneSet[])r1.obj1)+"\n");
//			stdout.println("r2: "+r2+"\n"+Arrays.toString((GeneSet[])r2.obj1)+"\n");
//			
//		}
//		if(flag){stdout.println("c");}
		
		Range r1=ranges[index];
		Range r2=(index>=ranges.length-1 ? null : ranges[index+1]);
		
		if(ranges[index].b>=loc2 || r2==null || r2.a>loc2){
			return (GeneSet[])r1.obj1;
		}
		
////		if(flag){stdout.println("e");}
//		if(ranges[index].b>=loc2 || (index==ranges.length-1) || ranges[index+1].a>loc2){
////			if(flag){
////				stdout.println("f");
////				stdout.println(ranges[index].b<=loc2);
////				stdout.println((index==ranges.length-1));
////				stdout.println(ranges[index+1].a>loc2);
////				stdout.println(".......");
////			}
//			return sets1;
//		}
		
		if(loc1>r1.b && loc2<r2.a){
			//No-man's land: Return closer of the bounding ranges.
			int dist1=loc1-r1.b;
			int dist2=r2.a-loc2;
			if(dist1>=dist2){
				return (GeneSet[])r1.obj1;
			}else{
				return (GeneSet[])r2.obj1;
			}
		}
		
//		assert(false) : "Test: This should be very rare, since it is slow.";
		
		//Otherwise, return all overlapping ranges.
		ArrayList<GeneSet> list=new ArrayList<GeneSet>(4);
		
		while(index<ranges.length && ranges[index].b<loc1){index++;} //Spin until in range
		for(; index<ranges.length && loc2>=ranges[index].a; index++){
//			if(flag){stdout.println("ADDED RANGE "+ranges[index]);}
			GeneSet[] gsa=(GeneSet[]) ranges[index].obj1;
			for(GeneSet gs : gsa){list.add(gs);}
		}
		return list.toArray(new GeneSet[list.size()]);
	}
	
	
	public static boolean isExonic(byte chrom, int point, int thresh, boolean isCoding){
		Range[] ranges=(isCoding ? Data.geneCodeAndExonRangeMatrix(chrom) : Data.exonRangeMatrix(chrom));
		return Search.containsPointBinary(point, ranges, thresh);
	}
	

	public static final String padFront(String num, int width, String symbol){
		String r=num;
		while(r.length()<width){r=symbol+r;}
		return r;
	}
	
	public static final String toBinaryString(long num, int width){
		String r=Long.toBinaryString(num);
		while(r.length()<width){r="0"+r;}
		return r;
	}
	
	public static final String toString(double[][] a){
		StringBuilder sb=new StringBuilder(256);
		sb.append("[\n");
		for(double[] b : a){
			sb.append(" ").append(Arrays.toString(b)).append(",\n");
		}
		sb.append("]");
		return sb.toString();
	}
	
	public static final <X> String toStringRecursive(Iterable<X> a){
		if(a==null){return "null";}
		StringBuilder sb=new StringBuilder(256);
		String prefix="";
		sb.append("[");
		for(X x : a){
			sb.append(toStringRecursive(a));
			if(x!=null && x instanceof Iterable<?>){
				sb.append("\n");
			}else{
				sb.append(", ");
			}
		}
		sb.append("]");
		return sb.toString();
	}
	
	public static final HashMap<String, Integer> geneNameToIdTable(){
		if(geneNameToIdTable==null){
			geneIdToNameTable();
			assert(geneIdToNameTable!=null);
			assert(geneNameToIdTable!=null);
		}
		return geneNameToIdTable;
	}
	
	public static final HashMap<Integer, String> geneIdToNameTable(){
		if(geneIdToNameTable==null){

			synchronized(GENEIDLOCK){
				if(geneIdToNameTable==null){

					//			TextFile tf=new TextFile(ROOT_GENE+"gene_names_36.3.txt");
					TextFile tf=new TextFile(ROOT_GENE+"gene_names_37.1.txt", false, false);
					String[] lines=tf.toStringLines();
					tf.close();
					HashMap<Integer, String> table=new HashMap<Integer, String>((lines.length*3)/2);
					for(String s : lines){
						if(!s.startsWith("#")){
							String[] line=s.split("\t", -1);
							//					assert(line.length==3) : "'"+s+"'";
							if(line.length>=3){

								int key=-1;
								try {
									key=Integer.parseInt(line[1]);
								} catch (NumberFormatException e) {
									System.err.println("Bad line: "+s);
									throw new RuntimeException(e);
								}

								table.put(key, (line[2]==null || line[2].length()==0) ? line[1] : line[2]);
							}
						}
					}

					geneIdToNameTable=table;
					
					HashMap<String, Integer> table2=new HashMap<String, Integer>((lines.length*3)/2);
					for(Integer id : geneIdToNameTable.keySet()){
						table2.put(geneIdToNameTable.get(id), id);
					}
					geneNameToIdTable=table2;
				}
			}
		}
		return geneIdToNameTable;
	}
	
	
	public static ChainLine[][] getChainLines(int from, int to){
		if(from==36 && to==37){
			if(chains36to37==null){
				chains36to37=ChainBlock.loadChainLines(ROOT_CHAIN+"hg18ToHg19.over.chain");
			}
			return chains36to37;
		}else if(from==37 && to==36){
			if(chains37to36==null){
				chains37to36=ChainBlock.loadChainLines(ROOT_CHAIN+"hg19ToHg18.over.chain");
			}
			return chains37to36;
		}
		throw new RuntimeException("Unknown chain file: "+from+" -> "+to);
	}
	
	
	public static final String toStringRecursive(Object a){
		return a==null ? "null" : a.toString();
	}
	
	public static boolean isBaited(Variation v){
		return isBaited(v, 0);
	}
	
	public static boolean isBaited(Variation v, int thresh){
		int mid=(v.beginLoc+v.endLoc)/2;
		int len=v.endLoc-v.beginLoc+1;
		return isBaited(v.chromosome, mid, len/2+thresh);
	}
	
	public static boolean isBaited(int chrom, int point, int thresh){
		if(BAITS==null){
			BAITS=(int[][][]) ReadWrite.readObject("UNDEFINED_ROOT"+"baits_"+BAIT_FILE+"_build"+GENOME_BUILD+".int3d");
		}
		return isBaited(point, BAITS[chrom], thresh);
	}
	
	/** Is this point within "thresh" of a bait */
	private static boolean isBaited(int point, int[][] baits, int thresh){

		if(baits==null || baits[0].length==0){return false;}
		
		int[] starts=baits[0];
		int[] stops=baits[1];
		int index=Arrays.binarySearch(stops, point);
		
		if(index>=0){return true;} //Hit inside a bait
		index=(-index)-1;
		
		if(index>=stops.length){return point<=(stops[stops.length-1]+thresh);}
		
//		if(index<0 || index>=stops.length){
//			System.err.println(point+" in "+starts[0]+", "+stops[stops.length-1]+" -> "+index+"/"+(stops.length-1));
//		}

		final int a=point-thresh;
		final int b=point+thresh;
		
		if(overlap(a, b, starts[index], stops[index])){return true;}
		for(int i=index+1; i<starts.length && b>=starts[i]; i++){
			if(overlap(a, b, starts[i], stops[i])){return true;}
		}
		for(int i=index-1; i>=0 && a<=stops[i]; i++){
			if(overlap(a, b, starts[i], stops[i])){return true;}
		}
		return false;
//		
//		return point>=(starts[index]-thresh) && point<=(stops[index]+thresh);
	}
	
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1; 
	}
	
	
	
	/** Note - this does not handle cases where the same instance uses different bait files. */
	public static int[][][] getBaits(){
		if(BAITS==null){
			synchronized(BAITLOCK){
				if(BAITS==null){
					BAITS=(int[][][]) ReadWrite.readObject("UNDEFINED_ROOT"+"baits_"+BAIT_FILE+"_build"+GENOME_BUILD+".int3d");
				}
			}
		}
		assert(BAITS!=null);
		return BAITS;
	}
	
	
	public static final synchronized void setGenome(int g){
		assert(g>0);
		if(genome_set_to==g){return;}
		if(genome_set_to<0){
			setGenome2(g);
		}else{
			throw new RuntimeException("Changing genomes is not currently supported.");
		}
	}
	
	private static final synchronized void setGenome2(int g){
		assert(genome_set_to!=g);
		GENOME_BUILD=g;
		genome_set_to=g;
		numChroms=-1;
		numBases=-1;
		numDefinedBases=-1;
		numContigs=-1;
		numScaffolds=-1;
		name=null;
		genomeSource=null;
		scaffoldPrefixes=false;
		long fastabytes=-1;
		long fastatime=-1;
		final int currentVersion=FastaToChromArrays.currentVersion();
		int version=0;
		
		if(GENOME_BUILD==FastaToChromArrays.LISTBUILD && FastaToChromArrays.SUMMARY_LIST!=null){
			for(int i=0; i<FastaToChromArrays.SUMMARY_LIST.size(); i++){
				final String s=FastaToChromArrays.SUMMARY_LIST.get(i);
				FastaToChromArrays.SUMMARY_LIST.set(i, null);
				if(s.charAt(0)=='#'){
					if(s.startsWith("#Version")){
						String[] split=s.split("\t");
						version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
					}
				}else{
					String[] split=s.split("\t");
					String a=split[0];
					String b=split[1];
					if(a.equalsIgnoreCase("chroms")){numChroms=(int)Long.parseLong(b);}
					else if(a.equalsIgnoreCase("bases")){numBases=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("defined")){numDefinedBases=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("contigs")){numContigs=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("scaffolds")){numScaffolds=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("interpad")){interScaffoldPadding=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("undefined")){}
					else if(a.equalsIgnoreCase("name")){name=b;}
					else if(a.equalsIgnoreCase("source")){genomeSource=b;}
					else if(a.equalsIgnoreCase("bytes")){fastabytes=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("last modified")){fastatime=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("scafprefixes")){scaffoldPrefixes=Tools.parseBoolean(b);}
					else{assert(version<currentVersion) : "In array: Unknown term "+s;}
				}
			}
			FastaToChromArrays.SUMMARY_LIST=null;
		}else{
			String s;
			TextFile tf=new TextFile(ROOT_GENOME+GENOME_BUILD+"/summary.txt", false, false);
			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.charAt(0)=='#'){
					if(s.startsWith("#Version")){
						String[] split=s.split("\t");
						version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
					}
				}else{
					String[] split=s.split("\t");
					String a=split[0];
					String b=split[1];
					if(a.equalsIgnoreCase("chroms")){numChroms=(int)Long.parseLong(b);}
					else if(a.equalsIgnoreCase("bases")){numBases=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("defined")){numDefinedBases=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("contigs")){numContigs=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("scaffolds")){numScaffolds=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("interpad")){interScaffoldPadding=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("undefined")){}
					else if(a.equalsIgnoreCase("name")){name=b;}
					else if(a.equalsIgnoreCase("source")){genomeSource=b;}
					else if(a.equalsIgnoreCase("bytes")){fastabytes=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("last modified")){fastatime=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("scafprefixes")){scaffoldPrefixes=Tools.parseBoolean(b);}
					else{assert(version<currentVersion) : "In file "+tf.name+": Unknown term "+s;}
				}
			}
			tf.close();
		}
		if(numScaffolds==-1){numScaffolds=numChroms;}
		
		if(version<currentVersion){
			assert(false) : "The index format has changed in this version of BBTools.  Please delete the /ref/ directory and re-index from the reference fasta," +
					" or use an older version of BBTools.";
			if(new File(ROOT_GENOME+GENOME_BUILD+"/info.txt").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/info.txt").delete();}
			if(new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt").delete();}
			if(new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz").delete();}
			sysout.println("Regenerating genome info in new format, version "+currentVersion+".");
			dna.FastaToChromArrays.writeInfo(GENOME_BUILD, numChroms, name, genomeSource, true, scaffoldPrefixes);
			genome_set_to=-1;
			setGenome2(g);
			return;
		}

		assert(numChroms>0) : "Genome "+g+": numChroms="+numChroms;
		assert(numBases>0) : "Genome "+g+": numBases="+numBases;
		assert(numDefinedBases>0) : "Genome "+g+": numDefinedBases="+numDefinedBases;
		assert(numBases>=numDefinedBases) : "Genome "+g+": numBases>numDefinedBases : "+numBases+">"+numDefinedBases;
		
		chromosomePlusMatrix=new ChromosomeArray[numChroms+1];
		chromLengths=new int[numChroms+1];
		chromDefinedBases=new int[numChroms+1];
		chromUndefinedBases=new int[numChroms+1];
		chromContigs=new int[numChroms+1];
		chromStartPad=new int[numChroms+1];
		chromScaffolds=new int[numChroms+1];
		
		scaffoldNames=new byte[numChroms+1][][];
		scaffoldLocs=new int[numChroms+1][];
		scaffoldLengths=new int[numChroms+1][];
		
		if(GENOME_BUILD==FastaToChromArrays.LISTBUILD && FastaToChromArrays.INFO_LIST!=null){
			for(int i=0; i<FastaToChromArrays.INFO_LIST.size(); i++){
				final String s=FastaToChromArrays.INFO_LIST.get(i);
				FastaToChromArrays.INFO_LIST.set(i, null);
				if(s.charAt(0)=='#'){
					if(s.startsWith("#Version")){
						String[] split=s.split("\t");
						version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
					}
				}else{
					assert(version==currentVersion);
					String[] split=s.split("\t");
					int chrom=Integer.parseInt(split[0]);
					chromScaffolds[chrom]=Integer.parseInt(split[1]);
					chromContigs[chrom]=(split.length>2 ? Integer.parseInt(split[2]) : -1);
					chromLengths[chrom]=Integer.parseInt(split[3]);
					chromDefinedBases[chrom]=Integer.parseInt(split[4]);
					chromUndefinedBases[chrom]=(split.length>5 ? Integer.parseInt(split[5]) : -1);
					chromStartPad[chrom]=(split.length>6 ? Integer.parseInt(split[6]) : -1);
					//						chromStopPad[chrom]=(split.length>7 ? Integer.parseInt(split[7]) : -1);
					
				}
			}
			FastaToChromArrays.INFO_LIST=null;
		}else{
			String s;
			TextFile tf=new TextFile(ROOT_GENOME+GENOME_BUILD+"/info.txt", false, false);
			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.charAt(0)=='#'){
					if(s.startsWith("#Version")){
						String[] split=s.split("\t");
						version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
					}
				}else{
					
					if(version>=currentVersion){
						String[] split=s.split("\t");
						int chrom=Integer.parseInt(split[0]);
						chromScaffolds[chrom]=Integer.parseInt(split[1]);
						chromContigs[chrom]=(split.length>2 ? Integer.parseInt(split[2]) : -1);
						chromLengths[chrom]=Integer.parseInt(split[3]);
						chromDefinedBases[chrom]=Integer.parseInt(split[4]);
						chromUndefinedBases[chrom]=(split.length>5 ? Integer.parseInt(split[5]) : -1);
						chromStartPad[chrom]=(split.length>6 ? Integer.parseInt(split[6]) : -1);
//						chromStopPad[chrom]=(split.length>7 ? Integer.parseInt(split[7]) : -1);
					}else{
						tf.close();
						if(new File(ROOT_GENOME+GENOME_BUILD+"/info.txt").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/info.txt").delete();}
						if(new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt").delete();}
						if(new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz").delete();}
						sysout.println("Regenerating genome info in new format.");
						dna.FastaToChromArrays.writeInfo(GENOME_BUILD, numChroms, name, genomeSource, true, scaffoldPrefixes);
						tf=new TextFile(ROOT_GENOME+GENOME_BUILD+"/info.txt", false, false);
					}
				}
				
			}
			
			tf.close();
		}
		
		String fname=ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz";
		boolean hasList=(GENOME_BUILD==FastaToChromArrays.LISTBUILD && FastaToChromArrays.SCAF_LIST!=null);
		if(!LOAD_SCAFFOLDS || (!hasList && !new File(fname).exists())){
			for(int i=0; i<scaffoldNames.length; i++){
				scaffoldNames[i]=new byte[][] {("chr"+i).getBytes()};
				scaffoldLocs[i]=new int[] {chromStartPad[i]<0 ? 0 : chromStartPad[i]};
				scaffoldLengths[i]=new int[] {chromLengths[i]};
				if(!LOAD_SCAFFOLDS){chromScaffolds[i]=1;}
				
				assert(chromScaffolds[i]==1) : "This appears to be an old index version.  " +
					"\nPlease regenerate it from the fasta file by rerunning this program,\nusing the ref=<reference file> and overwrite=true flags.\n"+i+", "+chromScaffolds[i];
			}
		}else{
			for(int chrom=0; chrom<scaffoldNames.length; chrom++){
				int num=chromScaffolds[chrom];
				assert(chrom==0 || num>=1) : chrom+", "+num+", "+Arrays.toString(chromScaffolds);
				if(num>0){
					scaffoldNames[chrom]=new byte[num][];
					scaffoldLocs[chrom]=new int[num];
					scaffoldLengths[chrom]=new int[num];
				}
			}
			int[] count=new int[numChroms+1];
			
			if(hasList){
				
				if(verbose){System.err.println("Fetching scaffold names from list:\n\n"+FastaToChromArrays.SCAF_LIST+"\n\n");}
				
				for(int i=0; i<FastaToChromArrays.SCAF_LIST.size(); i++){
					final String s=FastaToChromArrays.SCAF_LIST.get(i);
					
					if(verbose){System.err.println("Processing "+s);}
					
					FastaToChromArrays.SCAF_LIST.set(i, null);
					if(s.charAt(0)=='#'){
						if(verbose){System.err.println("Control string");}
						if(s.startsWith("#Version")){
							assert(version==currentVersion) : "Wrong index version; please delete /ref/genome/\n"+version+", "+currentVersion;
//							String[] split=s.split("\t");
//							version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
//							assert(version==currentVersion) : "Wrong version: "+version+", "+currentVersion;
						}
					}else{
						String[] split=s.split("\t");
						if(verbose){System.err.println("Split into "+Arrays.toString(split));}
						int chrom=Integer.parseInt(split[0]);
						int x=count[chrom];
						count[chrom]++;

						int scaffoldID=Integer.parseInt(split[1]);
						scaffoldLocs[chrom][x]=Integer.parseInt(split[2]);
						scaffoldLengths[chrom][x]=Integer.parseInt(split[3]);
						scaffoldNames[chrom][x]=split[4].getBytes();
						if(verbose){System.err.println("Set scaffoldNames["+chrom+"]["+x+" to "+(scaffoldNames[chrom][x]==null ? "null" : new String(scaffoldNames[chrom][x])));}
					}
				}
				FastaToChromArrays.SCAF_LIST=null;
			}else{

				String s;
				TextFile tf=new TextFile(fname, false, false);
				for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
					if(s.charAt(0)=='#'){
						if(s.startsWith("#Version")){
							assert(version==currentVersion) : "Wrong index version; please delete /ref/genome/\n"+version+", "+currentVersion;
//							String[] split=s.split("\t");
//							version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
//							assert(version==currentVersion) : "Wrong version: "+version+", "+currentVersion;
						}
					}else{
						String[] split=s.split("\t");
						int chrom=Integer.parseInt(split[0]);
						int x=count[chrom];
						count[chrom]++;

						int scaffoldID=Integer.parseInt(split[1]);
						scaffoldLocs[chrom][x]=Integer.parseInt(split[2]);
						scaffoldLengths[chrom][x]=Integer.parseInt(split[3]);
						scaffoldNames[chrom][x]=split[4].getBytes();
					}

				}

				tf.close();
			}
		}
		
//		for(int i=1; i<scaffoldNames.length; i++){
//			stdout.println(Arrays.toString(scaffoldLocs[i]));
//			stdout.println(Arrays.toString(scaffoldLengths[i]));
//			stdout.println(Arrays.toString(scaffoldNames[i]));
//		}
		
	}
	
	
//	public static String contigName(int x){return scaffoldName(x);}
//	
//	public static String scaffoldName(int x){
//		if(scaffoldNames==null){return "chr"+x;}
//		return scaffoldNames[x][0];
//	}
	
	public static HashMap<String, ScafLoc> scafNameTable(){
		
		if(GENOME_BUILD<0){
			assert(scaffoldNameTable==null);
			return null;
		}
		if(scaffoldNameTable!=null){return scaffoldNameTable;}
		synchronized(SCAFMAPLOCK){
			if(scaffoldNameTable!=null){return scaffoldNameTable;}
			scaffoldNameTable=new HashMap<String, ScafLoc>((int)Tools.min(2L*numScaffolds+10, 1000000000));
			for(int chrom=0; chrom<scaffoldNames.length; chrom++){
				if(scaffoldNames[chrom]!=null){
					for(int scafnum=0; scafnum<scaffoldNames[chrom].length; scafnum++){
						byte[] name=scaffoldNames[chrom][scafnum];
						if(name!=null){
							int loc=scaffoldLocs[chrom][scafnum];
							ScafLoc sc=new ScafLoc(new String(name), chrom, loc);
							scaffoldNameTable.put(sc.name, sc);
						}
					}
				}
			}
		}
		return scaffoldNameTable; 
	}


	public static ScafLoc getScafLoc(byte[] name){
		HashMap<String, ScafLoc> map=scafNameTable();
		if(map==null){return null;}
		return map.get(new String(name));
	}
	public static ScafLoc getScafLoc(String name){
		HashMap<String, ScafLoc> map=scafNameTable();
		if(map==null){return null;}
		return map.get(name);
	}

	public static byte[] scaffoldName(int chrom, int loc, int idx){return scaffoldNames[chrom][idx];}
	public static int scaffoldRelativeLoc(int chrom, int loc, int idx){return loc-scaffoldLocs[chrom][idx];}

	public static int scaffoldIndex(int chrom, int loc){
		int[] array=scaffoldLocs[chrom];
		if(array==null || array.length<2){return 0;}
		
		assert(interScaffoldPadding>0);
		loc=loc+interScaffoldPadding/2; //Puts it on closest scaffold if it is between scaffolds
		
		int idx=Arrays.binarySearch(array, loc);
		if(idx>=0){return idx;} //Perfect hit
		
		//Otherwise, return closest scaffold.
		int insertPoint=-1-idx;
		assert(insertPoint>=0 && insertPoint<=array.length);
		int r=max(0, insertPoint-1);
		assert(r>=0 && r<array.length);
		assert(r==0 || loc>array[r]);
		assert(r==array.length-1 || loc<array[r+1]);
		return r;
	}
	
	/** TODO: This can be made faster */
	public static boolean isSingleScaffold(int chrom, int loc1, int loc2){
		assert(loc2>=loc1);
		if(scaffoldLocs==null){return true;}
		assert(chrom>=0 && chrom<scaffoldLocs.length) : chrom+", "+scaffoldLocs.length;
		int[] array=scaffoldLocs[chrom];
		if(array==null || array.length<2){return true;}
		assert(interScaffoldPadding>0);
		
		int idx=Arrays.binarySearch(array, loc1+interScaffoldPadding);
		final int scaf;
		if(idx>=0){scaf=idx;} //Perfect hit
		else{
			int insertPoint=-1-idx;
			assert(insertPoint>=0 && insertPoint<=array.length);
			scaf=max(0, insertPoint-1);
			assert(scaf>=0 && scaf<array.length);
			assert(scaf==0 || loc1+interScaffoldPadding>array[scaf]);
			assert(scaf==array.length-1 || loc1+interScaffoldPadding<array[scaf+1]);
		}
		if(scaf==array.length-1){return true;}
		
		int lowerBound=array[scaf]-interScaffoldPadding;
		int upperBound=array[scaf+1];
		
		if(loc2<lowerBound || loc1>upperBound){return false;} //This could happen if a random read was generated in the start or stop padding.
		assert(scaf==0 || scaf==array.length-1 || (loc1>=lowerBound && loc1<upperBound)) : 
			"chrom="+chrom+", loc1="+loc1+", lowerBound="+lowerBound+", loc2="+loc2+", upperBound="+upperBound;
		return loc2<upperBound;
	}
	
	/** Returns overlap of these two points with the scaffold on which they are centered */
	public static int scaffoldOverlapLength(int chrom, int loc1, int loc2){
		assert(loc2>=loc1);
		int len=loc2-loc1+1;
		if(scaffoldLocs==null){return len;}
		int[] array=scaffoldLocs[chrom];
		if(array==null || array.length<2){return len;}
		assert(interScaffoldPadding>0);
		
		int mid=loc1+(interScaffoldPadding+len)/2;
		int idx=Arrays.binarySearch(array, mid);
		final int scaf;
		if(idx>=0){scaf=idx;} //Perfect hit
		else{
			int insertPoint=-1-idx;
			assert(insertPoint>=0 && insertPoint<=array.length);
			scaf=max(0, insertPoint-1);
			assert(scaf>=0 && scaf<array.length);
			assert(scaf==0 || mid>array[scaf]);
			assert(scaf==array.length-1 || mid<array[scaf+1]) : "\nscaf="+scaf+", array.length="+array.length+"\n"+
				"loc1="+loc1+", loc2="+loc2+", mid="+mid+", interScaffoldPadding="+interScaffoldPadding+", "+array[scaf]+", "+array[scaf+1]+"\n"+
				(loc1+interScaffoldPadding)+", "+array[scaf+1];
		}
		
		int lowerBound=array[scaf];
		int upperBound=lowerBound+scaffoldLengths[chrom][scaf];
//		assert(upperBound==array[scaf+1]) : lowerBound+", "+upperBound+", "+array[scaf+1]+", "+interScaffoldPadding; //This should fail.
		
		return Tools.overlapLength(loc1, loc2, lowerBound, upperBound);
	}
	
	public static void trimScaffoldNames(){
		if(scaffoldNames!=null){
			for(int i=0; i<scaffoldNames.length; i++){
				byte[][] matrix=scaffoldNames[i];
				if(matrix!=null){
					for(int j=0; j<matrix.length; j++){
						byte[] array=matrix[j];
						if(array!=null){
							for(int k=0; k<array.length; k++){
								if(Character.isWhitespace(array[k])){
									matrix[j]=Arrays.copyOf(array, k);
									break;
								}
							}
						}
					}
				}
			}
		}
	}
	
	public static final String findPath(String fname){
		assert(fname!=null);
		if(fname.startsWith("?")){//Look in standard locations
			fname=fname.substring(1);
		}else{//Use this as the literal path
			return fname;
		}
		String path=ROOT+fname;
		boolean vb=false;
		{
			File f=new File(path);
			if(!f.exists()){
				if(vb){System.err.println("Did not find "+fname+" at "+path);}
				f=new File(ROOT);
				String res=f.getParent();
				if(res.length()>0 && !res.endsWith("/")){res=res+"/";}
				res=res+"resources/"+fname;
				f=new File(res);
				if(f.exists()){path=res;}
				else{if(vb){System.err.println("Did not find "+fname+" at "+res);}}
			}
			if(!f.exists()){
				if(vb){System.err.println("Considering fixing "+path+"\n"+path.contains("/file:"));}
				if(path.contains("/file:")){
					String fixed=path.substring(path.lastIndexOf("/file:")+1);
					f=new File(fixed);
					if(f.exists()){path=fixed;}
					else{if(vb){System.err.println("Did not find "+fname+" at "+fixed);}}
				}
			}
			if(!f.exists()){
				if(vb){System.err.println("Considering getResource");}
				URL url=Primes.class.getResource("/"+fname);
				if(url!=null){
					String temp=url.toString();
					if(vb){System.err.println("Found URL "+temp);}
					f=new File(temp);
					//						if(f.exists()){fname=temp;}
					//						else{System.err.println("Did not find "+fname+" at "+temp);}
					path=temp;
				}
			}
			if(!f.exists() && !path.startsWith("jar:")){
				String hardlink="/global/projectb/sandbox/gaag/bbtools/resources/"+fname;
				f=new File(hardlink);
				if(f.exists()){path=hardlink;}
				else{if(vb){System.err.println("Did not find "+fname+" at "+hardlink);}}
			}
			if(!f.exists() && !path.startsWith("jar:")){
				System.err.println("Warning!  Cannot find "+fname+" "+path);
				return null;
			}
		}
		if(vb){System.err.println("Found "+fname+" at "+path);}
		return path;
	}
	
	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}

	public static final byte min(byte x, byte y){return x<y ? x : y;}
	public static final byte max(byte x, byte y){return x>y ? x : y;}
	
	public static final long min(long x, long y){return x<y ? x : y;}
	public static final long max(long x, long y){return x>y ? x : y;}
	
	public static final double min(double x, double y){return x<y ? x : y;}
	public static final double max(double x, double y){return x>y ? x : y;}
	
	public static final float min(float x, float y){return x<y ? x : y;}
	public static final float max(float x, float y){return x>y ? x : y;}
	
	public static int numChroms;
	public static long numBases;
	public static long numDefinedBases;
	public static int numContigs;
	public static int numScaffolds;
	public static int interScaffoldPadding;
	public static int[] chromLengths;
	public static int[] chromDefinedBases;
	public static int[] chromUndefinedBases;
	public static int[] chromContigs;
	public static int[] chromScaffolds;
	public static int[] chromStartPad;
	
	public static byte[][][] scaffoldNames;
	public static int[][] scaffoldLocs;
	/** Does NOT include interScaffoldPadding */
	public static int[][] scaffoldLengths;
	/** Should be true if scaffold names have extra prefixes (for BBSplitter mode), false otherwise */
	public static boolean scaffoldPrefixes;
	
	/** Allows translation of sam coordinates back to native coordinates */
	public static HashMap<String, ScafLoc> scaffoldNameTable;
	
	public static String genomeSource;
	public static String name;
	
	private static final GeneSet[][] geneSetMatrix=new GeneSet[63][];
	private static final Gene[][] geneMatrix=new Gene[63][];
	public static final Range[][] geneSetRangeMatrix=new Range[63][];
	public static final Range[][] geneTxRangeMatrix=new Range[63][];
	public static final Range[][] geneCodeRangeMatrix=new Range[63][];
	private static final Range[][] geneCodeAndExonRangeMatrix=new Range[63][];
	public static final Range[][] exonRangeMatrix=new Range[63][];
	public static HashMap<Integer, ArrayList<GeneSet>> geneIDTable;

	/** Ranges within genes and exons or within NEAR their ends */
	public static final Range[][] geneNearbyRangeMatrix=new Range[63][];

	public static ChromosomeArray[] chromosomePlusMatrix;
	public static ChromosomeArray[] chromosomeMinusMatrix;

	private static HashMap<Integer, String> geneIdToNameTable;
	private static HashMap<String, Integer> geneNameToIdTable;

	private static final HashMap<String, GeneSet>[] geneNameTable=new HashMap[63];
	private static final HashMap<String, Gene>[] transcriptNameTable=new HashMap[63];
	
	public static ChainLine[][] chains36to37;
	public static ChainLine[][] chains37to36;

	public static int[][][] BAITS;

	private static final int TX_RANGE=0;
	private static final int CODE_RANGE=1;
	private static final int CODEEXON_RANGE=2;
	private static final int EXON_RANGE=3;
	private static final int NEAR_RANGE=4;
	
	
	public static final int NEAR=200;

	public static boolean ENV=(System.getenv()!=null);
	public static boolean WINDOWS=(System.getenv().containsKey("OS") && System.getenv().get("OS").equalsIgnoreCase("Windows_NT"));
	public static int LOGICAL_PROCESSORS=Runtime.getRuntime().availableProcessors();
	private static String HOSTNAME;
	public static String HOSTNAME(){
		if(HOSTNAME==null){
			try {
				java.net.InetAddress localMachine = java.net.InetAddress.getLocalHost();
				HOSTNAME=localMachine.getHostName();
			} catch (UnknownHostException e) {
				// TODO Auto-generated catch block
//				e.printStackTrace();
				HOSTNAME="unknown";
			} catch (NullPointerException e) {
				// TODO Auto-generated catch block
//				e.printStackTrace();
				HOSTNAME="unknown";
			} catch (Throwable e) {
				HOSTNAME="unknown";
			}
		}
		return HOSTNAME;
	}
	
	
	/** Should be the same as ROOT_BASE but is found dynamically */
	public static final String ROOT;
	
	public static String ROOT_BASE;
	public static String ROOT_OUTPUT;
	public static String ROOT_REF;
	public static String ROOT_GENOME;
	public static String ROOT_INDEX;
	public static String ROOT_GENE;
	public static String ROOT_CHAIN;
	public static String ROOT_AFFY;
	public static String ROOT_TEMPDIR;
	
	static{
		ROOT=(new File(Data.class.getClassLoader().getResource(Data.class.getName().replace('.', '/') + ".class")
				.getFile()).getAbsolutePath().replace('\\', '/').replace("dna/Data.class", ""));
		setPath(WINDOWS ? "windows" : "unix");
		if(!WINDOWS || true){setPath("local");}
	}
	
	public static void setPath(String path){
//		System.err.println("***"+path);
		if(path.indexOf('\\')>=0){path=path.replace('\\', '/');}
		String mode=(path==null ? "null" : path.toLowerCase());
		boolean local=mode.equals("local") || mode.equals(".") || mode.equals("/.") || mode.equals("./");
		boolean win=mode.contains("windows");
		boolean unix=mode.contains("unix");
		
		if(local){
			ROOT_BASE="";
			ROOT_OUTPUT="";
			ROOT_REF="ref/";
			ROOT_GENOME=ROOT_REF+"genome/";
			ROOT_INDEX=ROOT_REF+"index/";
//			ROOT_GENE="D:/Data/ref/genes/";
//			ROOT_CHAIN="D:/Data/ref/chain/";
//			ROOT_AFFY="D:/Data/ref/affy/";
//			ROOT_TEMPDIR="";
		}else if(win){
			ROOT_BASE="C:/workspace/prune/";
			ROOT_OUTPUT="C:/workspace/";
			ROOT_REF="D:/Data/ref/";
			ROOT_GENOME=ROOT_REF+"genome/";
			ROOT_INDEX=ROOT_REF+"index/";
			ROOT_GENE="D:/Data/ref/genes/";
			ROOT_CHAIN="D:/Data/ref/chain/";
			ROOT_AFFY="D:/Data/ref/affy/";
			ROOT_TEMPDIR="C:/workspace/tempdir/";
		}else if(unix){
			ROOT_BASE="/house/homedirs/b/bushnell/prune/";
			ROOT_OUTPUT="/house/homedirs/b/bushnell/prune/";
			ROOT_REF="/house/homedirs/b/bushnell/Data/ref/";
			ROOT_GENOME=ROOT_REF+"genome/";
			ROOT_INDEX=ROOT_REF+"index/";
			ROOT_GENE="/house/homedirs/b/bushnell/Data/ref/genes/";
			ROOT_CHAIN="/house/homedirs/b/bushnell/Data/ref/chain/";
			ROOT_AFFY="/house/homedirs/b/bushnell/Data/ref/affy/";
			String s=System.getenv().get("TEMPDIR");
			ROOT_TEMPDIR=(s==null ? "" : s+"/");
		}else if(!"null".equals(mode)){
			if(!path.endsWith("/")){path=path+"/";}
			ROOT_BASE=path;
			ROOT_OUTPUT=path;
			ROOT_REF=path+"ref/";
			ROOT_GENOME=ROOT_REF+"genome/";
			ROOT_INDEX=ROOT_REF+"index/";
//			ROOT_GENE="D:/Data/ref/genes/";
//			ROOT_CHAIN="D:/Data/ref/chain/";
//			ROOT_AFFY="D:/Data/ref/affy/";
//			ROOT_TEMPDIR="";
		}else{
			ROOT_BASE=null;
			ROOT_OUTPUT=null;
			ROOT_REF=null;
			ROOT_GENOME=null;
			ROOT_GENE=null;
			ROOT_CHAIN=null;
			ROOT_AFFY=null;
			ROOT_TEMPDIR=null;
		}
	}
	
	public static final String VAR_FOLDER="VAR/";
	public static final String GENE_FOLDER="GENE/";
	
	public static int GENOME_BUILD=-1;
	private static int genome_set_to=-1;
	public static int DBSNP_BUILD=132;
	
	public static final byte CGI=0;
	public static final byte CORE_SOLID=1;
	public static final byte CORE_ILLUMINA=2;
	public static final byte OTOGEN_ILLUMINA=3;
	
	public static final boolean verbose=false;
	
	/** seqGene, knownGene, refGene, unionGene, seqRefGene, ccs */
	public static String GENE_MAP="seqRefGene"; 
	
	/** 30M, 50M, nimblegen, truseq_60M */
	public static String BAIT_FILE="truseq_60M"; 
	
	private static final String BAITLOCK=new String("BAITLOCK");
	private static final String GENEIDLOCK=new String("GENEIDLOCK");
	
	private static final String[] CHROMLOCKS=new String[256];
	
	static{
		for(int i=0; i<CHROMLOCKS.length; i++){
			CHROMLOCKS[i]=new String(i+"");
		}
	}
	
	private static final int INTERN_MAP_SIZE=(1<<20);
	private static final int INTERN_MAP_LIMIT=(1<<19);
	
	private static final HashMap<String, String> INTERNMAP=new HashMap<String, String>(INTERN_MAP_SIZE);
//	public static final void unloadInternMap(){
//		INTERNMAP=new HashMap<String, String>(INTERN_MAP_SIZE);
//	}
	
	private static String condense(String s){
		//TODO - finish this
		StringBuilder sb=new StringBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			sb.append('A');
		}
		return sb.toString();
	}
	
	public static final void intern(String[] s){
		if(s==null){return;}
		for(int i=0; i<s.length; i++){s[i]=intern(s[i]);}
	}
	public static String intern(String s){
		if(s==null || s.length()>25){return new String(s);}
//		calls++;
//		
//		if(s.length()>0 && s.charAt(0)!='?'){
//			s=condense(s);
//		}
		
		if(s.length()<2){
//			return s.intern();
			return forceIntern(s);
		}
		boolean acgtn=AminoAcid.containsOnlyACGTNQ(s);
		
		if(acgtn){
			if(s.length()<4){
//				return s.intern();
				return forceIntern(s);
			}
			if(s.length()>6){
				return new String(s);
			}
		}
		
		//Otherwise it is non-base string of length 2 to 20, or a base string of length 4 to 6.
		return forceIntern(s);
	}
	
	public static String forceIntern(String s){
		calls++;
		
//		if(s.length()<2){return s.intern();}
//		boolean acgtn=AminoAcid.containsOnlyACGTNQ(s);
//		
//		if(acgtn){
//			if(s.length()<4){return s.intern();}
//		}
		
		String old=INTERNMAP.get(s);
		if(old!=null){return old;}
		
		synchronized(INTERNMAP){
//			System.err.print(INTERNMAP.size()+"~"+calls+": "+s+", ");
			if(INTERNMAP.size()>INTERN_MAP_LIMIT){
				System.err.println("INTERNMAP overflow caused by "+s);
				INTERNMAP.clear();
			}
			if(INTERNMAP.containsKey(s)){return INTERNMAP.get(s);}
			s=new String(s);
			INTERNMAP.put(s, s);
		}
		return s;
	}
	static int calls=0;
	
	public static PrintStream sysout=System.err;//System.out;
	
	public static boolean CHROMC=false;
	public static boolean CHROMGZ=true;
	public static boolean LOAD_SCAFFOLDS=true;

//	private static final boolean GUNZIP=testExecute("gunzip --help");
//	private static final boolean GZIP=testExecute("gzip --help");
//	private static final boolean SAMTOOLS=testExecute("samtools --help");

	public static boolean GUNZIP(){return GUNZIP==0 ? GZIP() : GUNZIP>0;}
	public static boolean UNPIGZ(){return UNPIGZ==0 ? PIGZ() : UNPIGZ>0;}
	public static boolean GZIP(){
		if(GZIP==0 && !WINDOWS){
			synchronized(SUBPROCSYNC){
				if(GZIP==0){GZIP=testExecute("gzip --version");}
			}
		} 
		return GZIP>0;
	}
	public static boolean PIGZ(){
		if(PIGZ==0 && !WINDOWS){
			synchronized(SUBPROCSYNC){
				if(PIGZ==0){PIGZ=testExecute("pigz --version");}
			}
		} 
		return PIGZ>0;
	}
	public static boolean SAMTOOLS(){
		if(SAMTOOLS==0 && !WINDOWS){
			synchronized(SUBPROCSYNC){
				if(SAMTOOLS==0){SAMTOOLS=testExecute("samtools");}
			}
			System.err.println(SAMTOOLS>0 ? "Found samtools." : "Could not find samtools.");
		}
		return SAMTOOLS>0;
	}
	public static boolean SH(){
		if(SH==0 && !WINDOWS){
			synchronized(SUBPROCSYNC){
				if(SH==0){SH=testExecute("sh --version");}
			}
//			System.err.println(SH>0 ? "Found sh." : "Could not find sh.");
			if(SH<0){System.err.println("Could not find sh; won't launch I/O subprocesses.");}
		}
		return SH>0;
	}
	private static final String SUBPROCSYNC=new String("SUBPROCSYNC");
	private static final String SCAFMAPLOCK=new String("SCAFMAPLOCK");
	
	/* Set these to zero to enable or -1 to disable */
	private static int GUNZIP=-1;
	private static int UNPIGZ=0;
	private static int GZIP=0;
	private static int PIGZ=0;
	private static int SAMTOOLS=0;
	private static int SH=0;
	
	private static int testExecute(String s){
//		System.err.println("Testing "+s);
		try {
			Process p;
			p = Runtime.getRuntime().exec(s);
//			System.err.println("Got process.");
			while(p.getErrorStream().read()>-1){}
//			return p.exitValue()==0;
//			System.err.println("This system does has "+s+" installed.");
		} catch (IOException e) {
//			System.err.println("This system does not have "+s+" installed.");
			// TODO Auto-generated catch block
//			e.printStackTrace();
			return -1;
		}
		return 1;
	}
	
}
