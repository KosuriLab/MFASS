package driver;

import java.util.HashMap;

import dna.Data;
import dna.Exon;
import dna.Gene;
import dna.GeneSet;

public class CountNumberOfBaitedGenes {
	
	
	public static void main(String[] args){
		
		Data.GENOME_BUILD=37;
		Data.GENE_MAP="refGene";
		
		int genes=0;
		int partly=0;
		int fully=0;
		
		Data.getBaits();
		
		for(byte chrom=1; chrom<=24; chrom++){
			HashMap<String, GeneSet> table=Data.geneNameTable(chrom);
			
			
			
			for(GeneSet gs : table.values()){
				boolean fullybaited=false;
				boolean partlybaited=false;
				boolean exists=false;
				for(Gene g : gs.genes){
					if(!g.untranslated && g.codeLength>1){
						exists=true;
						fullybaited=(fullybaited||isFullyBaited(g));
						partlybaited=(partlybaited||isPartlyBaited(g));
					}
					if(fullybaited){break;}
				}
				if(fullybaited){fully++;}
				if(partlybaited){partly++;}
				if(exists){genes++;}
			}
			System.err.print(".");
		}
		System.err.println();
		System.err.println("Genes:         \t"+genes);
		System.err.println("Fully Baited:  \t"+fully);
		System.err.println("Partly Baited: \t"+partly);
		
	}
	
	
	public static boolean isFullyBaited(Gene g){
		
		int count=0;
		
		for(Exon e : g.exons){
			for(int i=e.a; i<=e.b; i++){
				if(i>=g.codeStart && i<=g.codeStop){
					if(!Data.isBaited(g.chromosome, i, 0)){return false;}
					count++;
				}
			}
		}
		
		if(count<1){return false;} //No content
		return true;
	}
	
	
	public static boolean isPartlyBaited(Gene g){
		
		int count=0;
		
		for(Exon e : g.exons){
			for(int i=e.a; i<=e.b; i++){
				if(i>=g.codeStart && i<=g.codeStop){
					if(Data.isBaited(g.chromosome, i, 0)){return true;}
					count++;
				}
			}
		}
		
		if(count<1){return false;} //No content
		return false;
	}
	
	
}
