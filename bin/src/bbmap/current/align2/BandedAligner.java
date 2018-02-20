package align2;

import java.util.Arrays;

import dna.AminoAcid;

/**
 * @author Brian Bushnell
 * @date Aug 5, 2013
 *
 */
public class BandedAligner {
	
	
	public static void main(String[] args){
		byte[] query=args[0].getBytes();
		byte[] ref=args[1].getBytes();
		int qstart=-1;
		int rstart=-1;
		int maxedits=big;
		int width=5;
		if(args.length>2){qstart=Integer.parseInt(args[2]);}
		if(args.length>3){rstart=Integer.parseInt(args[3]);}
		if(args.length>4){maxedits=Integer.parseInt(args[4]);}
		if(args.length>4){width=Integer.parseInt(args[5]);}
		
		BandedAligner ba=new BandedAligner(width);
		
		int edits;
		
		edits=ba.alignForward(query, ref, (qstart==-1 ? 0 : qstart), (rstart==-1 ? 0 : rstart), maxedits, true);
		System.out.println("Forward:    \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");
//		
//		edits=ba.alignForwardRC(query, ref, (qstart==-1 ? query.length-1 : qstart), (rstart==-1 ? 0 : rstart), maxedits, true);
//		System.out.println("ForwardRC:  \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
//		System.out.println("***********************\n");
		
		edits=ba.alignReverse(query, ref, (qstart==-1 ? query.length-1 : qstart), (rstart==-1 ? ref.length-1 : rstart), maxedits, true);
		System.out.println("Reverse:    \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");
		
//		edits=ba.alignReverseRC(query, ref, (qstart==-1 ? 0 : qstart), (rstart==-1 ? ref.length-1 : rstart), maxedits, true);
//		System.out.println("ReverseRC:  \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
//		System.out.println("***********************\n");
	}
	
	
	public BandedAligner(int width_){
		maxWidth=Tools.max(width_, 3)|1;
		assert(maxWidth>=3) : "width<3 : "+width_+" -> "+maxWidth;
		array1=new int[maxWidth+2];
		array2=new int[maxWidth+2];
		Arrays.fill(array1, big);
		Arrays.fill(array2, big);
//		for(int i=2; i<rows; i++){
//			matrix[i]=matrix[i-2];
//		}
		assert(big>maxWidth/2);
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public int alignForward(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		assert(big>maxEdits);
		if(verbose){System.err.println("alignForward("+new String(query)+", "+new String(ref)+", "+qstart+", "+rstart+", "+maxEdits+")");}
		if(query.length-qstart>ref.length-rstart){
			int x=alignForward(ref, query, rstart, qstart, maxEdits, exact);
			int temp=lastQueryLoc;
			lastQueryLoc=lastRefLoc;
			lastRefLoc=temp;
			if(verbose){
				System.out.println("Reversed.");
				System.out.println("Final state: lastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
						", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
			}
			return x;
		}
		int edits=0, row=0;
		lastRow=-1;
		lastEdits=0;
		lastOffset=0;
		
		final int width=Tools.min(maxWidth, (maxEdits*2)+1);
		final int halfWidth=width/2;
		final boolean inexact=!exact;
		
		int qloc=qstart;
		int rsloc=rstart-halfWidth;
		final int xlines=query.length-qstart;
		final int ylines=ref.length-rstart;
		final int len=Tools.min(xlines, ylines);
		if(verbose){System.err.println("xlines="+xlines+", ylines="+ylines+", len="+len);}
		if(len<1){
			if(false){
				throw new RuntimeException("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			}
			assert(false) : ("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			return 0;
		}

		Arrays.fill(array1, big);
		Arrays.fill(array2, big);
		arrayCurrent=array1;
		arrayPrev=array2;
		
		{
			if(verbose){System.err.println("\nFirst row.");}
			final byte q=query[qloc];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			edits=big;
			int mloc=1+(colStart-rsloc);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
//			assert(false) : mloc+", "+colStart+", "+rsloc;
			for(int col=colStart; col<colLimit; mloc++, col++){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int score=(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; prev=0; score="+score+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			row++; qloc++; rsloc++;
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}
		
		for(row=1; row<len; row++, qloc++, rsloc++){
//			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayCurrent));}
			arrayTemp=arrayCurrent;
			arrayCurrent=arrayPrev;
			arrayPrev=arrayTemp;
			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayPrev)+", current="+Arrays.toString(arrayCurrent));}
			final byte q=query[qloc];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			Arrays.fill(arrayCurrent, big);
			edits=big;
			int mloc=1+(colStart-rsloc);
			boolean forceDiag=(row==len-1);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colStart; col<colLimit; mloc++, col++){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int scoreUp=arrayPrev[mloc+1]+1;
				final int scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				final int scoreLeft=arrayCurrent[mloc-1]+1;
				final int score=(forceDiag || col==ref.length-1) ? scoreDiag : Tools.min(scoreUp, scoreDiag, scoreLeft);
				if(verbose){System.err.println("prev=min(s["+(mloc-1)+"]="+arrayCurrent[mloc-1]+", p["+(mloc)+"]="+arrayPrev[mloc]+", p["+(mloc+1)+"]="+arrayPrev[mloc+1]+")");}
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; up="+scoreUp+"; diag="+scoreDiag+"; left="+scoreLeft+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			if(edits>maxEdits){row++; break;}
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}

		lastRow=row-1;
		lastEdits=edits;
		lastQueryLoc=qloc-1;
		lastOffset=lastOffset(arrayCurrent, halfWidth);
		lastRefLoc=rsloc+halfWidth-lastOffset-1;
		while(lastRefLoc>=ref.length || lastQueryLoc>=query.length){lastRefLoc--; lastQueryLoc--;}
		if(verbose){
			System.out.println("\nFinal state: arrayCurrent="+Arrays.toString(arrayCurrent)+"\nlastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
					", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
		}
		return edits;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public int alignForwardRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		assert(big>maxEdits);
		if(verbose){System.err.println("alignForwardRC("+new String(query)+", "+new String(ref)+", "+qstart+", "+rstart+", "+maxEdits+")");}
		if(qstart+1>ref.length-rstart){
			int x=alignReverseRC(ref, query, rstart, qstart, maxEdits, exact);
			int temp=lastQueryLoc;
			lastQueryLoc=lastRefLoc;
			lastRefLoc=temp;
			if(verbose){
				System.out.println("Reversed.");
				System.out.println("Final state: lastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
						", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
			}
			return x;
		}
		int edits=0, row=0;
		lastRow=-1;
		lastEdits=0;
		lastOffset=0;
		
		final int width=Tools.min(maxWidth, (maxEdits*2)+1);
		final int halfWidth=width/2;
		final boolean inexact=!exact;
		
		int qloc=qstart;
		int rsloc=rstart-halfWidth;
		final int xlines=qstart+1;
		final int ylines=ref.length-rstart;
		final int len=Tools.min(xlines, ylines);
		if(verbose){System.err.println("xlines="+xlines+", ylines="+ylines+", len="+len);}
		if(len<1){
			if(false){
				throw new RuntimeException("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			}
			assert(false) : ("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			return 0;
		}

		Arrays.fill(array1, big);
		Arrays.fill(array2, big);
		arrayCurrent=array1;
		arrayPrev=array2;
		
		{
			if(verbose){System.err.println("\nFirst row.");}
			final byte q=AminoAcid.baseToComplementExtended[query[qloc]];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			edits=big;
			int mloc=1+(colStart-rsloc);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colStart; col<colLimit; mloc++, col++){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int score=(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; prev=0; score="+score+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			row++; qloc--; rsloc++;
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}
		
		for(row=1; row<len; row++, qloc--, rsloc++){
//			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayCurrent));}
			arrayTemp=arrayCurrent;
			arrayCurrent=arrayPrev;
			arrayPrev=arrayTemp;
			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayPrev)+", current="+Arrays.toString(arrayCurrent));}
			final byte q=AminoAcid.baseToComplementExtended[query[qloc]];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			Arrays.fill(arrayCurrent, big);
			edits=big;
			int mloc=1+(colStart-rsloc);
			boolean forceDiag=(row==len-1);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colStart; col<colLimit; mloc++, col++){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int scoreUp=arrayPrev[mloc+1]+1;
				final int scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				final int scoreLeft=arrayCurrent[mloc-1]+1;
				final int score=(forceDiag || col==ref.length-1) ? scoreDiag : Tools.min(scoreUp, scoreDiag, scoreLeft);
				if(verbose){System.err.println("prev=min(s["+(mloc-1)+"]="+arrayCurrent[mloc-1]+", p["+(mloc)+"]="+arrayPrev[mloc]+", p["+(mloc+1)+"]="+arrayPrev[mloc+1]+")");}
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; up="+scoreUp+"; diag="+scoreDiag+"; left="+scoreLeft+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			if(edits>maxEdits){row++; break;}
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}

		lastRow=row-1;
		lastEdits=edits;
		lastOffset=lastOffset(arrayCurrent, halfWidth);
		lastQueryLoc=qloc+1;
		lastRefLoc=rsloc+halfWidth-lastOffset-1;
		while(lastRefLoc>=ref.length || lastQueryLoc<0){lastRefLoc--; lastQueryLoc++;}
		if(verbose){
			System.out.println("\nFinal state: arrayCurrent="+Arrays.toString(arrayCurrent)+"\nlastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
					", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+", qloc="+qloc+"\n");
		}
		return edits;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public int alignReverse(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		assert(big>maxEdits);
		if(verbose){System.err.println("alignReverse("+new String(query)+", "+new String(ref)+", "+qstart+", "+rstart+", "+maxEdits+")");}
		if(qstart>rstart){
			int x=alignReverse(ref, query, rstart, qstart, maxEdits, exact);
			int temp=lastQueryLoc;
			lastQueryLoc=lastRefLoc;
			lastRefLoc=temp;
			if(verbose){
				System.out.println("Reversed.");
				System.out.println("Final state: lastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
						", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
			}
			return x;
		}
//		if(true){return big;}
		int edits=0, row=0;
		lastRow=-1;
		lastEdits=0;
		lastOffset=0;
		
		final int width=Tools.min(maxWidth, (maxEdits*2)+1);
		final int halfWidth=width/2;
		final boolean inexact=!exact;
		
		int qloc=qstart;
		int rsloc=rstart-halfWidth;
		final int xlines=qstart+1;
		final int ylines=rstart+1;
		final int len=Tools.min(xlines, ylines);
		if(verbose){System.err.println("xlines="+xlines+", ylines="+ylines+", len="+len);}
		if(len<1){
			if(false){
				throw new RuntimeException("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			}
			assert(false) : ("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			return 0;
		}

		Arrays.fill(array1, big);
		Arrays.fill(array2, big);
		arrayCurrent=array1;
		arrayPrev=array2;
		
		{
			if(verbose){System.err.println("\nFirst row.");}
			final byte q=query[qloc];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			edits=big;
			int mloc=1+width-(colLimit-rsloc);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colLimit-1; col>=colStart; mloc++, col--){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int score=(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; prev=0; score="+score+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			row++; qloc--; rsloc--;
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}
		
		for(row=1; row<len; row++, qloc--, rsloc--){
//			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayCurrent));}
			arrayTemp=arrayCurrent;
			arrayCurrent=arrayPrev;
			arrayPrev=arrayTemp;
			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayPrev)+", current="+Arrays.toString(arrayCurrent));}
			final byte q=query[qloc];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			Arrays.fill(arrayCurrent, big);
			edits=big;
			int mloc=1+width-(colLimit-rsloc);
			boolean forceDiag=(row==len-1);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colLimit-1; col>=colStart; mloc++, col--){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int scoreUp=arrayPrev[mloc+1]+1;
				final int scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				final int scoreLeft=arrayCurrent[mloc-1]+1;
				final int score=(forceDiag || col==0) ? scoreDiag : Tools.min(scoreUp, scoreDiag, scoreLeft);
				if(verbose){System.err.println("prev=min(s["+(mloc-1)+"]="+arrayCurrent[mloc-1]+", p["+(mloc)+"]="+arrayPrev[mloc]+", p["+(mloc+1)+"]="+arrayPrev[mloc+1]+")");}
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; up="+scoreUp+"; diag="+scoreDiag+"; left="+scoreLeft+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			if(edits>maxEdits){row++; break;}
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}

		lastRow=row-1;
		lastEdits=edits;
		lastOffset=lastOffset(arrayCurrent, halfWidth);
		lastQueryLoc=qloc+1;
		lastRefLoc=rsloc+halfWidth+lastOffset+1;
		while(lastRefLoc<0 || lastQueryLoc<0){lastRefLoc++; lastQueryLoc++;}
		if(verbose){
			System.out.println("\nFinal state: arrayCurrent="+Arrays.toString(arrayCurrent)+"\nlastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
					", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+", qloc="+qloc+", rsloc="+rsloc+"\n");
		}
		return edits;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public int alignReverseRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		assert(big>maxEdits);
		if(verbose){System.err.println("alignReverseRC("+new String(query)+", "+new String(ref)+", "+qstart+", "+rstart+", "+maxEdits+")");}
		if(query.length-qstart>rstart+1){
			int x=alignForwardRC(ref, query, rstart, qstart, maxEdits, exact);
			int temp=lastQueryLoc;
			lastQueryLoc=lastRefLoc;
			lastRefLoc=temp;
			if(verbose){
				System.out.println("Reversed.");
				System.out.println("Final state: lastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
						", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
			}
			return x;
		}
		int edits=0, row=0;
		lastRow=-1;
		lastEdits=0;
		lastOffset=0;
		
		final int width=Tools.min(maxWidth, (maxEdits*2)+1);
		final int halfWidth=width/2;
		final boolean inexact=!exact;
		
		int qloc=qstart;
		int rsloc=rstart-halfWidth;
		final int xlines=query.length-qstart;
		final int ylines=rstart+1;
		final int len=Tools.min(xlines, ylines);
		if(verbose){System.err.println("xlines="+xlines+", ylines="+ylines+", len="+len);}
		if(len<1){
			if(false){
				throw new RuntimeException("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			}
			assert(false) : ("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			return 0;
		}

		Arrays.fill(array1, big);
		Arrays.fill(array2, big);
		arrayCurrent=array1;
		arrayPrev=array2;
		
		{
			if(verbose){System.err.println("\nFirst row.");}
			final byte q=AminoAcid.baseToComplementExtended[query[qloc]];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			edits=big;
			int mloc=1+width-(colLimit-rsloc);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colLimit-1; col>=colStart; mloc++, col--){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int score=(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; prev=0; score="+score+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			row++; qloc++; rsloc--;
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}
		
		for(row=1; row<len; row++, qloc++, rsloc--){
//			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayCurrent));}
			arrayTemp=arrayCurrent;
			arrayCurrent=arrayPrev;
			arrayPrev=arrayTemp;
			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayPrev)+", current="+Arrays.toString(arrayCurrent));}
			final byte q=AminoAcid.baseToComplementExtended[query[qloc]];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			Arrays.fill(arrayCurrent, big);
			edits=big;
			int mloc=1+width-(colLimit-rsloc);
			boolean forceDiag=(row==len-1);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colLimit-1; col>=colStart; mloc++, col--){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int scoreUp=arrayPrev[mloc+1]+1;
				final int scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				final int scoreLeft=arrayCurrent[mloc-1]+1;
				final int score=(forceDiag || col==0) ? scoreDiag : Tools.min(scoreUp, scoreDiag, scoreLeft);
				if(verbose){System.err.println("prev=min(s["+(mloc-1)+"]="+arrayCurrent[mloc-1]+", p["+(mloc)+"]="+arrayPrev[mloc]+", p["+(mloc+1)+"]="+arrayPrev[mloc+1]+")");}
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; up="+scoreUp+"; diag="+scoreDiag+"; left="+scoreLeft+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			if(edits>maxEdits){row++; break;}
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}

		lastRow=row-1;
		lastEdits=edits;
		lastOffset=lastOffset(arrayCurrent, halfWidth);
		lastQueryLoc=qloc-1;
		lastRefLoc=rsloc+halfWidth+lastOffset+1;
		while(lastRefLoc<0 || lastQueryLoc>=query.length){lastRefLoc++; lastQueryLoc--;}
		if(verbose){
			System.out.println("\nFinal state: arrayCurrent="+Arrays.toString(arrayCurrent)+"\nlastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
					", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
		}
		return edits;
	}
	
	private static final void fillBig(int[] array){
		final int lim=array.length-1;
		for(int i=1; i<lim; i++){array[i]=big;}
	}
	
	/** Score is lastRow-edits */
	public int score(){
		return lastRow-lastEdits+1;
	}
	
	/** Position of min value in array (meaning the best alignment) relative to the middle of the array. */
	private int lastOffset(int[] array, int halfWidth){
		final int center=halfWidth+1;
		int minLoc=center;
		for(int i=1; i<=halfWidth; i++){
			if(array[center+i]<array[minLoc]){minLoc=center+i;}
			if(array[center-i]<array[minLoc]){minLoc=center-i;}
		}
		return center-minLoc;
	}
	
	private static final int penalizeOffCenter(int[] array, int halfWidth){
		final int center=halfWidth+1;
		int edits=array[center];
		for(int i=1; i<=halfWidth; i++){
			array[center+i]=Tools.min(big, array[center+i]+i);
			edits=Tools.min(edits, array[center+i]);
			array[center-i]=Tools.min(big, array[center-i]+i);
			edits=Tools.min(edits, array[center-i]);
		}
		return edits;
	}
	
	/** Final row aligned in last alignment. */
	public int lastRow;
	/** Final edits value in last alignment. */
	public int lastEdits;

	/** Position of min value in array (meaning the best alignment) relative to the middle of the array.
	 * Positive value is to the right (ref sequence longer than query), negative value left (ref shorter than query) */
	private int lastOffset;
	
	public int lastRefLoc;
	public int lastQueryLoc;
	
//	public int[][] matrix;
	private final int[] array1;
	private final int[] array2;
	private int[] arrayCurrent, arrayPrev, arrayTemp;
	public final int maxWidth;

	private static final int big=999;
	public static boolean verbose=false;
	/** Penalizes non-length-neutral alignments.  
	 * This Causes query-to-ref alignment to yield same score as ref-to-query alignment, which is useful for assertions.  */ 
	public static boolean penalizeOffCenter=true;
	
}
