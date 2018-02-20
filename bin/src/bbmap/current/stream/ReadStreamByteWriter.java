package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;

public class ReadStreamByteWriter extends ReadStreamWriter {

	public ReadStreamByteWriter(String fname_, boolean read1_, int bufferSize, boolean allowSubprocess_){
		this(fname_, null, read1_, bufferSize, false, false, false, false, false, false, false, false, allowSubprocess_);
	}
	
//	public ReadStreamByteWriter(String fname_, boolean read1_, int bufferSize, 
//			boolean outputSamFile, boolean outputBamFile, boolean fastq, boolean fasta, boolean sitesOnly, boolean attachment, boolean stdout){
//		this(fname_, null, read1_, bufferSize, outputSamFile, outputBamFile, fastq, fasta, sitesOnly, attachment, stdout, false);
//	}

	public ReadStreamByteWriter(String fname_, String qfname_, boolean read1_, int bufferSize, 
			boolean outputSamFile, boolean outputBamFile, boolean fastq, boolean fasta, boolean sitesOnly, boolean attachment, boolean stdout,
			boolean useSharedHeader, boolean allowSubprocess_){
		super(fname_, qfname_, read1_, bufferSize, 
				outputSamFile, outputBamFile, fastq, fasta, sitesOnly, attachment, stdout,
				useSharedHeader, false, buffered, allowSubprocess_);
	}

	public ReadStreamByteWriter(FileFormat ff, String qfname_, boolean read1_, int bufferSize, CharSequence header, boolean useSharedHeader){
		super(ff, qfname_, read1_, bufferSize, header, false, buffered, useSharedHeader);
	}
	
	@Override
	public void run() {
		try {
			run2();
		} catch (IOException e) {
			finishedSuccessfully=false;
//			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
	
	public void run2() throws IOException {
		
		if(!OUTPUT_SAM && !OUTPUT_FASTQ && !OUTPUT_FASTA && !OUTPUT_ATTACHMENT){
			if(OUTPUT_INTERLEAVED){
//				assert(false) : OUTPUT_SAM+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT+", "+OUTPUT_INTERLEAVED+", "+SITES_ONLY;
				myOutstream.write("#INTERLEAVED\n".getBytes());
			}
			if(SITES_ONLY){
				myOutstream.write(("#"+SiteScore.header()+"\n").getBytes());
			}else if(!OUTPUT_ATTACHMENT){
				myOutstream.write(("#"+Read.header()+"\n").getBytes());
			}
		}
		
		Job job=null;
		while(job==null){
			try {
				job=queue.take();
//				job.list=queue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		final ByteBuilder bb=new ByteBuilder(65000);
		final ByteBuilder bbq=(myQOutstream==null ? null : new ByteBuilder(65000));
		
		while(job!=null && !job.poison){

			final OutputStream abd=job.outstream;
			final OutputStream abc=myOutstream;
			
			if(!job.isEmpty()){
				
				if(myQOutstream!=null){
					bbq.setLength(0);
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								{
									bbq.append('\n');
									bbq.append('>');
									bbq.append(r.id);
									bbq.append('\n');
									if(r.bases!=null){toQualityB(r.quality, r.bases.length, bbq);}
									bbq.append('\n');
								}
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									bbq.append('\n');
									bbq.append('>');
									bbq.append(r2.id);
									bbq.append('\n');
									if(r2.bases!=null){toQualityB(r2.quality, r2.bases.length, bbq);}
									bbq.append('\n');
								}
							}
							if(bbq.length>=32768){
								myQOutstream.write(bbq.array, 0, bbq.length);
								bbq.setLength(0);
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								bbq.append('\n');
								bbq.append('>');
								bbq.append(r2.id);
								bbq.append('\n');
								if(r2.bases!=null){toQualityB(r2.quality, r2.bases.length, bbq);}
								bbq.append('\n');
							}
							if(bbq.length>=32768){
								myQOutstream.write(bbq.array, 0, bbq.length);
								bbq.setLength(0);
							}
						}
					}

//					if(bbq.length>0){
//						myQOutstream.write(bbq.array, 0, bbq.length);
//						bbq.setLength(0);
//					}
				}
//				assert(false) : OUTPUT_SAM+", "+SITES_ONLY+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT+"\n"+job.list.get(0).obj+"\n"+job.list.get(0);
				if(OUTPUT_SAM){
					assert(read1);
					for(final Read r : job.list){
						Read r2=(r==null ? null : r.mate);

						SamLine sl1=(r==null ? null : new SamLine(r, 0));
						SamLine sl2=(r2==null ? null : new SamLine(r2, 1));

						if(r!=null){
							
							if(verbose && r.numSites()>0){
								final Read clone=r.clone();
								for(SiteScore ss : r.sites){
									
									clone.setFromSite(ss);
									clone.setSecondary(true);
									SamLine sl=new SamLine(clone, 0);

									System.err.println("\n@************************************\n\n"+ss+"\n\n"+clone+"\n\n"+sl+"\n\n+************************************\n");
									
								}
							}
							
							assert(!ASSERT_CIGAR || !r.mapped() || sl1.cigar!=null) : r;
							sl1.toBytes(bb).append('\n');

							readsWritten++;
							basesWritten+=(r.bases!=null ? r.bases.length : 0);
							validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
							validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
							ArrayList<SiteScore> list=r.sites;
							if(OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
								final Read clone=r.clone();
								for(int i=1; i<list.size(); i++){
									SiteScore ss=list.get(i);
									clone.match=null;
									clone.setFromSite(ss);
									clone.setSecondary(true);
									
//									System.err.println(r.numericID+": "+(ss.match==null ? "null" : new String(ss.match)));
									
//									assert(false) : r.mapScore+"\n"+ss.header()+"\n"+r.sites+"\n";
									SamLine sl=new SamLine(clone, 0);
									assert(!sl.primary());
//									sl.setPrimary(false);
									

									assert(!ASSERT_CIGAR || sl.cigar!=null) : r;
									
									sl.toBytes(bb).append('\n');

//									readsWritten++;
//									basesWritten+=(r.bases!=null ? r.bases.length : 0);
//									validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
//									validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								}
							}
						}
						if(r2!=null){
							assert(!ASSERT_CIGAR || !r2.mapped() || sl2.cigar!=null) : r2;
							if(!SamLine.KEEP_NAMES && sl1!=null && ((sl2.qname==null) || !sl2.qname.equals(sl1.qname))){
								sl2.qname=sl1.qname;
							}
							sl2.toBytes(bb).append('\n');

							readsWritten++;
							basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
							validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
							validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
							
							ArrayList<SiteScore> list=r2.sites;
							if(OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
								final Read clone=r2.clone();
								for(int i=1; i<list.size(); i++){
									SiteScore ss=list.get(i);
									clone.match=null;
									clone.setFromSite(ss);
									clone.setSecondary(true);
//									assert(false) : r.mapScore+"\n"+ss.header()+"\n"+r.list+"\n";
									SamLine sl=new SamLine(clone, 0);
									assert(!sl.primary());
//									sl.setPrimary(false);
									
									assert(!ASSERT_CIGAR || sl.cigar!=null) : r2;
									
									sl.toBytes(bb).append('\n');

//									readsWritten++;
//									basesWritten+=(r.bases!=null ? r.bases.length : 0);
//									validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
//									validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								}
							}
						}
						if(bb.length>=32768){
							abd.write(bb.array, 0, bb.length);
							bb.setLength(0);
						}
						
					}
				}else if(SITES_ONLY){
					assert(read1);
					for(final Read r : job.list){
						Read r2=(r==null ? null : r.mate);
						
						if(r!=null && r.sites!=null){
							r.toSitesB(bb).append('\n');

							readsWritten++;
							basesWritten+=(r.bases!=null ? r.bases.length : 0);
							validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
							validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
						}
						if(r2!=null){
							r2.toSitesB(bb).append('\n');

							readsWritten++;
							basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
							validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
							validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
						}
						if(bb.length>=32768){
							abd.write(bb.array, 0, bb.length);
							bb.setLength(0);
						}
					}
				}else if(OUTPUT_FASTQ){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								r.toFastq(bb).append('\n');
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.bases.length : 0);
								validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
								validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									r2.toFastq(bb).append('\n');
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
									validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
								}
							}
							if(bb.length>=32768){
								abd.write(bb.array, 0, bb.length);
								bb.setLength(0);
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								r2.toFastq(bb).append('\n');
								readsWritten++;
								basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
								validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
								validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
							}
							if(bb.length>=32768){
								abd.write(bb.array, 0, bb.length);
								bb.setLength(0);
							}
						}
					}
				}else if(OUTPUT_FASTA){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								r.toFasta(bb).append('\n');
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.bases.length : 0);
								validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
								validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									r2.toFasta(bb).append('\n');
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
									validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
								}
							}
							if(bb.length>=32768){
								abd.write(bb.array, 0, bb.length);
								bb.setLength(0);
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								r2.toFasta(bb).append('\n');
								readsWritten++;
								basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
								validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
								validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
							}
							if(bb.length>=32768){
								abd.write(bb.array, 0, bb.length);
								bb.setLength(0);
							}
						}
					}
				}else if(OUTPUT_ATTACHMENT){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								if(r.obj==null){bb.append('.').append('\n');}
								else{bb.append(r.obj.toString()).append('.');}
								readsWritten++;
								validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									if(r2.obj==null){bb.append('.').append('\n');}
									else{bb.append(r2.obj.toString()).append('.');}
									readsWritten++;
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
								}
							}
							if(bb.length>=32768){
								abd.write(bb.array, 0, bb.length);
								bb.setLength(0);
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								if(r2!=null){
									if(r2.obj==null){bb.append('.').append('\n');}
									else{bb.append(r2.obj.toString()).append('.');}
									readsWritten++;
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
								}else{
									bb.append('.').append('\n');
								}
							}
							if(bb.length>=32768){
								abd.write(bb.array, 0, bb.length);
								bb.setLength(0);
							}
						}
					}
				}else{
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								r.toText(true, bb).append('\n');
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.bases.length : 0);
								validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
								validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									r2.toText(true, bb).append('\n');
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
									validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
								}
								
							}
							if(bb.length>=32768){
								abd.write(bb.array, 0, bb.length);
								bb.setLength(0);
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
//								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								if(r2!=null){
									r2.toText(true, bb).append('\n');
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
									validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
								}else{
									//TODO abd.print(".\n");
								}
							}
							if(bb.length>=32768){
								abd.write(bb.array, 0, bb.length);
								bb.setLength(0);
							}
						}
					}
				}
			}
			if(job.close){
				if(bb.length>0){
					abd.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
				assert(job.outstream!=null && job.outstream!=myOutstream);
				ReadWrite.finishWriting(null, job.outstream, fname, allowSubprocess); //TODO:  This should be job.fname
			}
			
			job=null;
			while(job==null){
				try {
					job=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		if(myOutstream!=null){
			if(bb.length>0){
				myOutstream.write(bb.array, 0, bb.length);
				bb.setLength(0);
			}
			ReadWrite.finishWriting(null, myOutstream, fname, allowSubprocess);
		}
		if(myQOutstream!=null){
			if(bbq.length>0){
				myQOutstream.write(bbq.array, 0, bbq.length);
				bbq.setLength(0);
			}
			ReadWrite.finishWriting(null, myQOutstream, qfname, allowSubprocess);
		}
		finishedSuccessfully=true;
	}

	private static final boolean buffered=true;
	private static final boolean ASSERT_CIGAR=false;
	private static final boolean verbose=false;
	
}
