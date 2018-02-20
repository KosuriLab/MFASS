package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;

public class ReadStreamStringWriter extends ReadStreamWriter {

	public ReadStreamStringWriter(String fname_, boolean read1_, int bufferSize, boolean allowSubprocess_){
		this(fname_, read1_, bufferSize, false, false, false, false, false, false, false, false, allowSubprocess_);
	}
	
//	public ReadStreamStringWriter(String fname_, boolean read1_, int bufferSize, 
//			boolean outputSamFile, boolean outputBamFile, boolean fastq, boolean fasta, boolean sitesOnly, boolean attachment, boolean stdout){
//		this(fname_, null, read1_, bufferSize, outputSamFile, outputBamFile, fastq, fasta, sitesOnly, attachment, stdout, false);
//	}

	public ReadStreamStringWriter(String fname_, boolean read1_, int bufferSize, 
			boolean outputSamFile, boolean outputBamFile, boolean fastq, boolean fasta, boolean sitesOnly, boolean attachment, boolean stdout,
			boolean useSharedHeader, boolean allowSubprocess_){
		super(fname_, null, read1_, bufferSize, 
				outputSamFile, outputBamFile, fastq, fasta, sitesOnly, attachment, stdout,
				useSharedHeader, true, true, allowSubprocess_);
	}

	public ReadStreamStringWriter(String fname_, String qfname_, boolean read1_, int bufferSize, 
			boolean outputSamFile, boolean outputBamFile, boolean fastq, boolean fasta, boolean sitesOnly, boolean attachment, boolean stdout,
			boolean useSharedHeader, boolean allowSubprocess_){
		super(fname_, qfname_, read1_, bufferSize, 
				outputSamFile, outputBamFile, fastq, fasta, sitesOnly, attachment, stdout,
				useSharedHeader, true, true, allowSubprocess_);
	}
	
	public ReadStreamStringWriter(FileFormat ff, String qfname_, boolean read1_, int bufferSize, CharSequence header, boolean useSharedHeader){
		super(ff, qfname_, read1_, bufferSize, header, true, true, useSharedHeader);
	}

	@Override
	public void run() {
		
		if(!OUTPUT_SAM && !OUTPUT_FASTQ && !OUTPUT_FASTA && !OUTPUT_ATTACHMENT){
			if(OUTPUT_INTERLEAVED){
//				assert(false) : OUTPUT_SAM+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT+", "+OUTPUT_INTERLEAVED+", "+SITES_ONLY;
				myWriter.print("#INTERLEAVED\n");
			}
			if(SITES_ONLY){
				myWriter.println("#"+SiteScore.header());
			}else if(!OUTPUT_ATTACHMENT){
				myWriter.println("#"+Read.header());
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
		
		while(job!=null && !job.poison){
//			System.err.println("Processing job "+job);
			if(!job.isEmpty()){
				
				if(myQWriter!=null){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								{
									CharSequence cs=(r.bases==null ? "\n" : toQualitySB(r.quality, r.bases.length).append('\n'));
									myQWriter.print('>');
									myQWriter.println(r.id);
									myQWriter.print(cs);
								}
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									CharSequence cs=(r2.bases==null ? "\n" : toQualitySB(r2.quality, r2.bases.length).append('\n'));
									myQWriter.print('>');
									myQWriter.println(r2.id);
									myQWriter.print(cs);
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								CharSequence cs=(r2.bases==null ? "\n" : toQualitySB(r2.quality, r2.bases.length).append('\n'));
								myQWriter.print('>');
								myQWriter.println(r2.id);
								myQWriter.print(cs);
							}
						}
					}
				}
//				assert(false) : OUTPUT_SAM+", "+SITES_ONLY+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT+"\n"+job.list.get(0).obj+"\n"+job.list.get(0);
				if(OUTPUT_SAM){
					assert(read1);
					for(final Read r : job.list){
						Read r2=(r==null ? null : r.mate);

						SamLine sl1=(r==null ? null : new SamLine(r, 0));
						SamLine sl2=(r2==null ? null : new SamLine(r2, 1));

						if(r!=null){
							job.writer.print(sl1.toText().append('\n'));

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
//									assert(false) : r.mapScore+"\n"+ss.header()+"\n"+r.list+"\n";
									SamLine sl=new SamLine(clone, 0);
									assert(!sl.primary());
//									sl.setPrimary(false);
									
									job.writer.print(sl.toText().append('\n'));

//									readsWritten++;
//									basesWritten+=(r.bases!=null ? r.bases.length : 0);
//									validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
//									validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								}
							}
						}
						if(r2!=null){
							if(!SamLine.KEEP_NAMES && sl1!=null && ((sl2.qname==null) || !sl2.qname.equals(sl1.qname))){
								sl2.qname=sl1.qname;
							}
							job.writer.print(sl2.toText().append('\n'));

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
									
									job.writer.print(sl.toText().append('\n'));

//									readsWritten++;
//									basesWritten+=(r.bases!=null ? r.bases.length : 0);
//									validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
//									validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								}
							}
						}
						
					}
				}else if(SITES_ONLY){
					assert(read1);
					for(final Read r : job.list){
						Read r2=(r==null ? null : r.mate);
						
						if(r!=null && r.sites!=null){
							job.writer.print(r.toSites().append('\n'));

							readsWritten++;
							basesWritten+=(r.bases!=null ? r.bases.length : 0);
							validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
							validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
						}
						if(r2!=null){
							job.writer.print(r2.toSites().append('\n'));

							readsWritten++;
							basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
							validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
							validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
						}
					}
				}else if(OUTPUT_FASTQ){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.print(r.toFastq().append('\n'));
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.bases.length : 0);
								validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
								validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.print(r2.toFastq().append('\n'));
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
									validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								job.writer.print(r2.toFastq().append('\n'));
								readsWritten++;
								basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
								validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
								validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
							}
						}
					}
				}else if(OUTPUT_FASTA){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.print(r.toFasta().append('\n'));
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.bases.length : 0);
								validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
								validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.print(r2.toFasta().append('\n'));
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
									validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								job.writer.print(r2.toFasta().append('\n'));
								readsWritten++;
								basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
								validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
								validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
							}
						}
					}
				}else if(OUTPUT_ATTACHMENT){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.println(r.obj==null ? "." : r.obj.toString());
								readsWritten++;
								validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.println(r2.obj==null ? "." : r2.obj.toString());
									readsWritten++;
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
//								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								if(r2!=null){
									job.writer.println(r2.obj==null ? "." : r2.obj.toString());
									readsWritten++;
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
								}else{
									job.writer.println(".");
								}
							}
						}
					}
				}else{
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.print(r.toText(true).append('\n'));
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.bases.length : 0);
								validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
								validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.bases.length : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.print(r2.toText(true).append('\n'));
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
									validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
								}
								
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
//								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								if(r2!=null){
									job.writer.print(r2.toText(true).append('\n'));
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.bases.length : 0);
									validReadsWritten+=(r2.valid() && r2.mapped() ? 1 : 0);
									validBasesWritten+=(r2.valid() && r2.mapped() && r2.bases!=null ? r2.bases.length : 0);
								}else{
									job.writer.print(".\n");
								}
							}
						}
					}
				}
			}
			if(job.close){
				assert(job.writer!=null && job.writer!=myWriter);
				ReadWrite.finishWriting(job.writer, job.outstream, fname, allowSubprocess);
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
		
		if(myWriter!=null){
			ReadWrite.finishWriting(myWriter, myOutstream, fname, allowSubprocess);
		}
		if(myQWriter!=null){
			ReadWrite.finishWriting(myQWriter, myQOutstream, qfname, allowSubprocess);
		}
		finishedSuccessfully=true;
	}
	
}
