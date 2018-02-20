ReformatReads readme by Brian Bushnell
Last updated February 18, 2014.
Please contact me at bbushnell@lbl.gov if you have any questions or encounter any errors.

This is currently a stub.

Change Log:

Moved processing phase from constructor into a new function.
Added pigz (parallel gzip) support, at suggestion of Rob Egan.
Added ftr/ftl (force trim to a certain base position), at request of Alicia Clum.
Added shellscript support for java -Xmx flag (Suggested by James Han).
Fixed stdout.fa.gz writing uncompressed via ReadStreamWriter.
Added scarf input support.
Major refactoring.
Output quality switch (qout) was being ignored and ASCII33 was always used; this has been fixed.
TrimRead.testOptimal() mode added, and made default when quality trimming is performed; old mode can be used with 'otf=f' flag.
Fixed bug in ByteBuilder for reads with no quality values.  Noted by Brian Foster.
All program message information now defaults to stderr.