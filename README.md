

# Alternative Splicing Analysis Tool Package (ASATP)

Alternative splicing analysis tool package (ASATP) includes a series of toolkits to analyze alternative splicing events, which could be used to 
 	detect and visualized alternative splicing events, 
 	check ORF changes, 
 	assess regulations of alternative splicing and 
 	do statistical analysis. 

## Tools in ASATP
- ASRecovist

Alternative splicing recognition and visualization tool, which is a program to detect alternative splicing events from a gene annotation and classified them into different types (CE, IR, MXE, A3SS, A5SS, AFE, ALE and other). Alternative splicing events will be showed in tables and graphs. 

- ASQuantityDiff

Alternative splicing quantity comparison between samples, which was used to detect expression regulations of alternative transcripts.

- ASAffectORF

Check AS event in CDS region, to identify the effect of alternative splicing events for ORFs.

- runASATP

Run Alternative splicing Analysis Tool Package, which was a combined pipeline including ASRecovist, ASQuantityDiff and ASAffectORF.

- asp2bit

Transform AS pattern to bit matrix.

- bit2asp

Transform bit to ASP code.

- splitGtf

Split a GTF file when it's too large to process. Then you can process them separately using multi-CPU to save time.


## How to install?

ASATP works under Linux. To use ASATP, you need install Perl (>5.10). Some Perl modules are needed:

- Bioperl
- Bio::Graphics
- GD::Image
- GD::SVG
- Statistics::R
- Math::BigInt

You can use CPAN command to install modules needed. For example: 


	perl -MCPAN -e 'install HTML::Template' 

	
To check whether ASATP work well, you can run commands like:

	cd Sample	

	perl ../runASATP.pl --gtf anno.gtf  --trExpFile isoform_exp.txt  --output output_svg --graph yes --graphFormat svg

	perl ../runASATP.pl --gtf anno.gtf  --trExpFile isoform_exp.txt  --output output_png --graph yes --graphFormat png


## How to use?

See "Full_documentation.pdf" in document folder.





