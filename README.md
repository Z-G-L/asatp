

# Alternative Splicing Analysis Tool Package (ASATP)

Alternative splicing analysis tool package (ASATP), including a series of toolkits to analyze alternative splicing events, which could be used to 
 	detect and visualized alternative splicing events, 
 	check ORF changes, 
 	assess regulations of alternative splicing and 
 	do statistical analysis. 
The users could either access this software from our webserver, or download it from GitHub.

## Tools in ASATP
- ASRecovist

Alternative splicing recognition and visualization tool, which is a program to detect alternative splicing events from a gene annotation and classified them into different types (CE, IR, MXE, A3SS, A5SS, AFE, ALE and other). Alternative splicing events will be showed in tables and graphs. 

- ASQuantityDiff

Alternative splicing quantity comparison between samples, which was used to detect expression regulations of alternative transcripts.

- ASAffectORF

Check AS event in CDS region, to identify the effect of alternative splicing event for ORFs.

- runASATP

Run Alternative splicing Analysis Tool Package, which was a combined pipeline including ASRecovist, ASQuantityDiff and ASAffectORF.

- asp2bit

Transform AS pattern to bit matrix.

- bit2asp

Transform bit to ASP code.

- splitGtf

Split a GTF file when it's too large to process. Then you can process them separately using multi-CPU to save time.


## How to install?

ASATP can be downloaded from https://github.com/Z-G-L/ASATP

To use ASATP, you need install Perl(>5.10). And some Perl modules are needed:

- Bioperl
- Bio::Graphics
- GD::Image
- GD::SVG
- Statistics::R
- Math::BigInt

You can use CPAN command to install modules needed. For example: 

perl -MCPAN -e 'install HTML::Template' 

## How to use?

Detail see Full_documentation.pdf in document fold.





