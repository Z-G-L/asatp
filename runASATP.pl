#!/usr/bin/env perl

=head2 NAME

Run Alternative splicing Analysis Tool Package.

=head2 SYNOPSIS

    Usage: perl runASATP.pl --gtf <gtf file> --trExpFile <transcript_expression_file> --output <output_folder> [--graph no] [--graphFormat SVG]
    Options:
        -help|h
        --gtf          STRING           Input gtf format file.
        --output       STRING           Output folder.
        --trExpFile    STRING           A file with expression levels of transcripts in different samples.
        --graph        no|yes           Create graph or not. default [no]
        --graphFormat  SVG|png          Graph format. default [SVG]

    Note:
    
        "--gtf" input file should contain CDS annotation.
    
        "--trExpFile" input file format (column separated by Tab):
        Gene    Transcript  Sample1 Sample2 ...
        g1  tr1 0.5 20  ...
        g1  tr2 53 19  ...

=head2 AUTHOR

Zhigang Li	lzg0063(at)126.com	 2014-10-30

=cut

################################################################################
#                             Options
################################################################################
BEGIN { use FindBin qw($Bin); use lib "$Bin"; }
use 5.010;
use strict;
#use warnings;
use Getopt::Long;
use Pod::Usage;
use lib::basic qw(logMsg logError logWarn logStd);

#use Data::Dumper;

#parameters
my ( $gtfFile, $expFile, $outputFolder ) = @ARGV;
my $outGraph    = 'no';    #whether output graph
my $graphFormat = 'svg';

#help
if ( !@ARGV ) {
    pod2usage( -noperldoc => 1, -verbose => 2 );
    exit(0);
}
Getopt::Long::GetOptions(
    "help|h" => sub { pod2usage( -noperldoc => 1, -verbose => 2 ); exit(0); },
    "gtf=s" => \$gtfFile,
    'graph=s'       => \$outGraph,
    'graphFormat=s' => \$graphFormat,
    'trExpFile=s'   => \$expFile,
    'output=s'      => \$outputFolder,
);

################################################################################
#                              Main
################################################################################
&logMsg("Starting...");

# check output fold
&lib::basic::checkFold( $outputFolder, 1 );




&logMsg("Running ASRecovist...");

my $asrecovistOut = $outputFolder . '/' . 'ASRecovist_out';


&logStd("perl $Bin/ASRecovist.pl --gtf $gtfFile --output $asrecovistOut --graph $outGraph --graphFormat $graphFormat");
my @arg = system("perl $Bin/ASRecovist.pl --gtf $gtfFile --output $asrecovistOut --graph $outGraph --graphFormat $graphFormat");

exit(1) if($? != 0);


########################################


&logMsg("Running ASAffectORF...");

my $asaffectOut = $outputFolder . '/' . 'ASAffectORF_out';

&logStd("perl $Bin/ASAffectORF.pl --gtf $gtfFile --asEvent $asrecovistOut/AS_event.xls --output $asaffectOut");
system("perl $Bin/ASAffectORF.pl --gtf $gtfFile --asEvent $asrecovistOut/AS_event.xls --output $asaffectOut");



exit(1) if($? != 0);



########################################



&logMsg("Running ASQuantityDiff...");

my $asquanOut = $outputFolder . '/' . 'ASQuantityDiff_out';

&logStd("perl $Bin/ASQuantityDiff.pl --asEvent $asaffectOut/ASAffectORF_event.xls --asGroup $asrecovistOut/AS_event_group.xls --trExpFile $expFile --output $asquanOut");
system(
    "perl $Bin/ASQuantityDiff.pl --asEvent $asaffectOut/ASAffectORF_event.xls --asGroup $asrecovistOut/AS_event_group.xls --trExpFile $expFile --output $asquanOut"
);


exit(1) if($? != 0);





&logMsg("Finishing...");

################################################################################
#                           Subroutines
################################################################################
