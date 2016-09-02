#!/usr/bin/env perl

=head2 NAME

Transform bit to ASP code.

=head2 SYNOPSIS

perl bit2asp.pl --bit bitCode --output <output_fold> [--graph no] [--graphFormat SVG]

    Options:
        -help|h
        --bit                           Bit code, e.g. 10001,--101
        --graph        no|yes           Create graph or not. default [no]
        --graphFormat  SVG|png          Graph format. default [SVG]
        --output       STRING           Output folder.

=head2 AUTHOR

Zhigang Li	lzg0063(at)126.com	 2015-10-16

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
use lib::graphic();
use lib::ASEvent();

#use Data::Dumper;

#parameters
my ($bitPair,$outputFolder);

my $outGraph    = 'no';    #whether output graph
my $graphFormat = 'svg';


#help
if ( !@ARGV ) {
        pod2usage( -noperldoc => 1, -verbose => 2 );
        exit(0);
}
Getopt::Long::GetOptions(
        "help|h" => sub { pod2usage( -noperldoc => 1, -verbose => 2 ); exit(0); },
        "bit=s" => \$bitPair,
        'output=s'      => \$outputFolder,
        'graph=s'       => \$outGraph,
        'graphFormat=s' => \$graphFormat,
);


################################################################################
#                              Main
################################################################################

# check output fold
&lib::basic::checkFold( $outputFolder, 1 );
my $outputGraphPath;
if ( $outGraph eq 'yes' ) {
    $outputGraphPath = $outputFolder . '/' . 'bin2asp_graph';
    &lib::basic::checkFold( $outputGraphPath, 1 );
}

my $oFile = $outputFolder.'/'.'text.txt';
open TX,">$oFile"||die;

my $check = &lib::ASEvent::bitCheck($bitPair);
if($check == 1){
        &logStd("\nThis ASP code is not an AS event. Please check it again!\n");
        print TX "\nThis ASP code is not an AS event. Please check it again!\n";
        exit;
}

my @part = split /\s*,\s*/, $bitPair;

my @tr1Bit = split //, $part[0];
my @tr2Bit = split //, $part[1];

&logError( "Error bit format!  Please try again!\n", __FILE__, __LINE__ ) if ( @part != 2 || @tr1Bit != @tr2Bit);

my @asEvent = (0, @tr1Bit-1);

my ($asp)=&lib::ASEvent::ASPatternScore(\@tr1Bit,\@tr2Bit,\@asEvent);

my $asType = &lib::ASEvent::ASTypeCheck($asp);




print TX "Bit        : $bitPair\n";
print TX "AS Pattern : $asp\n";
print TX "AS Type    : $asType\n";
close TX;

########################################

# create graph
if ( $outGraph eq 'yes' ) {
        my $asGraphFile = $outputGraphPath.'/asp_bit_graph';
        my %artificialGeneInfo;
        &lib::parsing::bit2geneInfo( $bitPair, \%artificialGeneInfo );
        my $geneId = "ASP: $asp   AS Type: $asType   Bit: $bitPair";
        &lib::graphic::bitGraph( $geneId, \%artificialGeneInfo, $graphFormat, $asGraphFile );
}










