#!/usr/bin/env perl

=head2 NAME

Transform AS pattern to bit matrix.

=head2 SYNOPSIS

perl asp2bit.pl --asp aspCode --output <output_fold> [--graph no] [--graphFormat SVG]

    Options:
        -help|h
        --asp                           Alternative splicing pattern code, e.g. 16
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
my ($asp,$outputFolder) = @ARGV;

my $outGraph    = 'no';    #whether output graph
my $graphFormat = 'svg';

#help
if ( !@ARGV ) {
        pod2usage( -noperldoc => 1, -verbose => 2 );
        exit(0);
}
Getopt::Long::GetOptions(
        "help|h" => sub { pod2usage( -noperldoc => 1, -verbose => 2 ); exit(0); },
        "asp=s" => \$asp,
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
    $outputGraphPath = $outputFolder . '/' . 'asp2bin_graph';
    &lib::basic::checkFold( $outputGraphPath, 1 );
}


my $oFile = $outputFolder.'/'.'text.txt';
open TX,">$oFile"||die;

my ( $first, $second ) = &lib::ASEvent::asp2bin($asp);

my $bit = $first . ',' . $second;


my $check = &lib::ASEvent::bitCheck($bit);
if($check == 1){
        &logStd("\nThis ASP code is not an AS event. Please check it again!\n");
        print TX "\nThis ASP code is not an AS event. Please check it again!\n";
        exit;
}

my $asType = &lib::ASEvent::ASTypeCheck($asp);



print TX "AS Pattern : $asp\n";
print TX "AS Type    : $asType\n";
print TX "Bit        : $bit\n";
close TX;


########################################

# create graph
if ( $outGraph eq 'yes' ) {
        my $asGraphFile = $outputGraphPath.'/asp_bit_graph';
        my %artificialGeneInfo;
        &lib::parsing::bit2geneInfo( $bit, \%artificialGeneInfo );
        my $geneId = "ASP: $asp   AS Type: $asType   Bit: $bit";
        &lib::graphic::bitGraph( $geneId, \%artificialGeneInfo, $graphFormat, $asGraphFile );
}

