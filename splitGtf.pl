#!/usr/bin/env perl

=head2 NAME

Split a gtf file when it's too large to process. 

=head2 SYNOPSIS

perl splitGtf.pl <file.gtf> <output_prefix>

=head2 AUTHOR

Zhigang Li	lzg0063(at)126.com	 2014-11-23

=cut

################################################################################
#                             Options
################################################################################
BEGIN { use FindBin qw($Bin); use lib "$Bin"; }
use 5.010;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use lib::basic qw(logMsg logError logWarn logStd);

#use Data::Dumper;

#parameters

my $geneNumPerFile = 10000;

#help
if ( !@ARGV ) {
    pod2usage( -noperldoc => 1, -verbose => 2 );
    exit(0);
}
Getopt::Long::GetOptions( "help|h" => sub { pod2usage( -noperldoc => 1, -verbose => 2 ); exit(0); }, );

################################################################################
#                              Main
################################################################################
&logMsg("Starting...");

my %fileHash;

open FI, "$ARGV[0]" || die;
while (<FI>) {
    next if ( $_ =~ /^#/ );
    next if ( $_ !~ /\S+/ );
    $_ =~ s/\R//;
    my @field = split /\t/, $_;
    if ( $field[2] eq 'exon' || $field[2] eq 'EXON' || $field[2] eq 'CDS' || $field[2] eq 'cds' ) {
        $field[8] =~ /gene_id "([^"]+)"/;
        my $id = $1;

        if ( exists $fileHash{$id} ) {
            $fileHash{$id} .= $_ . "\n";
        }
        else {
            $fileHash{$id} = $_ . "\n";
        }
    }
}
close FI;

my $num    = 0;
my $fileId = 1;

my $oFile = $ARGV[1] .'.'.$fileId.".gtf";
open OO, ">$oFile" || die;

foreach ( keys %fileHash ) {
    $num++;

    if ( $num > $geneNumPerFile ) {
        close OO;

        $fileId++;
        $num = 0;

        my $oFile = $ARGV[1] .'.'.$fileId.".gtf";
        open OO, ">$oFile" || die;

    }
    
    print OO $fileHash{$_}; 
}

close OO;

&logMsg("Finishing...");

################################################################################
#                           Subroutines
################################################################################


