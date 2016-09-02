#!/usr/bin/env perl

=head2 NAME

Alternative splicing recognition and visualization tool . 

=head2 SYNOPSIS

    Usage: perl ASRecovist.pl --gtf <gtf_fortmat_file> --output <output_fold> [--graph no] [--graphFormat SVG]
    Options:
        -help|h
        --gtf          STRING           Input gtf format file.
        --output       STRING           Output folder.
        --graph        no|yes           Create graph or not. default [no]
        --graphFormat  SVG|png          Graph format. default [SVG]

=head2 AUTHOR

Zhigang Li	lzg0063(at)126.com	 2014-3-4

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
#use Data::Dumper;
use lib::basic qw(logMsg logError logWarn logStd);
use lib::parsing();
use lib::ASEvent();
use lib::graphic();

#parameters
my ( $gtfFile, $outputFolder );
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
    'output=s'      => \$outputFolder,
    'graph=s'       => \$outGraph,
    'graphFormat=s' => \$graphFormat,
);

################################################################################
#                              Main
################################################################################
&logMsg('Starting...');

# check output fold
&lib::basic::checkFold( $outputFolder, 1 );
my $outputGraphPath;
if ( $outGraph eq 'yes' ) {
    $outputGraphPath = $outputFolder . '/' . 'ASRecovist_graph';
    &lib::basic::checkFold( $outputGraphPath, 1 );
}

########################################

# AS event info file
my $outDetailFile = $outputFolder . '/' . 'AS_event.xls';
open ODF, ">$outDetailFile" || die;
print ODF "Gene\tChromosome\tStrand\tTranscript1\tTranscript2\t";
print ODF "AS Pattern\tAS Event Type\tAS Event Span Unit\tAS Bit Span Unit\n";

# AS event group file
my $outGeneASGroup = $outputFolder . '/' . 'AS_event_group.xls';
open GEC, ">$outGeneASGroup" || die;
print GEC "Gene\tChromosome\tStrand\tAS Event Groups\tAS Pattern\tAS Event Type\t";
print GEC "AS Event Span Unit\tAS Bit Span Unit\tTranscript1\tTranscript2\n";

########################################

&logMsg('Input gtf file...');

my %gtfFileHash;

&lib::parsing::inputGtfFile( $gtfFile, \%gtfFileHash );

########################################

&logMsg('Detect AS event and output...');

# each gene
foreach my $geneId ( keys %gtfFileHash ) {
    my %geneASGroupInfo;    # AS groups in a gene

    # each pair of tr
    my $geneChr    = &lib::parsing::geneChr( $gtfFileHash{$geneId} );      # chr or scaffold of this gene
    my $geneStrand = &lib::parsing::geneStrand( $gtfFileHash{$geneId} );

    my @allTrId = sort { $a cmp $b } keys %{ $gtfFileHash{$geneId} };
    for ( my $i = 0 ; $i < @allTrId ; $i++ ) {                             # transcript 1
        for ( my $j = $i + 1 ; $j < @allTrId ; $j++ ) {                    # transcript 2
            my $tr1Id = $allTrId[$i];
            my $tr2Id = $allTrId[$j];

            # store AS event information
            # $asEvent{$asNumOfTr}=[$asp, $asType, $eventTr1Bit,$eventTr2Bit,$eventBitRegion,$geneStrand, $asSpanUnit, $aspScoreReverseTrId]
            # $aspScoreReverseTrId: if 0, bit string corresponding to ASP is tr1,tr2; if 1 is tr2,tr1
            my %asEvent;
            &lib::ASEvent::detectASEvent( $gtfFileHash{$geneId}, $tr1Id, $tr2Id, \%asEvent );

            # same AS event in a gene will be classified into the same group
            &lib::ASEvent::ASEventGroup( $tr1Id, $tr2Id, \%asEvent, \%geneASGroupInfo );

            #print Dumper( \%asEvent );

            # output AS event
            &outputASEvent( \*ODF, $geneId, $geneChr, $geneStrand, $tr1Id, $tr2Id, \%asEvent );
        }
    }

    # output AS group information
    &outputASGroup( \*GEC, $geneId, $geneChr, \%geneASGroupInfo );

    #print Dumper( \%geneASGroupInfo );

    ############################################################
    # create graph
    if ( $outGraph eq 'yes' ) {
        my $asGraphFile = $outputGraphPath . '/' . $geneId;
        my $groupSum    = scalar( keys %geneASGroupInfo );
        if ( $groupSum > 0 ) {
            &lib::graphic::createGraph( $geneId, $gtfFileHash{$geneId}, \%geneASGroupInfo, $graphFormat, $asGraphFile );
        }
    }

}
close ODF;
close GEC;

########################################
&logMsg("Get AS_event summary...");

my $eventSummary = $outputFolder . '/' . 'AS_event.summary.xls';

&eventSummary( $outDetailFile, $eventSummary );

########################################

&logMsg("Get AS_event_group summary...");

my $groupSummary = $outputFolder . '/' . 'AS_event_group.summary.xls';

&groupSummary( $outGeneASGroup, $groupSummary );

&logMsg('Finishing...');

################################################################################
#                           Subroutines
################################################################################

sub groupSummary {
    my ( $groupDiffFile, $groupSummary ) = @_;
    open EV, "$groupDiffFile" || die;

    my %summary;
    while (<EV>) {
        $_ =~ s/\s+$//;
        next if ( $_ eq '' );
        next if ( $_ =~ /^Gene\t/ );
        my @field = split /\t/, $_;

        if ( exists $summary{ $field[4] } ) {
            $summary{ $field[4] }->{'as'}++;
            $summary{ $field[4] }->{'gene'}->{ $field[0] } = 1;
            map( $summary{ $field[4] }->{'tr'}->{"$field[0]::$_"} = 1, split /,/, $field[8] );
            map( $summary{ $field[4] }->{'tr'}->{"$field[0]::$_"} = 1, split /,/, $field[9] );
        }
        else {
            $summary{ $field[4] } = {
                'gene' => { $field[0] => 1 },
                'tr'   => {},
                'as'   => 1,
                'type' => $field[5],
            };
            map( $summary{ $field[4] }->{'tr'}->{"$field[0]::$_"} = 1, split /,/, $field[8] );
            map( $summary{ $field[4] }->{'tr'}->{"$field[0]::$_"} = 1, split /,/, $field[9] );
        }
    }
    close EV;

    # output
    open SU, ">$groupSummary" || die;
    print SU "AS Pattern\tAS Event Type\tAS Event Group Num\tGene with AS Event\tTranscript with AS Event\n";

    foreach my $asp ( keys %summary ) {
        print SU $asp, "\t";
        print SU $summary{$asp}->{'type'}, "\t";
        print SU $summary{$asp}->{'as'},   "\t";
        print SU scalar( keys %{ $summary{$asp}->{'gene'} } ), "\t";
        print SU scalar( keys %{ $summary{$asp}->{'tr'} } ),   "\n";
    }

    close SU;

    return 0;
}

sub eventSummary {
    my ( $eventDiffFile, $eventSummary ) = @_;
    open EV, "$eventDiffFile" || die;

    my %summary;
    while (<EV>) {
        $_ =~ s/\s+$//;
        next if ( $_ eq '' );
        next if ( $_ =~ /^Gene\t/ );
        my @field = split /\t/, $_;

        if ( exists $summary{ $field[5] } ) {
            $summary{ $field[5] }->{'gene'}->{ $field[0] }          = 1;
            $summary{ $field[5] }->{'tr'}->{"$field[0]::$field[3]"} = 1;
            $summary{ $field[5] }->{'tr'}->{"$field[0]::$field[4]"} = 1;
            $summary{ $field[5] }->{'as'}++;
        }
        else {
            $summary{ $field[5] } = {
                'gene' => { $field[0] => 1 },
                'tr'   => { "$field[0]::$field[3]" => 1, "$field[0]::$field[4]" => 1 },
                'as'   => 1,
                'type' => $field[6],
            };
        }

    }
    close EV;

    # output
    open SU, ">$eventSummary" || die;
    print SU "AS Pattern\tAS Event Type\tAS Event Num\tGene with AS Event\tTranscript with AS Event\n";

    foreach my $asp ( keys %summary ) {
        print SU $asp, "\t";
        print SU $summary{$asp}->{'type'}, "\t";
        print SU $summary{$asp}->{'as'},   "\t";
        print SU scalar( keys %{ $summary{$asp}->{'gene'} } ), "\t";
        print SU scalar( keys %{ $summary{$asp}->{'tr'} } ),   "\n";
    }

    close SU;

    return 0;
}

# Title   : outputASEvent
sub outputASEvent {
    my ( $fileHandle, $geneId, $geneChr, $geneStrand, $tr1Id, $tr2Id, $asEventPairTrHR ) = @_;

    my @eventSum = keys %{$asEventPairTrHR};
    if ( @eventSum > 0 ) {
        foreach my $num (@eventSum) {
            print $fileHandle $geneId,     "\t";
            print $fileHandle $geneChr,    "\t";
            print $fileHandle $geneStrand, "\t";
            print $fileHandle $tr1Id,      "\t", $tr2Id, "\t";
            print $fileHandle $asEventPairTrHR->{$num}->[0], "\t";
            print $fileHandle $asEventPairTrHR->{$num}->[1], "\t";
            print $fileHandle $asEventPairTrHR->{$num}->[6], "\t";
            print $fileHandle $asEventPairTrHR->{$num}->[4], "\n";
        }
    }

    #else {
    #    print $fileHandle $geneId,     "\t";
    #    print $fileHandle $geneChr,    "\t";
    #    print $fileHandle $geneStrand, "\t";
    #    print $fileHandle $tr1Id,      "\t", $tr2Id, "\t";
    #    print $fileHandle 'None', "\t";
    #    print $fileHandle 'None', "\t";
    #    print $fileHandle 'None', "\t";
    #    print $fileHandle 'None', "\n";
    #}
    return 0;
}

sub outputASGroup {
    my ( $fileHandle, $geneId, $geneChr, $asGroupHR ) = @_;
    foreach my $identity ( keys %{$asGroupHR} ) {
        print GEC $geneId,  "\t";
        print GEC $geneChr, "\t";
        print GEC $asGroupHR->{$identity}->{strand},   "\t";
        print GEC $asGroupHR->{$identity}->{groupId},  "\t";
        print GEC $asGroupHR->{$identity}->{asp},      "\t";
        print GEC $asGroupHR->{$identity}->{astype},   "\t";
        print GEC $asGroupHR->{$identity}->{spanUnit}, "\t";
        print GEC $asGroupHR->{$identity}->{asRegion}, "\t";
        print GEC $asGroupHR->{$identity}->{tr1Ids},   "\t";
        print GEC $asGroupHR->{$identity}->{tr2Ids},   "\n";

    }
    return 0;
}
