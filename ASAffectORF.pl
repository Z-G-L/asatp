#!/usr/bin/env perl

=head2 NAME

Check AS event in CDS region.

=head2 SYNOPSIS

    Usage: perl ASAffectORF.pl --gtf <gtf file with CDS annotation>  --asEvent <AS event file> --output <output folder> 
    Options:
        -help|h
        --gtf          STRING           Input gtf format file with CDS annoation
        --output       STRING           Output folder.
        --asEvent      STRING           Output of program ASRecovist , i.e. "AS_event.xls"

=head2 AUTHOR

Zhigang Li	lzg0063(at)126.com	 2015-10-29

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
use lib::statistics;
use lib::parsing();
use lib::ASEvent();


#use Data::Dumper;

#parameters
my ( $gtfFile, $outputFolder, $asEventFile ) = @ARGV;

#help
if ( !@ARGV ) {
    pod2usage( -noperldoc => 1, -verbose => 2 );
    exit(0);
}
Getopt::Long::GetOptions(
    "help|h" => sub { pod2usage( -noperldoc => 1, -verbose => 2 ); exit(0); },
    'asEvent=s' => \$asEventFile,
    'gtf=s'     => \$gtfFile,
    'output=s'  => \$outputFolder,
);

################################################################################
#                              Main
################################################################################
&logMsg("Starting...");

# check output fold
&lib::basic::checkFold( $outputFolder, 1 );

########################################
&logMsg('Input gtf file...');

my %cdsInfoHash;

&lib::parsing::inputCDSInGtfFile( $gtfFile, \%cdsInfoHash );

########################################

&logMsg('Input AS event file...');
open ASE, "$asEventFile" || die;
my @asEvent = <ASE>;
close ASE;

########################################
&logMsg('Compare AS event and ORF info...');

my $cdsOutFile = $outputFolder . '/' . 'ASAffectORF_event.xls';
open CO, ">$cdsOutFile" || die;

$asEvent[0] =~ s/\s+$//;
my $title = $asEvent[0] . "\t" . "AS Event Location\tAS Event Frame Change\tTr2_vs_T1 ORF Diff Tag\n";
print CO $title;

# start column of the output, for summary
my @tmpCol = split /\t/, $asEvent[0];
my $outTableStartColumn = @tmpCol;    

for ( my $i = 1 ; $i < @asEvent ; $i++ ) {
    $asEvent[$i] =~ s/\s+$//;
    my @field = split /\t/, $asEvent[$i];
#next if($field[0] ne 'MGG_06174');
    # add column
    my $tag         = '';                                       # ORF change type
    my $location    = '';                                       # CDS UTR
    my $frameChange = '';                                       # 0 1 2 -

    # check locaion
    if ( exists $cdsInfoHash{ $field[0] } && exists $cdsInfoHash{ $field[0] }->{ $field[3] } ) {
        $location = &asLocation( $field[3], $field[7], $cdsInfoHash{ $field[0] }->{ $field[3] } );
    }
    if ( exists $cdsInfoHash{ $field[0] } && exists $cdsInfoHash{ $field[0] }->{ $field[4] } ) {
        my $tmp = &asLocation( $field[4], $field[7], $cdsInfoHash{ $field[0] }->{ $field[4] } );
        $location = ( $location eq '' ) ? ($tmp) : ( $location . ',' . $tmp );
    }
    $location = '-' if ( $location eq '' );                     # may be CDS is not annotated for some gene

    # check frame change
    $frameChange = &calFrameChange( $field[5], $field[8], $location );

    # check tag
    $tag = &checkTag( $cdsInfoHash{ $field[0] }->{ $field[3] }, $cdsInfoHash{ $field[0] }->{ $field[4] } );

    # output
    print CO $asEvent[$i], "\t", $location, "\t", $frameChange, "\t", $tag, "\n";
}

close CO;

########################################
&logMsg('Get summary...');

my $cdsOutSumFile = $outputFolder . '/' . 'ASAffectORF_event.summary.xls';
&eventSummary( $cdsOutFile, $outTableStartColumn, $cdsOutSumFile );

&logMsg("Finishing...");

################################################################################
#                           Subroutines
################################################################################

sub eventSummary {
    my ( $eventDiffFile, $outTableStartColumn, $eventSummary ) = @_;
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
                'tr'     => { "$field[0]::$field[3]" => 1, "$field[0]::$field[4]" => 1 },
                'as'     => 1,
                'type'   => $field[6],
                'in3Utr' => 0,
                'in5Utr' => 0,
                'inCds'  => 0,
                'frame0' => 0,
                'frame1' => 0,
                'frame2' => 0,
            };
        }

        # check location
        my $location = $field[$outTableStartColumn];
        foreach ( split /,/, $location ) {
            if ( $_ =~ /5UTR/ ) {
                $summary{ $field[5] }->{'in5Utr'}++;
            }
            elsif ( $_ =~ /3UTR/ ) {
                $summary{ $field[5] }->{'in3Utr'}++;
            }
            elsif ( $_ =~ /CDS/ ) {
                $summary{ $field[5] }->{'inCds'}++;
            }
        }

        # check frame
        my $frameChange = $field[ $outTableStartColumn + 1 ];
        if ( $frameChange eq '0' ) {
            $summary{ $field[5] }->{'frame0'}++;
        }
        elsif ( $frameChange eq '1' ) {
            $summary{ $field[5] }->{'frame1'}++;
        }
        elsif ( $frameChange eq '2' ) {
            $summary{ $field[5] }->{'frame2'}++;
        }

    }
    close EV;

    # output
    open SU, ">$eventSummary" || die;
    print SU "AS Pattern\tAS Event Type\tAS Event Num\tGene with AS Event\tTranscript with AS Event\t";
    print SU "AS Event In 5UTR\tAS Event In CDS\tAS Event In 3UTR\tFrameshift(0bp)\tFrameshift(1bp)\tFrameshift(2bp)\n";

    foreach my $asp ( keys %summary ) {
        print SU $asp, "\t";
        print SU $summary{$asp}->{'type'}, "\t";
        print SU $summary{$asp}->{'as'},   "\t";
        print SU scalar( keys %{ $summary{$asp}->{'gene'} } ), "\t";
        print SU scalar( keys %{ $summary{$asp}->{'tr'} } ),   "\t";
        print SU $summary{$asp}->{'in5Utr'}, "\t";
        print SU $summary{$asp}->{'inCds'},  "\t";
        print SU $summary{$asp}->{'in3Utr'}, "\t";
        print SU $summary{$asp}->{'frame0'}, "\t";
        print SU $summary{$asp}->{'frame1'}, "\t";
        print SU $summary{$asp}->{'frame2'}, "\n";
    }

    close SU;

    return 0;
}

sub checkTag {
    my ( $tr1CdsInfoHR, $tr2CdsInfoHR ) = @_;

    my $tag = '';

    my ( $start1, $end1, $strand1 ) = &lib::parsing::trStartEndStrand($tr1CdsInfoHR);
    my ( $start2, $end2, $strand2 ) = &lib::parsing::trStartEndStrand($tr2CdsInfoHR);

    my ( $os, $oe ) = &intersect( $start1, $end1, $start2, $end2 );

    # check whether ORF overlap with each other
    if ( $os < $oe && $oe - $os > 30 ) {

        # check frame
        my ( %cds1Frame, %cds2Frame );
        &cdsToFrame( $tr1CdsInfoHR, $start1, $end1, $strand1, \%cds1Frame );
        &cdsToFrame( $tr2CdsInfoHR, $start2, $end2, $strand2, \%cds2Frame );

        my ( $frameSameNum, $frameChangeNum ) = ( 0, 0 );

        for ( my $i = $os ; $i <= $oe ; $i++ ) {
            if ( exists $cds1Frame{$i} && exists $cds2Frame{$i} ) {
                #print STDERR "$i\t$cds1Frame{$i}\t$cds2Frame{$i}\n";
                ( $cds1Frame{$i} == $cds2Frame{$i} ) ? ( $frameSameNum++ ) : ( $frameChangeNum++ );
            }
        }
        if ( $frameSameNum == 0 && $frameChangeNum == 0 ) {

        }
        elsif ( $frameSameNum == 0 ) {
            $tag = ( $tag eq '' ) ? ("frame_full_change") : ( $tag . ',' . "frame_full_change" );
        }
        elsif ( $frameChangeNum == 0 ) {
            $tag = ( $tag eq '' ) ? ("frame_full_same") : ( $tag . ',' . "frame_full_same" );
        }
        else {
            $tag = ( $tag eq '' ) ? ("frame_part_same") : ( $tag . ',' . "frame_part_same" );
        }

        # check start and end
        my $startTag = '';
        my $stopTag  = '';
        if ( $strand1 ne '-' && $strand1 ne '-1' ) {
            if ( $start2 < $start1 ) {
                $startTag = ( $startTag eq '' ) ? ("start_gain") : ( $startTag . ',' . "start_gain" );
            }
            elsif ( $start2 > $start1 ) {
                $startTag = ( $startTag eq '' ) ? ("start_loss") : ( $startTag . ',' . "start_loss" );
            }
            else {
                $startTag = ( $startTag eq '' ) ? ("start_same") : ( $startTag . ',' . "start_same" );
            }

            if ( $end2 < $end1 ) {
                $stopTag = ( $stopTag eq '' ) ? ("stop_loss") : ( $stopTag . ',' . "stop_loss" );
            }
            elsif ( $end2 > $end1 ) {
                $stopTag = ( $stopTag eq '' ) ? ("stop_gain") : ( $stopTag . ',' . "stop_gain" );
            }
            else {
                $stopTag = ( $stopTag eq '' ) ? ("stop_same") : ( $stopTag . ',' . "stop_same" );
            }
        }
        else {
            if ( $end2 > $end1 ) {
                $startTag = ( $startTag eq '' ) ? ("start_gain") : ( $startTag . ',' . "start_gain" );
            }
            elsif ( $end2 < $end1 ) {
                $startTag = ( $startTag eq '' ) ? ("start_loss") : ( $startTag . ',' . "start_loss" );
            }
            else {
                $startTag = ( $startTag eq '' ) ? ("start_same") : ( $startTag . ',' . "start_same" );
            }

            if ( $start2 > $start1 ) {
                $stopTag = ( $stopTag eq '' ) ? ("stop_loss") : ( $stopTag . ',' . "stop_loss" );
            }
            elsif ( $start2 < $start1 ) {
                $stopTag = ( $stopTag eq '' ) ? ("stop_gain") : ( $stopTag . ',' . "stop_gain" );
            }
            else {
                $stopTag = ( $stopTag eq '' ) ? ("stop_same") : ( $stopTag . ',' . "stop_same" );
            }
        }

        $tag = ( $tag eq '' ) ? ($startTag) : ( $tag . ',' . $startTag );
        $tag = ( $tag eq '' ) ? ($stopTag)  : ( $tag . ',' . $stopTag );

    }
    else {
        $tag = 'CDS2 not span CDS1';
    }

    return $tag;
}

sub cdsToFrame {
    my ( $cdsInfoHR, $start, $end, $strand, $cdsPosiFrameHR ) = @_;
    my @cdsStartArray = keys %{$cdsInfoHR};

    if ( $strand ne '-' && $strand ne '-1' ) {    # strand +
    
        my $posi = 0;
        
        foreach my $cdsStart ( sort { $a <=> $b } @cdsStartArray ) {
            my $cdsEnd = $cdsInfoHR->{$cdsStart}->{end};
            for ( my $i = $cdsStart ; $i <= $cdsEnd ; $i++ ) {
                $posi++;
                $cdsPosiFrameHR->{$i} = $posi % 3;
            }
        }
    }
    else {                                        # strand -
    
        my $posi = 0;
        
        foreach my $cdsStart ( sort { $b <=> $a } @cdsStartArray ) {
            my $cdsEnd = $cdsInfoHR->{$cdsStart}->{end};

            for ( my $i = $cdsEnd ; $i >= $cdsStart ; $i-- ) {
                $posi++;
                $cdsPosiFrameHR->{$i} = $posi % 3;
            }
        }
    }

    return 0;
}

sub calFrameChange {
    my ( $asp, $binUnit, $location ) = @_;

    if ( $location !~ /\(CDS\),\S+\(CDS\)/ ) {
        return '-';
    }
    else {
        my ( $first, $second ) = &lib::ASEvent::asp2bin($asp);

        my $firstExonSumBP  = &exonSumBP( $first,  $binUnit );
        my $secondExonSumBP = &exonSumBP( $second, $binUnit );

        return abs( $secondExonSumBP - $firstExonSumBP ) % 3;
    }
}

sub exonSumBP {
    my ( $bit, $unit ) = @_;

    my $length = 0;

    my @bitArray  = split //,  $bit;
    my @unitArray = split /,/, $unit;

    if ( scalar(@bitArray) != scalar(@unitArray) ) {
        &logWarn("$bit/n$unit\n",__FILE__,__LINE__);
        &logError( "AS event bit length != span unit number!", __FILE__, __LINE__ );
    }

    for ( my $i = 0 ; $i < @bitArray ; $i++ ) {
        if ( $bitArray[$i] eq 1 ) {
            $unitArray[$i] =~ /(\d+)-(\d+)/;
            $length += ( $2 - $1 + 1 );
        }
    }
    return $length;
}

sub asLocation {
    my ( $trId, $spanUnit, $cdsInfoHR ) = @_;

    my $location = '';
    my ( $start, $end, $strand ) = &lib::parsing::trStartEndStrand($cdsInfoHR);

    $spanUnit =~ /(\d+)-(\d+)/;
    my $spanStart = $1;
    my $spanEnd   = $2;

    my ( $os, $oe ) = &intersect( $start, $end, $spanStart, $spanEnd );
    if ( $os <= $oe ) {
        $location = $trId . '(CDS)';
    }
    else {
        if ( $strand ne "-" && $strand ne '-1' ) {
            if ( $oe < $start ) {
                $location = $trId . '(5UTR)';
            }
            elsif ( $os > $end ) {
                $location = $trId . '(3UTR)';
            }
            else {
                &logError( "AS event region vs CDS region error", __FILE__, __LINE__ );
            }
        }
        else {
            if ( $oe < $start ) {
                $location = $trId . '(3UTR)';
            }
            elsif ( $os > $end ) {
                $location = $trId . '(5UTR)';
            }
            else {
                &logError( "AS event region vs CDS region error", __FILE__, __LINE__ );
            }
        }
    }
    return $location;

}

sub intersect {
    my ( $startA, $endA, $startB, $endB ) = @_;
    &logError( "Start > End in a range", __FILE__, __LINE__ ) if ( $startA > $endA || $startB > $endB );
    my ( $os, $oe );    #overlapStart, overlapEnd
    $os = ( $startA > $startB ) ? $startA : $startB;
    $oe = ( $endA > $endB )     ? $endB   : $endA;

    #check os oe
    return ( $os, $oe );
}

