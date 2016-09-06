#!/usr/bin/env perl

=head2 NAME

Alternative splicing quantity diff comparison betwwen samples.

=head2 SYNOPSIS

    Usage: perl ASQuantityDiff.pl --asEvent <AS_event.xls> --asGroup <AS_event_group.xls> --trExpFile <transcript_expression_file> --output <output_fold>
    Options:
        -help|h
        --output       STRING           Output folder.
        --asGroup      STRING           Output of program ASRecovist, i.e. "AS_event_group.xls"
        --asEvent      STRING           Output of program ASRecovist or ASAffectORF, i.e. "AS_event.xls" or ASAffectORF_event.xls
        --trExpFile    STRING           A file with expression levels of transcripts in different samples.
        --qvalue       FLOAT            q-vlaue cutoff [default: 0.05]
        --expCutoff    FLOAT            Expression level cutoff. 
                                        A transcript will be considered to be not expressed if its expression level less than this cutoff.
        
    Note:
        "--trExpFile" input file format (column separated by Tab):
        Gene    Transcript  Sample1 Sample2 ...
        g1  tr1 0.5 20  ...
        g1  tr2 53 19  ...

=head2 AUTHOR

Zhigang Li  lzg0063(at)126.com   2013-10-22

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

use Data::Dumper;

#parameters
my ( $asGroupFile, $asEventFile, $expFile, $outputFolder );
my $qvalueCutoff = 0.05;
my $expCutoff    = 0;

#help
if ( !@ARGV ) {
    pod2usage( -noperldoc => 1, -verbose => 2 );
    exit(0);
}
Getopt::Long::GetOptions(
    "help|h" => sub { pod2usage( -noperldoc => 1, -verbose => 2 ); exit(0); },
    'asGroup=s'   => \$asGroupFile,
    'asEvent=s'   => \$asEventFile,
    'trExpFile=s' => \$expFile,
    'output=s'    => \$outputFolder,
    'qvalue=f'    => \$qvalueCutoff,
    'expCutoff=f' => \$expCutoff,

);

################################################################################
#                              Main
################################################################################
&logMsg("Starting...");

# check output fold
&lib::basic::checkFold( $outputFolder, 1 );

########################################

&logMsg("Input transcript expresion file...");

my %trExp;
my @samples;

&inputExp( $expFile, \%trExp, \@samples );

########################################

&logMsg("Process AS_event...");

my $eventDiffFile = $outputFolder . '/' . 'ASExpDiff_event.xls';

my $outTableStartColumnEvent = &processEvent( \%trExp, \@samples, $asEventFile, $eventDiffFile );

########################################

&logMsg("Process AS_event_event_group...");

my $groupDiffFile = $outputFolder . '/' . 'ASExpDiff_event_group.xls';

my $outTableStartColumnGroup = &processGroup( \%trExp, \@samples, $asGroupFile, $groupDiffFile );

########################################

&logMsg("Get ASExpDiff_event summary...");

my $eventSummary = $outputFolder . '/' . 'ASExpDiff_event.summary.xls';

&eventSummary( $eventDiffFile, $eventSummary, $qvalueCutoff, $expCutoff, $outTableStartColumnEvent );

########################################

&logMsg("Get ASExpDiff_event_group summary...");

my $groupSummary = $outputFolder . '/' . 'ASExpDiff_event_group.summary.xls';

&groupSummary( $groupDiffFile, $groupSummary, $qvalueCutoff, $expCutoff, $outTableStartColumnGroup );

&logMsg("Finishing...");

################################################################################
#                           Subroutines
################################################################################

sub groupSummary {
    my ( $groupDiffFile, $groupSummary, $qvalueCutoff, $expCutoff, $outTableStartColumn ) = @_;
    open EV, "$groupDiffFile" || die;

    my %summary;
    while (<EV>) {
        $_ =~ s/\s+$//;
        next if ( $_ eq '' );
        next if ( $_ =~ /^Gene\t/ );
        my @field = split /\t/, $_;

        my $pair = $field[$outTableStartColumn] . "::" . $field[ $outTableStartColumn + 1 ];
        if ( exists $summary{$pair}->{ $field[4] } ) {
            $summary{$pair}->{ $field[4] }->{'as'}++;
            $summary{$pair}->{ $field[4] }->{'gene'}->{ $field[0] } = 1;
            map( $summary{$pair}->{ $field[4] }->{'tr'}->{"$field[0]::$_"} = 1, split /,/, $field[8] );
            map( $summary{$pair}->{ $field[4] }->{'tr'}->{"$field[0]::$_"} = 1, split /,/, $field[9] );

        }
        else {
            $summary{$pair}->{ $field[4] } = {
                'gene'       => { $field[0] => 1 },
                'tr'         => {},
                'as'         => 1,
                'type'       => $field[5],
                'asdiff'     => 0,
                'genediff'   => {},
                'trdiff'     => {},
                's1specific' => 0,
                's2specific' => 0,
            };
            map( $summary{$pair}->{ $field[4] }->{'tr'}->{"$field[0]::$_"} = 1, split /,/, $field[8] );
            map( $summary{$pair}->{ $field[4] }->{'tr'}->{"$field[0]::$_"} = 1, split /,/, $field[9] );
        }

        # check exp diff
        if (   $field[ $outTableStartColumn + 7 ] ne '-'
            && $field[ $outTableStartColumn + 7 ] <= $qvalueCutoff
            && $field[ $outTableStartColumn + 2 ] >= $expCutoff
            && $field[ $outTableStartColumn + 3 ] >= $expCutoff
            && $field[ $outTableStartColumn + 4 ] >= $expCutoff
            && $field[ $outTableStartColumn + 5 ] >= $expCutoff )
        {
            $summary{$pair}->{ $field[4] }->{'asdiff'}++;
            $summary{$pair}->{ $field[4] }->{'genediff'}->{ $field[0] } = 1;
            map( $summary{$pair}->{ $field[4] }->{'trdiff'}->{"$field[0]::$_"} = 1, split /,/, $field[8] );
            map( $summary{$pair}->{ $field[4] }->{'trdiff'}->{"$field[0]::$_"} = 1, split /,/, $field[9] );
        }

        # check sample specific as event
        $summary{$pair}->{ $field[4] }->{'s1specific'}++
          if ( ( $field[ $outTableStartColumn + 2 ] > $expCutoff && $field[ $outTableStartColumn + 3 ] > $expCutoff )
            && ( $field[ $outTableStartColumn + 4 ] <= $expCutoff || $field[ $outTableStartColumn + 5 ] <= $expCutoff ) );
        $summary{$pair}->{ $field[4] }->{'s2specific'}++
          if ( ( $field[ $outTableStartColumn + 2 ] <= $expCutoff || $field[ $outTableStartColumn + 3 ] <= $expCutoff )
            && ( $field[ $outTableStartColumn + 4 ] > $expCutoff && $field[ $outTableStartColumn + 5 ] > $expCutoff ) );
    }
    close EV;

    # output
    open SU, ">$groupSummary" || die;
    print SU "Sample1\tSample2\tAS Pattern\tAS Event Type\tAS Event Group Num\tGene with AS Event\tTranscript with AS Event\t";
    print SU "AS Event Group QuantityDiff\tGene with AS Event Group QuantityDiff\tTranscript with AS Event Group QuantityDiff\t";
    print SU "Sample1 Specific AS Event Group\tSample2 Specific AS Event Group\n";
    foreach my $sp ( keys %summary ) {
        my @pair = split /\:\:/, $sp;

        foreach my $asp ( keys %{ $summary{$sp} } ) {
            print SU $pair[0], "\t", $pair[1], "\t", $asp, "\t";
            print SU $summary{$sp}->{$asp}->{'type'}, "\t";
            print SU $summary{$sp}->{$asp}->{'as'},   "\t";
            print SU scalar( keys %{ $summary{$sp}->{$asp}->{'gene'} } ), "\t";
            print SU scalar( keys %{ $summary{$sp}->{$asp}->{'tr'} } ),   "\t";
            print SU $summary{$sp}->{$asp}->{'asdiff'}, "\t";
            print SU scalar( keys %{ $summary{$sp}->{$asp}->{'genediff'} } ), "\t";
            print SU scalar( keys %{ $summary{$sp}->{$asp}->{'trdiff'} } ),   "\t";
            print SU $summary{$sp}->{$asp}->{'s1specific'}, "\t";
            print SU $summary{$sp}->{$asp}->{'s2specific'}, "\n";
        }
    }
    close SU;

    return 0;
}

sub eventSummary {
    my ( $eventDiffFile, $eventSummary, $qvalueCutoff, $expCutoff, $outTableStartColumn ) = @_;
    open EV, "$eventDiffFile" || die;

    my %summary;
    while (<EV>) {
        $_ =~ s/\s+$//;
        next if ( $_ eq '' );
        next if ( $_ =~ /^Gene\t/ );
        my @field = split /\t/, $_;

        my $pair = $field[$outTableStartColumn] . "::" . $field[ $outTableStartColumn + 1 ];
        if ( exists $summary{$pair}->{ $field[5] } ) {
            $summary{$pair}->{ $field[5] }->{'gene'}->{ $field[0] }          = 1;
            $summary{$pair}->{ $field[5] }->{'tr'}->{"$field[0]::$field[3]"} = 1;
            $summary{$pair}->{ $field[5] }->{'tr'}->{"$field[0]::$field[4]"} = 1;
            $summary{$pair}->{ $field[5] }->{'as'}++;
        }
        else {
            $summary{$pair}->{ $field[5] } = {
                'gene' => { $field[0] => 1 },
                'tr'         => { "$field[0]::$field[3]" => 1, "$field[0]::$field[4]" => 1 },
                'as'         => 1,
                'type'       => $field[6],
                'asdiff'     => 0,
                'genediff'   => {},
                'trdiff'     => {},
                's1specific' => 0,
                's2specific' => 0,
            };
        }

        # check exp diff
        if (
               $field[ $outTableStartColumn + 7 ] ne '-'
            && $field[ $outTableStartColumn + 7 ] <= $qvalueCutoff
            && $field[ $outTableStartColumn + 2 ] >= $expCutoff
            && $field[ $outTableStartColumn + 3 ] >= $expCutoff
            && $field[ $outTableStartColumn + 4 ] >= $expCutoff
            && $field[ $outTableStartColumn + 5 ] >= $expCutoff
          )
        {
            $summary{$pair}->{ $field[5] }->{'asdiff'}++;
            $summary{$pair}->{ $field[5] }->{'genediff'}->{ $field[0] }          = 1;
            $summary{$pair}->{ $field[5] }->{'trdiff'}->{"$field[0]::$field[3]"} = 1;
            $summary{$pair}->{ $field[5] }->{'trdiff'}->{"$field[0]::$field[4]"} = 1;
        }

        # check sample specific as event
        $summary{$pair}->{ $field[5] }->{'s1specific'}++
          if ( ( $field[ $outTableStartColumn + 2 ] > $expCutoff && $field[ $outTableStartColumn + 3 ] > $expCutoff )
            && ( $field[ $outTableStartColumn + 4 ] <= $expCutoff || $field[ $outTableStartColumn + 5 ] <= $expCutoff ) );
        $summary{$pair}->{ $field[5] }->{'s2specific'}++
          if ( ( $field[ $outTableStartColumn + 2 ] <= $expCutoff || $field[ $outTableStartColumn + 3 ] <= $expCutoff )
            && ( $field[ $outTableStartColumn + 4 ] > $expCutoff && $field[ $outTableStartColumn + 5 ] > $expCutoff ) );
    }
    close EV;

    # output
    open SU, ">$eventSummary" || die;
    print SU "Sample1\tSample2\tAS Pattern\tAS Event Type\tAS Event Num\tGene with AS Event\tTranscript with AS Event\t";
    print SU "AS Event QuantityDiff\tGene with AS Event QuantityDiff\tTranscript with AS Event QuantityDiff\t";
    print SU "Sample1 Specific AS Event\tSample2 Specific AS Event\n";
    foreach my $sp ( keys %summary ) {
        my @pair = split /\:\:/, $sp;

        foreach my $asp ( keys %{ $summary{$sp} } ) {
            print SU $pair[0], "\t", $pair[1], "\t", $asp, "\t";
            print SU $summary{$sp}->{$asp}->{'type'}, "\t";
            print SU $summary{$sp}->{$asp}->{'as'},   "\t";
            print SU scalar( keys %{ $summary{$sp}->{$asp}->{'gene'} } ), "\t";
            print SU scalar( keys %{ $summary{$sp}->{$asp}->{'tr'} } ),   "\t";
            print SU $summary{$sp}->{$asp}->{'asdiff'}, "\t";
            print SU scalar( keys %{ $summary{$sp}->{$asp}->{'genediff'} } ), "\t";
            print SU scalar( keys %{ $summary{$sp}->{$asp}->{'trdiff'} } ),   "\t";
            print SU $summary{$sp}->{$asp}->{'s1specific'}, "\t";
            print SU $summary{$sp}->{$asp}->{'s2specific'}, "\n";
        }
    }
    close SU;

    return 0;
}

sub processGroup {
    my ( $trExpHR, $samplesAR, $groupFile, $groupDiffFile ) = @_;

    my @asGroupFile;
    open AG, "$groupFile" || die;
    @asGroupFile = <AG>;
    close AG;

    # cal pvalue
    my $title;
    my %storePQvalue;    # store pvalue and calculate qvalue
    my %storeExpData;    # store four exp data, for output step

    # title
    $asGroupFile[0] =~ s/\s+$//;
    $title = $asGroupFile[0];
    $title .= "\tSample1\tSample2\tTr1_Sample1_Exp\tTr2_Sample1_Exp\tTr1_Sample2_Exp\tTr2_Sample2_Exp\tp-value\tq-value";

    # start column of the output, for summary
    my @tmpCol = split /\t/, $asGroupFile[0];
    my $outTableStartColumn = @tmpCol;

    # sample group pair compare
    for ( my $s1 = 0 ; $s1 < @{$samplesAR} ; $s1++ ) {
        for ( my $s2 = $s1 + 1 ; $s2 < @{$samplesAR} ; $s2++ ) {

            for ( my $i = 1 ; $i < @asGroupFile ; $i++ ) {
                $asGroupFile[$i] =~ s/\s+$//;
                next if ( $asGroupFile[$i] eq '' );
                my @field = split /\t/, $asGroupFile[$i];

                # check exp info
                # it has been checked in &processEvent()

                # stat
                my $n11 = &expSum( $samplesAR->[$s1], $field[0], $field[8], $trExpHR );    # tr1 list, sample 1
                my $n21 = &expSum( $samplesAR->[$s1], $field[0], $field[9], $trExpHR );    # tr2 list, sample 1
                my $n12 = &expSum( $samplesAR->[$s2], $field[0], $field[8], $trExpHR );    # tr1, sample 2
                my $n22 = &expSum( $samplesAR->[$s2], $field[0], $field[9], $trExpHR );    # tr2, sample 2

                # store four exp data, for output step
                my $samplePair = $samplesAR->[$s1] . '::' . $samplesAR->[$s2];
                $storeExpData{$samplePair}->[$i] = [ $n11, $n21, $n12, $n22 ];

                # if four number contain >= 2 zero, pass
                my $zeroNum = 0;
                $zeroNum++ if ( $n11 == 0 );
                $zeroNum++ if ( $n21 == 0 );
                $zeroNum++ if ( $n12 == 0 );
                $zeroNum++ if ( $n22 == 0 );
                next       if ( $zeroNum >= 2 );

                # store good data to do stat test

                push @{ $storePQvalue{$samplePair}->{'data'} }, ( $n11, $n21, $n12, $n22 );
                push @{ $storePQvalue{$samplePair}->{'i'} }, $i;
            }

        }
    }

    # cal pvalue and qvalue
    # output statistical test
    open OE, ">$groupDiffFile" || die;
    print OE $title, "\n";

    foreach my $sp ( keys %storeExpData ) {
        if ( exists $storePQvalue{$sp}->{'data'} ) {
            &lib::statistics::ChiSquareTest( $storePQvalue{$sp}->{'data'}, $storePQvalue{$sp} );
        }

        #print Dumper(\%storeQvalue);
        my @pair = split /\:\:/, $sp;
        my $j = 0;    # q value array position
        for ( my $i = 1 ; $i < @asGroupFile ; $i++ ) {
            $asGroupFile[$i] =~ s/\s+$//;
            next if ( $asGroupFile[$i] eq '' );

            my @field = split /\t/, $asGroupFile[$i];

            print OE $asGroupFile[$i], "\t";
            print OE $pair[0], "\t", $pair[1], "\t";
            print OE join( "\t", @{ $storeExpData{$sp}->[$i] }[ 0 .. 3 ] ), "\t";

            if ( exists $storePQvalue{$sp}->{'data'} && $storePQvalue{$sp}->{'i'}->[$j] == $i ) {

                print OE $storePQvalue{$sp}->{'pvalue'}->[$j], "\t";
                print OE $storePQvalue{$sp}->{'qvalue'}->[$j], "\n";
                $j++;
            }
            else {
                print OE "-\t-\n";
            }
        }
    }

    close OE;

    return $outTableStartColumn;

}

sub expSum {
    my ( $sample, $geneId, $trIdList, $trExpHR ) = @_;
    my @trs = split /,/, $trIdList;
    my $sum = 0;
    foreach my $trId (@trs) {
        $sum += $trExpHR->{$sample}->{ $geneId . '::' . $trId };
    }
    return $sum;
}

sub processEvent {
    my ( $trExpHR, $samplesAR, $eventFile, $eventDiffFile ) = @_;

    my @asEventFile;
    open AE, "$eventFile" || die;
    @asEventFile = <AE>;
    close AE;

    # cal pvalue
    my $title;

    # store pvalue and calculate qvalue
    my %storePQvalue;

    # store exp data, for output step
    my %storeExpData;

    # title
    $asEventFile[0] =~ s/\s+$//;
    $title = $asEventFile[0];
    $title .= "\tSample1\tSample2\tTr1_Sample1_Exp\tTr2_Sample1_Exp\tTr1_Sample2_Exp\tTr2_Sample2_Exp\tp-value\tq-value";

    # start column of the output, for summary
    my @tmpCol = split /\t/, $asEventFile[0];
    my $outTableStartColumn = @tmpCol;

    # sample pair compare
    for ( my $s1 = 0 ; $s1 < @{$samplesAR} ; $s1++ ) {
        for ( my $s2 = $s1 + 1 ; $s2 < @{$samplesAR} ; $s2++ ) {

            for ( my $i = 1 ; $i < @asEventFile ; $i++ ) {
                $asEventFile[$i] =~ s/\s+$//;
                next if ( $asEventFile[$i] eq '' );
                my @field = split /\t/, $asEventFile[$i];

                # check exp info
                if (   !exists $trExpHR->{ $samplesAR->[$s1] }->{"$field[0]::$field[3]"}
                    || !exists $trExpHR->{ $samplesAR->[$s1] }->{"$field[0]::$field[4]"}
                    || !exists $trExpHR->{ $samplesAR->[$s2] }->{"$field[0]::$field[3]"}
                    || !exists $trExpHR->{ $samplesAR->[$s2] }->{"$field[0]::$field[4]"} )
                {
                    &logError(
                        "No expression information for sample "
                          . $samplesAR->[$s1] . '/'
                          . $samplesAR->[$s2]
                          . " transcript $field[0]::$field[3]/$field[0]::$field[4]",
                        __FILE__, __LINE__
                    );
                }

                # stat
                my $n11 = $trExpHR->{ $samplesAR->[$s1] }->{"$field[0]::$field[3]"};    # tr1, sample 1
                my $n21 = $trExpHR->{ $samplesAR->[$s1] }->{"$field[0]::$field[4]"};    # tr2, sample 1
                my $n12 = $trExpHR->{ $samplesAR->[$s2] }->{"$field[0]::$field[3]"};    # tr1, sample 2
                my $n22 = $trExpHR->{ $samplesAR->[$s2] }->{"$field[0]::$field[4]"};    # tr2, sample 2

                # store four exp data, for output step
                my $samplePair = $samplesAR->[$s1] . '::' . $samplesAR->[$s2];
                $storeExpData{$samplePair}->[$i] = [ $n11, $n21, $n12, $n22 ];

                # if four number contain >= 2 zero, pass
                my $zeroNum = 0;
                $zeroNum++ if ( $n11 == 0 );
                $zeroNum++ if ( $n21 == 0 );
                $zeroNum++ if ( $n12 == 0 );
                $zeroNum++ if ( $n22 == 0 );
                next       if ( $zeroNum >= 2 );

                # store good data to do stat test
                push @{ $storePQvalue{$samplePair}->{'data'} }, ( $n11, $n21, $n12, $n22 );
                push @{ $storePQvalue{$samplePair}->{'i'} }, $i;
            }

        }
    }

    # cal pvalue and qvalue
    # output statistical test
    open OE, ">$eventDiffFile" || die;
    print OE $title, "\n";

    foreach my $sp ( keys %storeExpData ) {

        # check whether there are data to do test
        if ( exists $storePQvalue{$sp}->{'data'} ) {
            &lib::statistics::ChiSquareTest( $storePQvalue{$sp}->{'data'}, $storePQvalue{$sp} );
        }

        #print Dumper(\%storeQvalue);
        my @pair = split /\:\:/, $sp;
        my $j = 0;    # q value array position
        for ( my $i = 1 ; $i < @asEventFile ; $i++ ) {
            $asEventFile[$i] =~ s/\s+$//;
            next if ( $asEventFile[$i] eq '' );

            my @field = split /\t/, $asEventFile[$i];

            print OE $asEventFile[$i], "\t";
            print OE $pair[0], "\t", $pair[1], "\t";
            print OE join( "\t", @{ $storeExpData{$sp}->[$i] }[ 0 .. 3 ] ), "\t";

            # check whether this line was tested
            if ( exists $storePQvalue{$sp}->{'data'} && $storePQvalue{$sp}->{'i'}->[$j] == $i ) {
                print OE $storePQvalue{$sp}->{'pvalue'}->[$j], "\t";
                print OE $storePQvalue{$sp}->{'qvalue'}->[$j], "\n";
                $j++;
            }
            else {
                print OE "-\t-\n";
            }
        }
    }

    close OE;

    return $outTableStartColumn;
}

sub inputExp {
    my ( $expFile, $trExpHR, $samplesAR ) = @_;

    open EXP, "$expFile" || die;
    my $line = <EXP>;
    $line =~ s/\R//g;
    $line =~ s/\s+$//g;
    my @sample = split /\t/, $line;

    # exp file should contain more than 2 samples
    if ( @sample <= 3 ) {
        &logError( "Transcript expression file should contain >= two samples!", __FILE__, __LINE__ );
    }

    @{$samplesAR} = @sample[ 2 .. $#sample ];

    while ( $line = <EXP> ) {
        $line =~ s/\R//g;
        $line =~ s/\s+$//g;
        my @field = split /\t/, $line;

        for ( my $i = 2 ; $i < @field ; $i++ ) {
            $trExpHR->{ $sample[$i] }->{"$field[0]::$field[1]"} = $field[$i];
        }
    }

    close EXP;
    return 0;
}

