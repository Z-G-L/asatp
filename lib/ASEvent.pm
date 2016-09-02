#
#
package lib::ASEvent;

=head1 NAME



=head1 SYNOPSIS



=head1 AUTHOR

Zhigang Li	lzg0063(at)126.com	2013-12-31

=head1 DESCRIPTION



=cut

use 5.010;
use strict;

#use warnings;
use lib::basic qw(logMsg logError logWarn logStd);
use lib::math;
#use Data::Dumper;


require Exporter;

our @ISA         = qw(Exporter);
our @EXPORT      = qw();
our @EXPORT_OK   = qw();
our %EXPORT_TAGS = ();
our $VERSION     = 0.1;

################################################################################

=head2 asp2bin

    Title   : asp2bin
    Function: Transform ASP code to spanning bin matrix
    Usage   : &asp2bin($asp)
    Returns : ($tr1Bin,$tr2Bin)
    Args    : 
    Des     :
    
=cut

sub asp2bin {
    my ($asp) = @_;
    my ( $first, $second );
    if ( $asp =~ /^\d+[FT]\d+$/ ) {
        my ( $aa, $bb ) = split /[FT]/, $asp;

        my $binvalue = &lib::math::dec2bin($aa);
        #print STDERR $binvalue;
        
        my $length = length($binvalue);
        if ( $length % 2 == 1 ) {
            # to omit 10--,1101
            # to omit 101,-11
            if ( $bb != 1 && $bb != ($length+1)/2) {
                $binvalue = '0' . $binvalue;
                $length   = length($binvalue);
            }
        }

        $first = substr( $binvalue, 0, $bb );
        $second = substr( $binvalue, $bb );

        my ( $lf, $ls ) = ( length($first), length($second) );

        if ( $asp =~ /F/ ) {
            if ( $lf < $ls ) {
                $first = ( '-' x abs( $ls - $lf ) ) . $first;
            }
            elsif ( $lf > $ls ) {
                $second = ( '-' x abs( $ls - $lf ) ) . $second;
            }
            else {
                &logError( "Error ASP format!  Please try again!\n", __FILE__, __LINE__ );
            }

            $first  .= '1';
            $second .= '1';
        }
        elsif ( $asp =~ /T/ ) {

            if ( $lf < $ls ) {
                $first = $first . ( '-' x abs( $ls - $lf ) );
            }
            elsif ( $lf > $ls ) {
                $second = $second . ( '-' x abs( $ls - $lf ) );
            }
            else {
                &logError( "Error ASP format!  Please try again!\n", __FILE__, __LINE__ );
            }

            $first  = '1' . $first;
            $second = '1' . $second;
        }

    }
    elsif ( $asp =~ /^\d+$/ ) {

        my $binvalue = &lib::math::dec2bin($asp);

        my $length = length($binvalue);
        if ( $length % 2 == 1 ) {
            $binvalue = '0' . $binvalue;
            $length   = length($binvalue);
        }

        $first = substr( $binvalue, 0, $length / 2 );
        $second = substr( $binvalue, $length / 2 );

        $first  = '1' . $first . '1';
        $second = '1' . $second . '1';

    }
    else {
        &logError( "Error ASP format!  Please try again!\n", __FILE__, __LINE__ );
    }

    return ( $first, $second );
}

################################################################################

=head2 bitCheck

    Title   : bitCheck
    Function: Check whether the bit pair is a AS span unit.
    Usage   : 
    Returns : 
    Args    : $bitPair: e.g. 10001,--101
    Des     :
    
=cut

sub bitCheck {
    my ($bitPair) = @_;
    my @part = split /,/, $bitPair;
    my @bit1 = split //,  $part[0];
    my @bit2 = split //,  $part[1];

    if ( @bit1 != @bit2 ) {
        return 1;
    }

    for ( my $i = 1 ; $i < @bit1 - 1 ; $i++ ) {
        if ( $bit1[$i] eq '1' && $bit2[$i] eq '1' ) {
            return 1;
        }
    }

    for ( my $i = 0 ; $i < @bit1 - 1 ; $i++ ) {
        if ( $bit1[$i] eq $bit1[ $i + 1 ] && $bit2[$i] eq $bit2[ $i + 1 ] ) {
            return 1;
        }
    }

    return 0;
}

################################################################################

=head2 asSpanUnit

    Title   : asSpanUnit
    Function: Get the start and end of the span unit.
    Usage   : 
    Returns : 
    Args    : $eventBitRegion: It has removed exon-exon type.
    Des     : 
    
=cut

sub asSpanUnit {
    my ( $eventBitRegion, $strand ) = @_;

    my @region = split /,/, $eventBitRegion;

    @region = reverse(@region) if ( $strand eq '-' || $strand eq '-1' );

    $region[0] =~ /(\d+)-\d+/;
    my $start = $1;

    $region[$#region] =~ /\d+-(\d+)/;
    my $end = $1;
    return $start . '-' . $end;
}

################################################################################

=head2 ASEventGroup

    Title   : ASEventGroup
    Function: Same AS event in a gene will be classified into the same group
    Usage   : &lib::ASEvent::ASEventGroup($tr1Id, $tr2Id, \%asEvent, \%geneASInfo )
    Returns : 0
    Args    : %asEvent: output of subfunction detectASEvent
              %geneASInfo = {
                      $groupIdentify=>{
                                'groupId'      => $geneASGroupId,
                                'tr1Ids'       => $tr1Id,
                                'tr2Ids'       => $tr2Id,
                                'tr1BitString' => $asEventPairTrHR->{$asNumOfTr}->[2],
                                'tr2BitString' => $asEventPairTrHR->{$asNumOfTr}->[3],
                                'asRegion'     => $asEventPairTrHR->{$asNumOfTr}->[4],
                                'asp'          => $asEventPairTrHR->{$asNumOfTr}->[0],
                                'astype'       => $asEventPairTrHR->{$asNumOfTr}->[1],
                                'strand'       => $asEventPairTrHR->{$asNumOfTr}->[5],
                                'spanUnit'     => $asEventPairTrHR->{$asNumOfTr}->[6],
                      }
              }
    Des     :
    
=cut

sub ASEventGroup {

    my ( $tr1Id, $tr2Id, $asEventPairTrHR, $asGroupHR ) = @_;

    my $geneASGroupIdNum = keys %{$asGroupHR};

    # process each AS event
    foreach my $asNumOfTr ( keys %{$asEventPairTrHR} ) {    # each AS event number
        my $groupIdentify = $asEventPairTrHR->{$asNumOfTr}->[0] . "_" . $asEventPairTrHR->{$asNumOfTr}->[4];

        # get AS groups in a gene
        my $geneASGroupId;
        if ( !exists $asGroupHR->{$groupIdentify} ) {
            $geneASGroupIdNum++;
            $geneASGroupId = 'as_g' . $geneASGroupIdNum;
            $asGroupHR->{$groupIdentify} = {
                'groupId'      => $geneASGroupId,
                'tr1Ids'       => $tr1Id,
                'tr2Ids'       => $tr2Id,
                'tr1BitString' => $asEventPairTrHR->{$asNumOfTr}->[2],
                'tr2BitString' => $asEventPairTrHR->{$asNumOfTr}->[3],
                'asRegion'     => $asEventPairTrHR->{$asNumOfTr}->[4],    # as event region, including span exon-exon
                'asp'          => $asEventPairTrHR->{$asNumOfTr}->[0],
                'astype'       => $asEventPairTrHR->{$asNumOfTr}->[1],
                'strand'       => $asEventPairTrHR->{$asNumOfTr}->[5],
                'spanUnit'     => $asEventPairTrHR->{$asNumOfTr}->[6],    # start and end of bit region
            };
        }
        else {
            my ( $tmpTr1Id, $tmpTr2Id );
            if ( $asEventPairTrHR->{$asNumOfTr}->[2] eq $asGroupHR->{$groupIdentify}->{tr1BitString} ) {
                $tmpTr1Id = $tr1Id;
                $tmpTr2Id = $tr2Id;
            }
            elsif ( $asEventPairTrHR->{$asNumOfTr}->[2] eq $asGroupHR->{$groupIdentify}->{tr2BitString} ) {
                $tmpTr1Id = $tr2Id;
                $tmpTr2Id = $tr1Id;

            }
            else {
                &logError( "AS in the same group wit different bit string:\n $tr1Id\t $tr2Id\n", __FILE__, __LINE__ );
            }

            my $tii = $asGroupHR->{$groupIdentify}->{tr1Ids} . ',';
            $asGroupHR->{$groupIdentify}->{tr1Ids} .= ',' . $tmpTr1Id if ( $tii !~ /$tmpTr1Id,/ );
            $tii = $asGroupHR->{$groupIdentify}->{tr2Ids} . ',';
            $asGroupHR->{$groupIdentify}->{tr2Ids} .= ',' . $tmpTr2Id if ( $tii !~ /$tmpTr2Id,/ );
        }
    }

    return 0;
}

################################################################################

=head2 spanUnitForGraph

    Title   : spanUnitForGraph
    Function: Combining the span unit which belongs the same exon.
    Usage   : &lib::ASEvent::extendASPSpanUnit($spanUnit,$bitString, $strand, \@combine)
    Returns : 
    Args    : $spanUnit: "8455-8813,8253-8454,7554-8252,7295-7553,6693-7294"
              $bitString: 10001
              $strand: 1,-1,+,-
              @combine: the exon region only
    Des     :
    
=cut

sub spanUnitForGraph {
    my ( $spanUnit, $bitString1, $bitString2, $strand, $oCombineUnit ) = @_;
    my @region = split /,/, $spanUnit;
    my @bit1   = split //,  $bitString1;
    my @bit2   = split //,  $bitString2;
    if ( $strand eq '-' || $strand eq '-1' ) {
        @bit1   = reverse(@bit1);
        @bit2   = reverse(@bit2);
        @region = reverse(@region);
    }

    # delete the end span if it's exon-exon type, for better graph
    if ( $bit1[0] eq 1 && $bit2[0] eq 1 ) {
        $region[0] =~ /(\d+)-(\d+)/;
        if ( $2 - $1 > 10 ) {
            $region[0] = ( $2 - 14 ) . '-' . $2;
        }
    }
    if ( $bit1[$#bit1] eq 1 && $bit2[$#bit2] eq 1 ) {
        $region[$#region] =~ /(\d+)-(\d+)/;
        if ( $2 - $1 > 10 ) {
            $region[$#region] = $1 . '-' . ( $1 + 14 );
        }
    }

    # extend span unit and combine units in the same exons
    &logError( "The number of bit != region:\n@region\n@bit1", __FILE__, __LINE__ ) if ( @region != @bit1 );
    my ( $tmpStart, $tmpEnd );
    my $connect = 0;
    for ( my $i = 0 ; $i < @bit1 ; $i++ ) {
        if ( $bit1[$i] eq 1 && $connect eq 0 ) {
            $region[$i] =~ /(\d+)-(\d+)/;
            $tmpStart = $1;
            $tmpEnd   = $2;
            $connect  = 1;
        }
        elsif ( $bit1[$i] eq 1 && $connect eq 1 ) {
            $region[$i] =~ /(\d+)-(\d+)/;
            $tmpEnd = $2;
        }

        if ( $connect == 1 ) {
            if ( $bit1[$i] eq 0 || $i == $#bit1 ) {
                push @{$oCombineUnit}, $tmpStart . '-' . $tmpEnd;
            }

            $connect = 0 if ( $bit1[$i] eq 0 );
        }
    }
    return 0;
}

################################################################################

=head2 detectAlternativeSplicingEvent

    Title   : detectAlternativeSplicingEvent
    Function: Detect AS event between transcrips pairs.
    Usage   : &detectAlternativeSplicingEvent(\%geneInfo,$tr1Id,$tr2Id,\%oTrPairAsEventsHR);
    Returns : 0
    Args    : %geneInfo: gene informations from gtf file. %geneInfo = (trId=>{startOfExon=>{featuresOfGtf=>Values}}), 
              $tr1Id,$tr2Id: transcript ids compared
              %oTrPairAsEventsHR: All AS evernts between a transcript pair. 
                        %oTrPairAsEventsHR:  {
                                                $asEventNumber=>[$asp, 
                                                                $asType, 
                                                                $eventTr1Bit,
                                                                $eventTr2Bit, 
                                                                $eventBitRegion, 
                                                                $geneStrand],
                                                                $asSpanUnit)
                                                }
    Des     :
    
=cut

sub detectASEvent {
    my ( $geneInfoHashRef, $tr1Id, $tr2Id, $oTrPairAsEventsHR ) = @_;
    my $geneStrand = &lib::parsing::geneStrand($geneInfoHashRef);
    my @tr1BitNum;       # bit number of tr1 in comparison unit, exon=1, intron=0, extragenic=-
    my @tr2BitNum;       # # bit number of tr2 in comparison unit, exon=1, intron=0, extragenic=-
    my @bitNumRegion;    # comparison unit between tr1 and tr2
    &transcripPairCompare( $geneInfoHashRef, $tr1Id, $tr2Id, \@tr1BitNum, \@tr2BitNum, \@bitNumRegion );

    # if the gene is on the negative strand
    if ( $geneStrand eq '-' || $geneStrand eq '-1' ) {
        @tr1BitNum    = reverse(@tr1BitNum);
        @tr2BitNum    = reverse(@tr2BitNum);
        @bitNumRegion = reverse(@bitNumRegion);
    }

    ########################################
    # detect AS event and type
    my $eventStart         = 0;             # AS event start bit/region
    my $eventEnd           = 0;             # AS event end bit/region
    my $overlapBit         = 0;             # check each tr overlap region: exon-exon, exon-intron, or others
    my $trOverlapRegionNum = @tr1BitNum;    # the number of tr overlap regions/bits
    my @allASEvent;
    for ( my $i = 0 ; $i < $trOverlapRegionNum ; $i++ ) {
        $overlapBit++ if ( $tr1BitNum[$i] eq '1' || $tr1BitNum[$i] eq '-' );
        $overlapBit++ if ( $tr2BitNum[$i] eq '1' || $tr2BitNum[$i] eq '-' );
        if ( $overlapBit == 2 ) {
            if ( $eventEnd != $eventStart ) {
                $eventEnd = $i;

                # at least one side is exon-exon if the span unit is AS event.
                # (exon--,--exon) type will be removed
                if (   ( $tr1BitNum[$eventStart] eq '1' && $tr2BitNum[$eventStart] eq '1' )
                    || ( $tr1BitNum[$eventEnd] eq '1' && $tr2BitNum[$eventEnd] eq '1' ) )
                {
                    push @allASEvent, [ $eventStart, $eventEnd ];    # each AS event start and end bit/region
                }
            }

            $eventStart = $i;
            $eventEnd   = $eventStart;
        }
        elsif ( $overlapBit == 1 ) {
            $eventEnd++;
        }
        $overlapBit = 0;
    }

    my $asEventNumber = 0;    # the number of AS event in a gene
    for ( my $i = 0 ; $i < @allASEvent ; $i++ ) {
        $asEventNumber = $i + 1;

        #print Dumper(\@allASEvent );
        # Alternative splicing pattern
        # $aspScoreReverseTrId: if 0, bit string corresponding to ASP is tr1,tr2; if 1 is tr2,tr1
        my ( $asp, $asUnitStart, $asUnitEnd, $aspScoreReverseTrId ) = &ASPatternScore( \@tr1BitNum, \@tr2BitNum, $allASEvent[$i] );
        my $asType = &ASTypeCheck($asp);    # Alternative splicing pattern

        # return AS event info
        my $eventTr1Bit = join( '',  @tr1BitNum[ $allASEvent[$i]->[0] .. $allASEvent[$i]->[1] ] );
        my $eventTr2Bit = join( '',  @tr2BitNum[ $allASEvent[$i]->[0] .. $allASEvent[$i]->[1] ] );
        my $eventRegion = join( ',', @bitNumRegion[ $allASEvent[$i]->[0] .. $allASEvent[$i]->[1] ] );
        my $asASPRegion = join( ',', @bitNumRegion[ $asUnitStart .. $asUnitEnd ] );
        my $asSpanUnit = &asSpanUnit( $asASPRegion, $geneStrand );
        $oTrPairAsEventsHR->{$asEventNumber} =
          [ $asp, $asType, $eventTr1Bit, $eventTr2Bit, $eventRegion, $geneStrand, $asSpanUnit, $aspScoreReverseTrId ];
    }

    return 0;
}

################################################################################

=head2 ASPatternScore

    Title   : ASPatternScore
    Function: Calculate alternative splicing pattern score.
    Usage   : 
    Returns : ASP score
    Args    : 
    
=cut

sub ASPatternScore {
    my ( $tr1BitAR, $tr2BitAR, $asStartEndAR ) = @_;
    my $aspScoreReverseTr = 0;
    my ( $start, $end ) = @{$asStartEndAR};
    $start++ if ( $tr1BitAR->[$start] eq 1 && $tr2BitAR->[$start] eq 1 );
    $end--   if ( $tr1BitAR->[$end] eq 1   && $tr2BitAR->[$end] eq 1 );
    &logError( "AS region: start > end ", __FILE__, __LINE__ ) if ( $start > $end );

    # select which is the first/second bit array to calculate asp code
    my $firstBitString;
    my $secondBitString;
    for ( my $i = 0 ; $i < 2 ; $i++ ) {
        if ( $tr1BitAR->[ $start + $i ] eq 1 ) {
            $firstBitString  = join( '', @{$tr1BitAR}[ $start .. $end ] );
            $secondBitString = join( '', @{$tr2BitAR}[ $start .. $end ] );
            last;
        }
        elsif ( $tr2BitAR->[ $start + $i ] eq 1 ) {
            $firstBitString  = join( '', @{$tr2BitAR}[ $start .. $end ] );
            $secondBitString = join( '', @{$tr1BitAR}[ $start .. $end ] );
            $aspScoreReverseTr = 1;
            last;
        }
        &logError( "Bit number error: @{$tr1BitAR}\t@{$tr2BitAR}", __FILE__, __LINE__ ) if ( $i == 1 );
    }

    # check extracgenic AS
    my $extType = '';
    if ( $firstBitString =~ /-/ || $secondBitString =~ /-/ ) {
        my $tmp1String = join( "", @{$tr1BitAR}[ $asStartEndAR->[0] .. $asStartEndAR->[1] ] );
        my $tmp2String = join( "", @{$tr2BitAR}[ $asStartEndAR->[0] .. $asStartEndAR->[1] ] );
        if ( $tmp1String =~ /^-/ || $tmp2String =~ /^-/ ) {
            $firstBitString =~ s/-//g;
            $secondBitString =~ s/-//g;
            $extType = 'F' . length($firstBitString);
        }
        elsif ( $tmp1String =~ /-$/ || $tmp2String =~ /-$/ ) {
            $firstBitString =~ s/-//g;
            $secondBitString =~ s/-//g;
            $extType = 'T' . length($firstBitString);
        }
    }

    my $aspBitString = $firstBitString . $secondBitString;
    # abp score
    my $asp = &lib::math::bin2dec($aspBitString);    
    $asp .= $extType;
    return ( $asp, $start, $end, $aspScoreReverseTr );
}

################################################################################

=head2 ASType

    Title   : ASType
    Function: Define patterns of alternative splicing events.
    Usage   : 
    Returns : AS type
    Args    : 
    
=cut

sub ASTypeCheck {
    my ($aspScore) = @_;
    my $asType;

    if ( $aspScore eq '16' ) {
        $asType = 'CE';    #Cassette exon (skipped exon)
    }
    elsif ( $aspScore eq '2' ) {
        $asType = 'IR';    #Intron retention
    }
    elsif ( $aspScore eq '8' ) {
        $asType = 'A5SS';    #Alternative 5' sites
    }
    elsif ( $aspScore eq '4' ) {
        $asType = 'A3SS';    #Alternative 3' sites
    }
    elsif ( $aspScore eq '258' ) {
        $asType = 'MXE';     #Mutually exclusive exons
    }
    elsif ( $aspScore eq '34F4' ) {
        $asType = 'AFE';     #Alternative first exon
    }
    elsif ( $aspScore eq '17T2' ) {
        $asType = 'ALE';     #Alternative last exon
    }
    else {
        $asType = 'complex';
    }
    return ($asType);
}

################################################################################

=head2 transcripPairCompare

    Title   : transcripPairCompare
    Function: Compare a pair of transcripts and output overlap regions.
    Usage   : &transcripPairCompare( \%geneInfo, $tr1Id, $tr2Id, \@tr1BitNum, \@tr2BitNum, \@bitNumRegion  )
    Returns : 0
    Args    : %geneInfo: gene information. It's same as $gtfFileHash{$geneId}
              $tr1Id, $tr2Id: transcript id
              @tr1BitNum, @tr2BitNum: tr1 and tr2 overlap bits/regions
              @bitNumRegion: genome position of overlap bits/regions
    Des     :
    
=cut

sub transcripPairCompare {
    my ( $geneInfoHashRef, $tr1Id, $tr2Id, $oTr1BitNumAR, $oTr2BitNumAR, $oBitNumRegionAR, ) = @_;

    my ( %boundariesInTr1, %boundariesInTr2, %exonInTr1, %exonInTr2 );

    # %exonInTr1, all positions of exon
    # %boundariesInTr1, boundaries of exons and introns
    my ( $tr1Start, $tr1End ) = &exonIntronBoundary( $geneInfoHashRef->{$tr1Id}, \%boundariesInTr1, \%exonInTr1 );
    my ( $tr2Start, $tr2End ) = &exonIntronBoundary( $geneInfoHashRef->{$tr2Id}, \%boundariesInTr2, \%exonInTr2 );
    my @uniqueBoundary = &uniqueBoundaryOfATranscipt( \%boundariesInTr1, \%boundariesInTr2 );

    # check whether the two tr span overlap
    my $os = ( $tr1Start > $tr2Start ) ? $tr1Start : $tr2Start;
    my $oe = ( $tr1End > $tr2End )     ? $tr2End   : $tr1End;
    return 0 if ( $os >= $oe );

    #compare boudnaries, a boudary is in exon, intron and extragenic region will be expressed as 1,0 and -
    my @tmpTr1Num;
    my @tmpTr2Num;
    foreach my $b (@uniqueBoundary) {    # region end position
        my $tr1Num = 0;                  # tr1 boudnaries number, which will be translate into bit number
        my $tr2Num = 0;                  # tr2 boudnaries number

        # check tr1
        if ( exists $exonInTr1{$b} ) {    # in exon of tr1
            $tr1Num = 1;
        }
        elsif ( $b > $tr1End || $b < $tr1Start ) {    #in extragenic region, others is in intron of tr1
            $tr1Num = '-';
        }

        # check tr2
        if ( exists $exonInTr2{$b} ) {
            $tr2Num = 1;
        }
        elsif ( $b > $tr2End || $b < $tr2Start ) {
            &logError( "Maybe include comparision of extragenic-extragenic!", __FILE__, __LINE__ ) if ( $tr1Num eq '-' );
            $tr2Num = '-';
        }
        push @tmpTr1Num, $tr1Num;
        push @tmpTr2Num, $tr2Num;

        #print Dumper(\@tmpCompare);
    }

    # simplify @compare and @uniqueBoundary;
    # @compare: combine ends of each region as a number, e.g. 2,2 as 2; @uniqueBoundary: start,end --> start-end
    for ( my $i = 0 ; $i < @tmpTr1Num ; $i += 2 ) {
        &logError( "Two ends of a region are not the same!\t$tr1Id, $tr2Id\n@tmpTr1Num\n@tmpTr2Num\n@uniqueBoundary", __FILE__, __LINE__ )
          if ( $tmpTr1Num[$i] ne $tmpTr1Num[ $i + 1 ] || $tmpTr2Num[$i] ne $tmpTr2Num[ $i + 1 ] );
        my $region = $uniqueBoundary[$i] . '-' . $uniqueBoundary[ $i + 1 ];
        push @{$oTr1BitNumAR},    $tmpTr1Num[$i];
        push @{$oTr2BitNumAR},    $tmpTr2Num[$i];
        push @{$oBitNumRegionAR}, $region;
    }

    return 0;

}

################################################################################

=head2 uniqueBoundaryOfATranscipt

    Title   : uniqueBoundaryOfATranscipt
    Function: Unique all exon and intron boudaries in a transcript.
    Usage   : &uniqueBoundaryOfATranscipt(\%boundariesInTr1,\%boundariesInTr2)
    Returns : 
    Args    : 
    Des     :
    
=cut

sub uniqueBoundaryOfATranscipt {
    my ( $tr1Boundaries, $tr2Boundaries ) = @_;
    my %unique = ( %{$tr1Boundaries}, %{$tr2Boundaries} );
    my @tmpUniqueBoundarySorted = sort { $a <=> $b } keys %unique;

    #remove the most first and last position of transcripts, they are the first bp out of gene
    my @uniqueBoundarySorted = @tmpUniqueBoundarySorted[ 1 .. $#tmpUniqueBoundarySorted - 1 ];

    for ( my $u = 0 ; $u < @uniqueBoundarySorted ; $u++ ) {
        $uniqueBoundarySorted[$u] =~ s/\.5|\.3//;
    }

    return (@uniqueBoundarySorted);
}

################################################################################

=head2 exonIntronBoundary

    Title   : exonIntronBoundary
    Function: Set exon bouddaries as 1, whereas intron boudaries as 0.
    Usage   : ($trStart, $trEnd) = &exonIntronBoundary(\%transciptInformation, \%output)
    Returns : 
    Args    : 
    Des     :
    
=cut

sub exonIntronBoundary {
    my ( $trInfo, $boundaryNumHR, $oExonPosiHR ) = @_;
    my ( $trStart, $trEnd );
    my @exonStarts = sort { $a <=> $b } keys %{$trInfo};    #all start positions of a transcript
    for ( my $i = 0 ; $i < @exonStarts ; $i++ ) {
        my $start = $exonStarts[$i];                        #start and end of an exon in a transcript
        my $end   = $trInfo->{$start}->{end};

        # exon base as 1
        map( $oExonPosiHR->{$_} = 1, ( $start .. $end ) );

        # check transcript start and end
        $trStart = $start if ( $i == 0 );
        $trEnd   = $end   if ( $i == $#exonStarts );

        # check exon boudaries  as 1
        my $eFivePrimer  = $start . '.5';
        my $eThreePrimer = $end . '.3';
        $boundaryNumHR->{$eFivePrimer}  = 1;
        $boundaryNumHR->{$eThreePrimer} = 1;

        # check intron boundary as 0
        # other exons, include the head and end of transcript(which are filtered later)
        my $iFivePrimer  = ( $end + 1 ) . '.5';
        my $iThreePrimer = ( $start - 1 ) . '.3';
        $boundaryNumHR->{$iThreePrimer} = 0;
        $boundaryNumHR->{$iFivePrimer}  = 0;
    }
    return ( $trStart, $trEnd );

}

################################################################################

1;
__END__
