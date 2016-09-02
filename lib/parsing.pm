#
#
package lib::parsing;

=head1 NAME



=head1 SYNOPSIS



=head1 AUTHOR

Zhigang Li  lzg0063(at)126.com  2013-11-28

=head1 DESCRIPTION



=cut

use 5.010;
use strict;
#use warnings;
use lib::basic qw(logMsg logError logWarn logStd);

require Exporter;

our @ISA         = qw(Exporter);
our @EXPORT      = qw();
our @EXPORT_OK   = qw();
our %EXPORT_TAGS = ();
our $VERSION     = 0.1;

################################################################################

=head2 trStartEndStrand

    Title   : trStartEndStrand
    Function: Transcripts gene start, end postion and strand.
    Usage   : ($start,$end, $chr,$strand) = &trStartEndStrand(\%trInfo)
    Returns : Gene start and end postion.
    Args    : %trInfo: gene information in gtf file. if %gtfInfo is the output of subroutine &inputGtfFile(),
                       \%trInfo = $gtfInfo{$gId}->{$trId}.
    Des     :
    
=cut

sub trStartEndStrand {
    my ($trInfoHR) = @_;
    my ( $start, $end, $strand );    #gene start and end position

    foreach my $exonStart ( keys %{$trInfoHR} ) {
        my $exonEnd = $trInfoHR->{$exonStart}->{end};
        if ( defined $start ) {
            $start = ( $start > $exonStart ) ? $exonStart : $start;
            $end   = ( $end < $exonEnd )     ? $exonEnd   : $end;
        }
        else {
            $start  = $exonStart;
            $end    = $exonEnd;
            $strand = $trInfoHR->{$exonStart}->{strand};
        }
    }
    return ( $start, $end, $strand );
}

################################################################################

=head2 inputCDSInGtfFile

    Title   : inputCDSInGtfFile
    Function: Parse CDS info a gtf file into a hash
    Usage   : $sig = &inputGtfFile($in,\%out)
    Returns : 0: ok
              1: no feature
    Args    : $in : a gtf file
              %out:   all CDS feature of gtf file
    Des     : %out = ( 
                        geneId => {
                                transcriptId => {
                                        cdsStart = {
                                            'seqid' => 'Chr1',
                                            'source' => 'TAIR10',
                                            'type' => 'gene',
                                            'start' => '3631',
                                            'end' => '5899',
                                            'score' => '.',
                                            'strand' => '+',
                                            'phase' => '.',
                                            'attributes' => {
                                                    'ID' => 'AT1G01010',
                                                    'Name' => 'AT1G01010',
                                                    'Note' => 'protein_coding_gene'
                                            }
                                        }
                                }
                            } 
                        )
    
=cut

sub inputCDSInGtfFile {
    my ( $gtfFile, $gtfHashRef ) = @_;

    open GF, "$gtfFile" || die "Cannot open gtf file $!\n";

    #parse each line
    while ( my $line = <GF> ) {
        last if ( $line =~ /^##FASTA/ );
        next if ( $line =~ /^#/ );
        my %oneLineFeature;
        my $sig = &gtfFeature( $line, \%oneLineFeature );
        next if ( $sig == 1 );
        next if ( $oneLineFeature{type} ne 'CDS' && $oneLineFeature{type} ne 'cds' );    #select exon only
        my $geneId    = $oneLineFeature{attributes}->{gene_id};
        my $trId      = $oneLineFeature{attributes}->{transcript_id};
        my $exonStart = $oneLineFeature{start};
        $gtfHashRef->{$geneId}->{$trId}->{$exonStart} = \%oneLineFeature;
    }
    close GF;

    return 0;
}

################################################################################

=head2 inputMultiGtfFile

    Title   : inputMultiGtfFile
    Function: Parse gtf files. To distinguish different samples, transcript id will be added a prefix/label.
    Usage   : $sig = &inputGtfFile(\%files,\%out)
    Returns : 0: ok
              1: no feature
    Args    : %files: gtf files
              %out:   all exon feature of gtf file
    Des     : %files =  ( $sampleLabel => $gtfFilePath)
              %out = ( 
                        geneId => {
                                transcriptId => {
                                        exonStart = {
                                            'seqid' => 'Chr1',
                                            'source' => 'TAIR10',
                                            'type' => 'gene',
                                            'start' => '3631',
                                            'end' => '5899',
                                            'score' => '.',
                                            'strand' => '+',
                                            'phase' => '.',
                                            'attributes' => {
                                                    'ID' => 'AT1G01010',
                                                    'Name' => 'AT1G01010',
                                                    'Note' => 'protein_coding_gene'
                                            }
                                        }
                                }
                            } 
                        )
    
=cut

sub inputMultiGtfFile {
    my ( $gtfFiles, $gtfHashRef ) = @_;

    foreach my $label ( keys %{$gtfFiles} ) {

        open GF, $gtfFiles->{$label} || die "Cannot open gtf file $!\n";

        #parse each line
        while ( my $line = <GF> ) {
            last if ( $line =~ /^##FASTA/ );
            next if ( $line =~ /^#/ );
            my %oneLineFeature;
            my $sig = &gtfFeature( $line, \%oneLineFeature );
            next if ( $sig == 1 );
            next if ( $oneLineFeature{type} ne 'exon' );    #select exon only
            my $geneId = $oneLineFeature{attributes}->{gene_id};
            my $trId   = $oneLineFeature{attributes}->{transcript_id};
            $trId = $label . '::' . $trId;
            $oneLineFeature{attributes}->{transcript_id} = $trId;
            my $exonStart = $oneLineFeature{start};
            $gtfHashRef->{$geneId}->{$trId}->{$exonStart} = \%oneLineFeature;
        }
        close GF;
    }

    return 0;
}

################################################################################

=head2 bit2geneInfo

    Title   : bit2geneInfo
    Function: Transform bits to an artificial gene information, whose format is the same as inputGtfFile:
    Usage   : bit2geneInfo( $bit, \%artificialGeneInfo );
    Returns : 
    Args    : $bit: e.g. 10001,--101
    Des     :
    
=cut

sub bit2geneInfo {
    my ( $bitPair, $geneInfoHR ) = @_;
    my @part = split /,/, $bitPair;

    my $unit = 50;

    for ( my $k = 0 ; $k < @part ; $k++ ) {

        my $trId = 'transcript' . $k;
        my @bit = split //, $part[$k];

        # connect the same exon
        my ( $tmpStart, $tmpEnd );
        my $connect = 0;
        for ( my $i = 0 ; $i < @bit ; $i++ ) {
            if ( $bit[$i] eq 1 && $connect eq 0 ) {
                $tmpStart = $i * $unit + 1;
                $tmpEnd   = $i * $unit + $unit;
                $connect  = 1;
            }
            elsif ( $bit[$i] eq 1 && $connect eq 1 ) {
                $tmpEnd = $i * $unit + $unit;
            }

            if ( $connect == 1 ) {
                if ( $bit[$i] eq 0 || $i == $#bit ) {
                    $geneInfoHR->{$trId}->{$tmpStart} = {
                        'seqid'  => 'ASP',
                        'source' => 'ASP',
                        'type'   => 'gene',
                        'start'  => $tmpStart,
                        'end'    => $tmpEnd,
                        'score'  => '.',
                        'strand' => '+',
                        'phase'  => '.',

                    };
                }

                $connect = 0 if ( $bit[$i] eq 0 );
            }
        }
    }
}

################################################################################

=head2 geneChr

    Title   : geneChr
    Function: 
    Usage   : 
    Returns : 
    Args    : 
    Des     :
    
=cut

sub geneChr {
    my ($geneInfoHR) = @_;
    my $chr;
    my ($strand);    #gene start and end position

    foreach my $trId ( keys %{$geneInfoHR} ) {
        foreach my $exonStart ( keys %{ $geneInfoHR->{$trId} } ) {
            $chr = $geneInfoHR->{$trId}->{$exonStart}->{seqid};
            last;
        }
        last;
    }
    return ($chr);
}

################################################################################

=head2 geneStrand

    Title   : geneStrand
    Function: Get gene strand.
    Usage   : ($strand) = &geneStrand(\%geneInfo)
    Returns : $strand eq '+','-', or '0'
    Args    : %geneInfo: gene information in gtf file, which is values of the output of subroutine &inputGtfFile().
    Des     :
    
=cut

sub geneStrand {
    my ($geneInfoHR) = @_;
    my ($strand);    #gene start and end position

    foreach my $trId ( keys %{$geneInfoHR} ) {
        foreach my $exonStart ( keys %{ $geneInfoHR->{$trId} } ) {
            $strand = $geneInfoHR->{$trId}->{$exonStart}->{strand};
            last;
        }
        last;
    }
    return ($strand);
}

################################################################################

=head2 geneStartEnd

    Title   : geneStartEnd
    Function: Get gene start, end postion and strand.
    Usage   : ($start,$end) = &geneStartEnd(\%geneInfo)
    Returns : Gene start and end postion.
    Args    : %geneInfo: gene information in gtf file, which is values of the output of subroutine &inputGtfFile().
    Des     :
    
=cut

sub geneStartEnd {
    my ($geneInfoHR) = @_;
    my ( $start, $end );    #gene start and end position

    foreach my $trId ( keys %{$geneInfoHR} ) {
        foreach my $exonStart ( keys %{ $geneInfoHR->{$trId} } ) {
            my $exonEnd = $geneInfoHR->{$trId}->{$exonStart}->{end};
            if ( defined $start ) {
                $start = ( $start > $exonStart ) ? $exonStart : $start;
                $end   = ( $end < $exonEnd )     ? $exonEnd   : $end;
            }
            else {
                $start = $exonStart;
                $end   = $exonEnd;
            }

        }
    }
    return ( $start, $end );
}

################################################################################

=head2 inputGtfFile

    Title   : inputGtfFile
    Function: Parse a gtf file into a hash
    Usage   : $sig = &inputGtfFile($in,\%out)
    Returns : 0: ok
              1: no feature
    Args    : $in : a gtf file
              %out:   all exon feature of gtf file
    Des     : %out = ( 
                        geneId => {
                                transcriptId => {
                                        exonStart = {
                                            'seqid' => 'Chr1',
                                            'source' => 'TAIR10',
                                            'type' => 'gene',
                                            'start' => '3631',
                                            'end' => '5899',
                                            'score' => '.',
                                            'strand' => '+',
                                            'phase' => '.',
                                            'attributes' => {
                                                    'ID' => 'AT1G01010',
                                                    'Name' => 'AT1G01010',
                                                    'Note' => 'protein_coding_gene'
                                            }
                                        }
                                }
                            } 
                        )
    
=cut

sub inputGtfFile {
    my ( $gtfFile, $gtfHashRef ) = @_;

    open GF, "$gtfFile" || die "Cannot open gtf file $!\n";
    #parse each line
    while ( my $line = <GF> ) {
        last if ( $line =~ /^##FASTA/ );
        next if ( $line =~ /^#/ );
        my %oneLineFeature;
        my $sig = &gtfFeature( $line, \%oneLineFeature );
        next if ( $sig == 1 );
        next if ( $oneLineFeature{type} ne 'exon' && $oneLineFeature{type} ne 'EXON' );    #select exon only
        my $geneId    = $oneLineFeature{attributes}->{gene_id};
        my $trId      = $oneLineFeature{attributes}->{transcript_id};
        my $exonStart = $oneLineFeature{start};
        $gtfHashRef->{$geneId}->{$trId}->{$exonStart} = \%oneLineFeature;
    }
    close GF;
    return 0;
}

################################################################################

=head2 gtfFeature

    Title   : gtfFeature
    Function: Parse one line of gff file
    Usage   : $sig = &gtfFeature($in,\%out)
    Returns : 0: feature is ok
              1: no feature
    Args    : $in : one line of gff file
              \%out:    feature output
    Des     : \% = {
                        'seqid' => 'Chr1'
                        'source' => 'TAIR10',
                        'type' => 'gene',
                        'start' => '3631',
                        'end' => '5899',
                        'score' => '.',
                        'strand' => '+',
                        'phase' => '.',
                        'attributes' => {
                                       'ID' => 'AT1G01010',
                                       'Name' => 'AT1G01010',
                                       'Note' => 'protein_coding_gene'
                        }
    
=cut

sub gtfFeature {
    my ( $gffOneLine, $featureHashRef ) = @_;

    if ( $gffOneLine =~ /(?:^\s*#)|(?:^\s*$)/ ) {
        return 1;
    }

    $gffOneLine =~ s/\R//g;
    $gffOneLine =~ s/#[\s\S]+$//;
    my @field = split /\s*\t\s*/, $gffOneLine;

    # check
    if ( $gffOneLine !~ /\S/ ) {
        return 1;
    }
    if ( @field != 9 ) {
        &logError( "Error line: \n\t$gffOneLine\tnot 9 columns\nExit 1", __FILE__, __LINE__ );
    }
    if ( $field[3] > $field[4] ) {
        &logError( "$gffOneLine\tstart > end", __FILE__, __LINE__ );
    }

    # features of a line
    %{$featureHashRef} = (
        'seqid'      => $field[0],
        'source'     => $field[1],
        'type'       => $field[2],
        'start'      => $field[3],
        'end'        => $field[4],
        'score'      => $field[5],
        'strand'     => $field[6],
        'phase'      => $field[7],
        'attributes' => {},
    );

    # attributes
    my @tags = split /\s*;\s*/, $field[8];
    foreach my $att (@tags) {
        my ( $attName, $attValue );
        if ( $att =~ /^\s*(\S+)\s+['"](.+)['"]\s*$/ ) {    #gff2 or gtf
            ( $attName, $attValue ) = ( $1, $2 );
            $attValue =~ s/^['"]|['"]$//g;
        }
        elsif ( $att =~ /=/ ) {                            #gff3
            ( $attName, $attValue ) = split /=/, $att;
        }
        else {
            &logError( "Error line: \n\t$gffOneLine\nnot 9 columns\nExit 1", __FILE__, __LINE__ );
        }

        # set attributes
        $attName =~ s/^\s+|\s+$//g;
        $attValue =~ s/^\s+|\s+$//g;
        $featureHashRef->{attributes}->{$attName} = $attValue;
    }

    # return
    return 0;
}

1;
__END__
