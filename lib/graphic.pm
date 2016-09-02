#
#
package lib::graphic;

=head1 NAME



=head1 SYNOPSIS



=head1 AUTHOR

Zhigang Li	lzg0063(at)126.com	2014-3-3

=head1 DESCRIPTION



=cut

use 5.010;
use strict;
#use warnings;
use lib::basic qw(logMsg logError logWarn pbar);
use lib::ASEvent();
use lib::parsing();
use Bio::Graphics;
use Bio::SeqFeature::Generic;
#use Data::Dumper;

require Exporter;

our @ISA         = qw(Exporter);
our @EXPORT      = qw();
our @EXPORT_OK   = qw();
our %EXPORT_TAGS = ();
our $VERSION     = 0.1;

=head2 bitGraph

    Title   : bitGraph
    Function: Create a artificial gene graph. Similar as createGraph function.
    Usage   : 
    Returns : 
    Args    : 
    Des     :
    
=cut

sub bitGraph {

        my ( $geneId, $geneInfoHR, $picFormat, $picFile ) = @_;

        my ( $geneStart, $geneEnd ) = &lib::parsing::geneStartEnd($geneInfoHR);

        ############################################################

        #create panel
        my $panel = &createPanel( $geneStart, $geneEnd, $picFormat );

        ############################################################

        #create a scale track
        my $full_length = Bio::SeqFeature::Generic->new(
                -start => $geneStart - 10,
                -end   => $geneEnd + 10,
        );

        $panel->add_track(
                $full_length,
                -glyph     => 'arrow',
                -tick      => 2,
                -fgcolor   => 'black',
                -northeast => 1,
        );

        ############################################################

        #create gene track
        my $track = $panel->add_track(
                -glyph        => 'segments',
                -label        => 1,
                -strand_arrow => 1,
                -key          => "$geneId",
                -bgcolor      => 'orange',

                #-bump         => 1,
        );

        foreach my $trId ( sort { $a cmp $b } keys %{$geneInfoHR} ) {
                my $feature = Bio::SeqFeature::Generic->new( -display_name => $trId, );
                foreach my $exonStart ( keys %{ $geneInfoHR->{$trId} } ) {
                        #my $strand = 0;
                        #$strand = ( $geneInfoHR->{$trId}->{$exonStart}->{strand} eq '+' ) ? (1) : (-1);
                        my $subFeature = Bio::SeqFeature::Generic->new(
                                -start  => $exonStart,
                                -end    => $geneInfoHR->{$trId}->{$exonStart}->{end},
                                -strand => 0,
                        );
                        $feature->add_sub_SeqFeature( $subFeature, 'EXPAND' );
                }
                $track->add_feature($feature);
        }


        ############################################################
        #create pic
        $picFile .= '.' . $picFormat;
        &createPic( $panel, $picFormat, $picFile );
        $panel->finished;

}

################################################################################

=head2 createGraph

    Title   : 
    Function: 
    Usage   : 
    Returns : 
    Args    : 
    Des     :
    
=cut

sub createGraph {
        my ( $geneId, $geneInfoHR, $geneASInfoHR, $picFormat, $picFile ) = @_;

        my ( $geneStart, $geneEnd ) = &lib::parsing::geneStartEnd($geneInfoHR);

        ############################################################

        #create panel
        my $panel = &createPanel( $geneStart, $geneEnd, $picFormat );

        ############################################################

        #create a scale track
        my $full_length = Bio::SeqFeature::Generic->new(
                -start => $geneStart - 100,
                -end   => $geneEnd + 100,
        );

        $panel->add_track(
                $full_length,
                -glyph     => 'arrow',
                -tick      => 2,
                -fgcolor   => 'black',
                -northeast => 1,
        );

        ############################################################

        #create gene track
        my $track = $panel->add_track(
                -glyph        => 'segments',
                -label        => 1,
                -strand_arrow => 1,
                -key          => "Gene: $geneId",
                -bgcolor      => 'orange',

                #-bump         => 1,
        );

        foreach my $trId ( sort { $a cmp $b } keys %{$geneInfoHR} ) {
                my $feature = Bio::SeqFeature::Generic->new( -display_name => $trId, );
                foreach my $exonStart ( keys %{ $geneInfoHR->{$trId} } ) {
                        my $strand = 0;
                        $strand = ( $geneInfoHR->{$trId}->{$exonStart}->{strand} eq '+' ) ? (1) : (-1);
                        my $subFeature = Bio::SeqFeature::Generic->new(
                                -start  => $exonStart,
                                -end    => $geneInfoHR->{$trId}->{$exonStart}->{end},
                                -strand => $strand,
                        );
                        $feature->add_sub_SeqFeature( $subFeature, 'EXPAND' );
                }
                $track->add_feature($feature);
        }

        ############################################################

        # create AS event tracks
        my $asEventInGeneNum = 0;    # the number of AS event group in a gene, for draw color of AS track
        foreach my $geneASId ( sort { $geneASInfoHR->{$a}->{astype} cmp $geneASInfoHR->{$b}->{astype} } keys %{$geneASInfoHR} ) {    #each AS event
                $asEventInGeneNum++;

                #track color
                my $trackColor = '';
                if ( $asEventInGeneNum % 2 == 0 ) {
                        $trackColor = "gainsboro";
                }
                else {
                        $trackColor = "seashell";
                }

                #track feature
                my $keyName =
                    ' ' x 5
                  . 'GroupId: '
                  . $geneASInfoHR->{$geneASId}->{"groupId"}
                  . ' ' x 5 . "ASP="
                  . $geneASInfoHR->{$geneASId}->{"asp"}
                  . ' ' x 5
                  . 'AS_Type: '
                  . $geneASInfoHR->{$geneASId}->{"astype"};

                my $track = $panel->add_track(
                        -glyph        => 'segments',
                        -label        => 1,
                        -key          => $keyName,
                        -connector    => 'hat',
                        -strand_arrow => 1,
                        -tkcolor      => $trackColor,
                        -height       => 8,
                        -bgcolor      => '#99E64D',     # Grass Green
                        -bump         => 1,

                        #-font=>'Times New Roman',
                );

                my $asStrand = ( $geneASInfoHR->{$geneASId}->{strand} eq '+' ) ? (1) : (-1);
                my $asEndExonExtendLen = 10;            # For the first and the last exon region, only show a small fragment

                # tr1
                my $display1 = $geneASInfoHR->{$geneASId}->{tr1Ids};
                $display1 =~ s/,/, /g;
                my $asFeature1 = Bio::SeqFeature::Generic->new( -display_name => $display1 );

                my @tr1ASExon;
                &lib::ASEvent::spanUnitForGraph(
                        $geneASInfoHR->{$geneASId}->{"asRegion"},
                        $geneASInfoHR->{$geneASId}->{"tr1BitString"},
                        $geneASInfoHR->{$geneASId}->{"tr2BitString"},
                        $asStrand, \@tr1ASExon
                );

                for ( my $i = 0 ; $i < @tr1ASExon ; $i++ ) {
                        $tr1ASExon[$i] =~ /(\d+)-(\d+)/;
                        my $asSubFeature = Bio::SeqFeature::Generic->new(
                                -start  => $1,
                                -end    => $2,
                                -strand => 0,
                        );
                        $asFeature1->add_sub_SeqFeature( $asSubFeature, 'EXPAND' );
                }
                $track->add_feature($asFeature1);

                # tr2
                my $display2 = $geneASInfoHR->{$geneASId}->{tr2Ids};
                $display2 =~ s/,/, /g;
                my $asFeature2 = Bio::SeqFeature::Generic->new( -display_name => $display2 );

                my @tr2ASExon;
                &lib::ASEvent::spanUnitForGraph(
                        $geneASInfoHR->{$geneASId}->{"asRegion"},
                        $geneASInfoHR->{$geneASId}->{"tr2BitString"},
                        $geneASInfoHR->{$geneASId}->{"tr1BitString"},
                        $asStrand, \@tr2ASExon
                );

                for ( my $i = 0 ; $i < @tr2ASExon ; $i++ ) {
                        $tr2ASExon[$i] =~ /(\d+)-(\d+)/;
                        my $asSubFeature = Bio::SeqFeature::Generic->new(
                                -start  => $1,
                                -end    => $2,
                                -strand => 0,
                        );
                        $asFeature2->add_sub_SeqFeature( $asSubFeature, 'EXPAND' );
                }
                $track->add_feature($asFeature2);
        }

        ############################################################
        #create pic
        $picFile .= '.' . $picFormat;
        &createPic( $panel, $picFormat, $picFile );
        $panel->finished;

}
################################################################################

=head2 createPic

    Title   : createPic
    Function: 
    Usage   : 
    Returns : 
    Args    : 
    Des     :
    
=cut

sub createPic {
        my ( $panel, $picFormat, $picFile ) = @_;
        open PF, ">$picFile" || die "cannot open pic file: $picFile\n";
        if ( $picFormat eq 'svg' || $picFormat eq 'SVG' ) {
                print PF $panel->svg();
        }
        elsif ( $picFormat eq 'png' || $picFormat eq 'SVG' ) {
                binmode PF;
                print PF $panel->png();
        }
        close PF;

}

################################################################################

=head2 createPanel

    Title   : createPanel
    Function: 
    Usage   : 
    Returns :
    Args    : 
    Des     :
    
=cut

sub createPanel {
        my ( $geneStart, $geneEnd, $picFormat ) = @_;
        my $panel;
        if ( $picFormat eq 'svg' || $picFormat eq 'SVG' ) {
                $panel = Bio::Graphics::Panel->new(
                        -start       => $geneStart - 100,
                        -stop        => $geneEnd + 100,
                        -width       => 800,
                        -pad_left    => 10,
                        -pad_right   => 10,
                        -key_style   => 'between',
                        -image_class => 'GD::SVG',
                );
        }
        elsif ( $picFormat eq 'png' || $picFormat eq 'PNG' ) {
                $panel = Bio::Graphics::Panel->new(
                        -start     => $geneStart - 100,
                        -stop      => $geneEnd + 100,
                        -width     => 800,
                        -pad_left  => 10,
                        -pad_right => 10,
                        -key_style => 'between',
                );
        }
        else {
                &logError( "Graph format '$picFormat' is not supported!", __FILE__, __LINE__ );
        }
        return $panel;
}

1;
__END__
