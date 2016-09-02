#
#
package lib::math;

=head1 NAME



=head1 SYNOPSIS



=head1 AUTHOR

Zhigang Li  lzg0063(at)126.com  2015-10-24

=head1 DESCRIPTION



=cut

use 5.010;
use strict;

#use warnings;
use lib::basic qw(logMsg logError logWarn pbar);

require Exporter;

our @ISA         = qw(Exporter);
our @EXPORT      = qw();
our @EXPORT_OK   = qw();
our %EXPORT_TAGS = ();
our $VERSION     = 0.1;



sub bin2dec {
    my ($bin) = @_;
    # larger than cutoff, use function to convert
    if ( length($bin) < 31 ) {
        return oct("0b$bin");
    }
    else {
        use Math::BigInt;
        return Math::BigInt->new("0b$bin");
    }
}

sub dec2bin {
    my $dec       = shift;

        # larger than cutoff, use function to convert
    if ( $dec < 4294967290 ) {
        return sprintf( "%b", $dec );
    }
    else {
        use Math::BigInt;
        my $i = Math::BigInt->new("$dec");
        return substr( $i->as_bin(), 2 );
    }
}

1;
__END__
