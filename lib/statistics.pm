#
#
package lib::statistics;

=head1 NAME



=head1 SYNOPSIS



=head1 AUTHOR

Zhigang Li	lzg0063(at)126.com	2013-10-24

=head1 DESCRIPTION



=cut

use 5.010;
use strict;
#use warnings;
use lib::basic qw(logMsg logError logWarn pbar);
use Statistics::R;

require Exporter;

our @ISA         = qw(Exporter);
our @EXPORT      = qw();
our @EXPORT_OK   = qw();
our %EXPORT_TAGS = ();
our $VERSION     = 0.1;







sub ChiSquareTest {
    my ( $data, $store ) = @_;

    my $R = Statistics::R->new();
    $R->set( 'x', $data );

    my $cmds = <<"EOF";
    p=c()
    for( i in c(1:(length(x)/4)) )
    {
        s = (i-1) * 4 + 1
        p[i] = chisq.test(matrix(c(x[s],x[s+1],x[s+2],x[s+3]),nrow=2),correct=T)\$p.value
    }
    print(p)
    q = p.adjust(p, method="BH")
EOF
    $R->run($cmds);
    $store->{'pvalue'} = $R->get('p');
    $store->{'qvalue'} = $R->get('q');
    $R->stop();

    return 0;
}


1;
__END__