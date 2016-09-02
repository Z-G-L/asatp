#
#
package lib::basic;

=head1 NAME

    A simple log package. All log show in STDERR.

=head1 SYNOPSIS


=head1 AUTHOR

Zhigang Li  lzg0063(at)126.com  2012-8-24

=head1 DESCRIPTION

version 0.2:    
   Change package name from 'warning' to log.

=cut

use strict;

require Exporter;

our @ISA         = qw(Exporter);
our @EXPORT      = qw(logMsg logError logWarn pbar logStd);
our @EXPORT_OK   = qw();
our %EXPORT_TAGS = ();
our $VERSION     = 0.2;

################################################################################

=head2 checkFold

    Title   : checkFold
    Function: 
    Usage   : 
    Returns : 
    Args    : $force: 0 use 'rmdir' delete fold; if 1 use 'rm -rf'
    Des     :
    
=cut

sub checkFold {
    my ( $folder, $force ) = @_;
    if ( -e $folder ) {
        if ( $force == 0 ) {
            rmdir $folder;
            &logStd("\nWARNING:   Output folder '$folder' has existed! Please check it again.\n") if ( -e $folder );
            exit(1);
        }
        else {
            &logStd("\nWARNING:   Output folder '$folder' has existed! It will be deleted!.\n") if ( -e $folder );
            system("rm -rf $folder");
        }
        mkdir( $folder, 0755 ) || die;
    }
    else {
        mkdir( $folder, 0755 ) || die;
    }

}

################################################################################

=head2 logStd

    Title   : logStd
    Function: output in STDERR
    Usage   : &logStd($string)
    Returns : In STDERR
    Args    : a string

=cut

sub logStd {
    my ($string) = @_;
    print STDERR "$string\n";
}

=head2 logMsg

    Title   : logMsg
    Function: output warning in STDERR
    Usage   : &logMsg($string)
    Returns : In STDERR
    Args    : a string
    
=cut

sub logMsg {
    my ($string) = @_;
    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime();
    my $format_time = sprintf( "%d-%d-%d %d:%d:%d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
    print STDERR "\n", $string, "\t\t", $format_time, "\n";
}

=head2 logError

    Title   : logError
    Function: output Error and exit
    Usage   : &logError($string,__FILE__,__LINE__)
    Returns : In STDERR
    Args    : a string
    
=cut

sub logError {
    my ( $string, $fileName, $lineNum ) = @_;
    print STDERR "!! [ ERROR ] $fileName\tline: $lineNum\n";
    print STDERR "!! $string\n";
    print STDERR "!! ** Exit 1 **\n";
    exit 1;
}

=head2 logWarn

    Title   : logWarn
    Function: output warning
    Usage   : %logWarn($string,__FILE__,__LINE__)
    Returns : In STDERR
    Args    : a string
    
=cut

sub logWarn {
    my ( $string, $fileName, $lineNum ) = @_;
    print STDERR "!! [ WARNING ] $string $fileName\tline: $lineNum\n";
}

1;

__END__
