#!/usr/bin/perl -w

=pod

=head1 NAME

filter.pl

=head1 SYNOPSIS

perl filter.pl START_TIME [END_TIME]

=cut

sub usage {
    perldoc $0;
}

my $line = 0;
my $time_start = shift or die usage();
my $time_end   = shift;
while(<>) {
    /^#/ && next;
    / *([^ ]*) / && do
            {
                if($time_start < $1) {
                    next if(defined($time_end) && $time_end < $1);
                    print $_;
                };
            };
}

