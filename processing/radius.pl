#!/usr/bin/perl -w

=pod

=head1 NAME

perl radius.pl

=head1 DESCRIPTION

C<radius.pl>

=head1 AUTHORS

Michael J. Beer <michael@ubeer.org>

=cut

while (<>) {
    do {
        print $_;
        next;
    } if /^ *#/;
    my @parts = split '  *', $_;
    my $no_particles = $#parts / 4;
    my ($x, $y);
    my @radiuses = ();
    for(my $i = 0; $i < $no_particles; $i++) {
        $x = $parts[1 + $i] + 0.0;
        $y = $parts[1 + 2 * $no_particles + $i] + 0.0;
        push @radiuses, sqrt($x * $x + $y *$y);
    }
    print $parts[0] . ' ';
    print join ' ', @radiuses;
    print "\n";
}
