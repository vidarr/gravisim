#!/usr/bin/perl -w

=pod

=head1 NAME

normalize.pl [STEADY_PARTICLE_INDEX]

=head1 DESCRIPTION

C<normalize.pl> reades a coordinate vector from stdin.
It will transform all coordinates in a way such that
one of the position of one particle will appear fixed,
i.e. its coordinates be (0, 0).
The particle which will be fixed to a constant position
is chosen by STEADY_PARTICLE_INDEX.
If STEADY_PARTICLE_INDEX is not given, the first particle
will be chosen.

=head1 AUTHORS

Michael J. Beer <michael@ubeer.org>

=cut

my $steady_particle_index = shift || 0;

while (<>) {
    do {
        print $_;
        next;
    } if /^ *#/;
    my @parts = split '  *', $_;
    my $no_particles = ($#parts - 1) / 4;
    my $x = $parts[1 + $steady_particle_index] + 0.0;
    my $y = $parts[2 * $no_particles + 1 + $steady_particle_index] + 0.0;
    @new_parts = ($parts[0]);
    push @new_parts, map {$_ - $x} @parts[1 .. 2 * $no_particles];
    push @new_parts, map {$_ - $y} @parts[2 * $no_particles + 1 .. 4 * $no_particles];
    print join ' ', @new_parts;
    print "\n";
}
