#!/usr/bin/perl -w
while (<>) {
    do {
        print $_;
        next;
    } if /^ *#/;
    my @parts = split '  *', $_;
    my $no_particles = ($#parts - 1) / 4;
    my $x = $parts[1] + 0.0;
    my $y = $parts[2 * $no_particles + 1] + 0.0;
    @new_parts = ($parts[0]);
    push @new_parts, map {$_ - $x} @parts[1..2 * $no_particles];
    push @new_parts, map {$_ - $y} @parts[2 * $no_particles + 1..3 * $no_particles];
    push @new_parts, @parts[2 * $no_particles + 1 .. $#parts];
    print join ' ', @new_parts;
    print "\n";
}
