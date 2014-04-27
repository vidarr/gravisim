#!/usr/bin/bash
: <<EOD

=pod

=head1 NAME

   orbit_10.sh

=head1 SYNOPSIS

   orbit_10.sh DATA_FILE OUTPUT_FILE TIME_END [TITLE]

=cut

EOD

if [ $# -lt 3 ]; then
    perldoc $0;
    exit 1;
fi

DATA_FILE=$1
TITLE=$4

gnuplot <<EOC
set t png size 2000,2000
set o "$2"
set title "$TITLE time_end = $3"
plot '$DATA_FILE' u 2:22, '$DATA_FILE' u 3:23 w l, '$DATA_FILE' u 4:24 w l, '$DATA_FILE' u 5:25 w l, '$DATA_FILE' u 6:26 w l, '$DATA_FILE' u 7:27 w l, '$DATA_FILE' u 8:28 w l, '$DATA_FILE' u 9:29 w l, '$DATA_FILE' u 10:30 w l, '$DATA_FILE' u 11:31 w l
exit
EOC
