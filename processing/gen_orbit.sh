#!/usr/bin/bash
: <<EOD

=pod

=head1 NAME

   bash gen_orbit.sh DATA_FILE TITLE TOTAL_NUMBER OF PARTICLES [NUMBER OF PARTICLE TO PLOT]

=head1 SYNOPSIS

   bash gen_orbit.sh 10.csv particle_1 1009  10 > orbit-10.plt
   If NUMBER_OF_PARTICLES_TO_PLOT is omitted, all particles will be plotted

=cut

EOD

#OUT_CMD=cat
OUT_CMD=gnuplot

if [ $# -lt 3 ]; then
    perldoc $0;
    exit 1;
fi

function generate_plot_line_cmd()
{
    FIRST_COLUMN=$((1 + $PARTICLE_NO_TO_PLOT ))
    SECOND_COLUMN=$(( $TOTAL_NO * 2 ))
    SECOND_COLUMN=$(( $SECOND_COLUMN + $FIRST_COLUMN ))
    PLOT_LINE_CMD="'$DATA_FILE' u ${FIRST_COLUMN}:${SECOND_COLUMN}"
}

DATA_FILE=$1
TITLE=$2
TOTAL_NO=$3
PARTICLE_NO_TO_PLOT=$4

PLOT_CMD=""
if [ -z "$PARTICLE_NO_TO_PLOT" ]; then
    PARTICLE_NO_TO_PLOT=0
    generate_plot_line_cmd
    PLOT_CMD=$PLOT_LINE_CMD
    for PARTICLE_NO_TO_PLOT in `seq $TOTAL_NO`; do
        generate_plot_line_cmd
        PLOT_CMD="$PLOT_CMD ,$PLOT_LINE_CMD"
    done
else
    echo "Particle no $PARTICLE_NO_TO_PLOT out of $TOTAL_NO - first col no $FIRST_COLUMN , second col no $SECOND_COLUMN"
    generate_plot_line_cmd
    PLOT_CMD=$PLOT_LINE_CMD
fi


$OUT_CMD <<EOC
set t png size 6000,6000
# set xrange [-2000000000:2000000000]
# set yrange [-40000000:40000000]
set o "$TITLE.png"
set title "$TITLE"
plot $PLOT_CMD
exit
EOC
