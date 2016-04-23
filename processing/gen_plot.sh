#!/usr/bin/bash
: <<EOD

=pod

=head1 NAME

   bash gen_plot.sh DATA_FILE TITLE START_INDEX END_INDEX

=head1 SYNOPSIS

=cut

EOD

# OUT_CMD=cat
OUT_CMD=gnuplot

if [ $# -lt 3 ]; then
    perldoc $0;
    exit 1;
fi


function generate_plot_line_cmd()
{
    COLUMN=$(( $INDEX + 1))
    PLOT_LINE_CMD="'$DATA_FILE' u 1:${COLUMN}"
}

DATA_FILE=$1
TITLE=$2
START=$3
END=$4

PLOT_CMD=""
INDEX=$START
START=$(($START + 1))
generate_plot_line_cmd
PLOT_CMD=$PLOT_LINE_CMD
for INDEX in `seq $START $END`; do
    generate_plot_line_cmd
    PLOT_CMD="$PLOT_CMD ,$PLOT_LINE_CMD"
done


$OUT_CMD <<EOC
set t png size 6000,6000
set o "$TITLE.png"
set title "$TITLE"
plot $PLOT_CMD
exit
EOC
