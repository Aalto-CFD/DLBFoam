#!/bin/bash
# Plot chemistry solve buffer times per rank using gnuplot

OUTFILE="rankbased_solve.png"
TMPDIR=$(mktemp -d)
DATAFILE="$TMPDIR/all_ranks.dat"
MEANFILE="$TMPDIR/mean.dat"

# Collect processor directories (sorted)
PROC_DIRS=($(ls -d processor* 2>/dev/null | sort))
NPROCS=${#PROC_DIRS[@]}

if [ "$NPROCS" -eq 0 ]; then
    echo "No processor* directories found."
    exit 1
fi

# Determine number of time steps from processor0 (skip header, drop last 2 rows)
NLINES=$(tail -n +2 "${PROC_DIRS[0]}/loadBal/cpu_solve.out" | wc -l)
SIZE=$(( NLINES - 2 ))

if [ "$SIZE" -le 0 ]; then
    echo "Not enough data in cpu_solve.out (need > 2 lines after header)."
    exit 1
fi

# Extract solve_buffer column (column 5) for each rank into separate tmp files
RANK_FILES=()
for rank in "${PROC_DIRS[@]}"; do
    rfile="$TMPDIR/${rank}.dat"
    tail -n +2 "${rank}/loadBal/cpu_solve.out" | head -n "$SIZE" | awk '{print NR, $5}' > "$rfile"
    RANK_FILES+=("$rfile")
done

# Compute mean solve_buffer across all ranks
awk -v nprocs="$NPROCS" '
{
    sum[$1] += $2
    count[$1]++
}
END {
    for (i = 1; i <= length(sum); i++)
        print i, sum[i] / nprocs
}
' "${RANK_FILES[@]}" | sort -k1,1n > "$MEANFILE"

# Build gnuplot script
GPSCRIPT="$TMPDIR/plot.gp"

# Start plot commands: one line per rank
{
    echo "set terminal pngcairo size 1200,600 enhanced font 'Arial,10'"
    echo "set output '$OUTFILE'"
    echo "set xlabel 'Number of iterations'"
    echo "set ylabel 'Chemistry CPU time [s]'"
    echo "set key top right box off"
    echo "set xrange [0:*]"
    echo "set yrange [0:*]"
    echo "set border 15"
    echo "set tics in"
    echo "set mxtics"
    echo "set mytics"
    echo ""

    # Color palette for ranks (cycles automatically)
    COLORS=('#e41a1c' '#377eb8' '#4daf4a' '#984ea3' '#ff7f00' '#a65628' '#f781bf' '#999999')

    # Build plot command
    PLOT_CMD="plot "
    FIRST=1
    for i in "${!RANK_FILES[@]}"; do
        rfile="${RANK_FILES[$i]}"
        color="${COLORS[$((i % ${#COLORS[@]}))]}"
        if [ "$FIRST" -eq 1 ]; then
            PLOT_CMD+="'$rfile' using 1:2 with lines lw 1.5 lc rgb '$color' title 'Rank $i', \\"
            FIRST=0
        else
            PLOT_CMD+="'$rfile' using 1:2 with lines lw 1.5 lc rgb '$color' title 'Rank $i', \\"
        fi
        echo "$PLOT_CMD"
        PLOT_CMD=""
    done
    echo "'$MEANFILE' using 1:2 with lines lw 2.5 lc rgb 'black' dt 2 title 'Mean'"

} > "$GPSCRIPT"

gnuplot "$GPSCRIPT"

if [ $? -eq 0 ]; then
    echo "Plot saved to $OUTFILE"
else
    echo "gnuplot failed."
    exit 1
fi

rm -rf "$TMPDIR"
