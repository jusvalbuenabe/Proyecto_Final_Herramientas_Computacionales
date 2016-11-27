reset

	set terminal epslatex size 3.5,2.62 standalone color colortext
	set output 'plot.tex'

	set autoscale
	set lmargin 6
	set rmargin 3
	set xlabel '\small $T/T_{D}$'
	set ylabel '$c_{N}/k_{B}$' offset 3.5
	set format '\scriptsize $%.1f$'

	set key at 2.8,2.3
    set key spacing 0.8
    set key samplen 2
    set key box 3
    set key width -4
    set key height 0.3
    set key reverse Left
    
    set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
	set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
    show grid
    
    set yr [0:3.2]
    #set xtics (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
	plot	sprintf('specifHeat.dat') using 1:2 title '\scriptsize Debye' with line lw 3, \
			sprintf('specifHeat.dat') using 1:3 title '\scriptsize $T_{E}=0.5\ T_{D}$' with line lw 3, \
			sprintf('specifHeat.dat') using 1:4 title '\scriptsize $T_{E}=0.75\ T_{D}$' with line lw 4 lc rgb "#f0e442" dashtype 12, \
			sprintf('specifHeat.dat') using 1:5 title '\scriptsize $T_{E}=1.0\ T_{D}$' with line lw 3, \
			sprintf('specifHeat.dat') using 1:6 title '\scriptsize $T_{E}=1.5\ T_{D}$' with line lw 3 lc rgb "#dc143c"
			

