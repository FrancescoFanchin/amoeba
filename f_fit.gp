set xlabel "tempo (s)"
set ylabel "V / # cellule"
set xrange[0:80000]

set title " E=0.41"
plot  "result.dat" using 1:5 t 'f_{teo}' w l linestyle 1 linecolor 1, "f_exp" t "f_{exp}" pt 2
pause -1