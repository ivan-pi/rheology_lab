set term epscairo
set output 'tok_gel4.eps'
set xrange [0:110]
set title '0,5 % Gellan, 2. interval'

set xlabel 'Strizna hitrost [1/s]'
set ylabel 'Strizna napetost [Pa]'

k = 0.00000105080913703538
vr0 = 0.000000000396648477393813
n = 0.0396249008723081
f(x) = (vr0**n + (k*x)**n)**(1/n)

set key bottom

plot 'list.dat' u 1:2 title 'Meritve', f(x) title 'Model'
reset
set term wxt
