#set term epscairo
#set output 'cross_gel4.eps'
set xrange [0.009:110]
set title '0,5 % Gellan, 2. interval'

set xlabel 'Strizna hitrost [1/s]'
set ylabel 'Viskoznost [Pa s]'

k = 0.21212564925826
vr0 = 1.27708990806004
vrinf = 0
n = 0.661890668379998
f(x) = (vr0-vrinf)/(1+(k*x)**n)+vrinf

tau0 =4*10e-10
vinf =1.03e-6
m = 0.0396

g(x) = ((tau0)**m+(vinf*x)**m)**(1.0/m)/x

set key top
#set logscale xy

plot 'list.dat' u 1:2 title 'Meritve', f(x) title 'Cross', g(x) title 'Potencna'
reset
#set term wxt