# This is a GNUPLOT script to plot an orbit in the phase space
eps=1
x11=0

if(!exist("term")) term=x11;

if(term==eps){
set term post eps enhanced color "Symbol" 24
set out "tmp.eps"
set xlabel "{/TimesNewRomanPS-ItalicMT=24 x_n}{/Symbol=24 \256}" offset 1.0,0.0
set ylabel "{/TimesNewRomanPS-ItalicMT=24 x_{n+1}}{/Symbol=24 \256}" offset 1.2,1.0
width=4
}else{
set term x11
if(!exist("width")){width=0.5;}
}

unset key
set size square

PERIOD = 2*pi

mapFrom  = 0
mapTo    = 100
timeFrom = mapFrom * PERIOD
timeTo   = mapTo * PERIOD

p \
"pp.orbit" u 2:($1>timeFrom)?(($1<timeTo)?$3:1/0):1/0 w l lw width lc rgb "black",\
"pp.poin"  u 2:($1>timeFrom)?(($1<timeTo)?$3:1/0):1/0 w p pt 7 lc rgb "red"

if(term==0) pause -1;
#EOF
