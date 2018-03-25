# GNUPLOT script to plot an orbit in a phase space
eps=1
x11=0

# Terminal options
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

# Input file
if(!exist("file")) file="01.pt";

# Plotting options
unset key
set size square

set xr [-2:2]
set yr [-2:2]

PERIOD = 2*pi
mapFrom  = 0
mapTo    = 100
timeFrom = mapFrom * PERIOD
timeTo   = mapTo * PERIOD

p \
file.".pp.orbit" u 2:($1>timeFrom)?(($1<timeTo)?(($4==0||$4==1)?$3:1/0):1/0):1/0 w l lw width lc rgb "red",\
file.".pp.orbit" u 2:($1>timeFrom)?(($1<timeTo)?(($4==2||$4==3)?$3:1/0):1/0):1/0 w l lw width lc rgb "blue",\
file.".pp.poin"  u 2:($1>timeFrom)?(($1<timeTo)?$3:1/0):1/0 w p pt 7 lc rgb "black"

if(term==0) pause -1;
#EOF
