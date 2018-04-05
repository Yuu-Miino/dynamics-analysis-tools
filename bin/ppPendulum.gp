# GNUPLOT script to plot an orbit in a phase space
eps=1
x11=0

# Terminal options
if(!exist("term")) term=x11;
if(term==eps){
set term post eps enhanced color "Symbol" 24
set out "tmp.eps"
set xlabel "{/TimesNewRomanPS-ItalicMT=24 x}{/Symbol=24 \256}" offset 1.0,0.0
set ylabel "{/TimesNewRomanPS-ItalicMT=24 y}{/Symbol=24 \256}" offset 1.2,1.0
width=4
}else{
set term x11
if(!exist("width")){width=0.5;}
}

# Input file
if(!exist("file")) file="07.pt.bf2.pt";

# Plotting options
unset key
set size square

set xr [-pi-0.2:pi+0.2]
set yr [-2:2]
mod(x,y) = (x-floor(x/y)*y)

p \
file.".pp.orbit" u (mod($2+pi,2*pi)-pi):3 every 10 w l lw 0.5 lc rgb "black",\
file.".pp.poin"  u (mod($2+pi,2*pi)-pi):3 w p pt 7 lc rgb "red"

if(term==0) pause -1;
#EOF
