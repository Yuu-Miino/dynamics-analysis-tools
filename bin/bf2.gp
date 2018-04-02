# GNUPLOT script to plot an orbit in a phase space
eps=1
x11=0

# Terminal options
if(!exist("term")) term=x11;
if(term==eps){
set term post eps enhanced color "Symbol" 24
set out "tmp.eps"
set xlabel "{/TimesNewRomanPS-ItalicMT=24 B}{/Symbol=24 \256}" offset 1.0,0.0
set ylabel "{/TimesNewRomanPS-ItalicMT=24 B_{/Symbol= 0}}{/Symbol=24 \256}" offset 1.2,1.0
width=4
}else{
set term x11
if(!exist("width")){width=0.5;}
}

# Input file
filesI="\
07.pt.1-2(+1.0e-03).bf2.cont \
07.pt.1-2(-1.0e-03).bf2.cont \
08.pt.2-1(-5.0e-04).bf2.cont \
"
fileNew="\
08.pt.2-1(-5.0e-04).bf2.cont \
"


# Plotting options
unset key
set size square

#set xr [-2:2]
#set yr [-2:2]

p \
for [infile in filesI] infile u 3:2 w l lw width lc rgb "blue",\
for [infile in fileNew] infile u 3:2 w l lw width lc rgb "cyan",\


if(term==0) pause -1;
#EOF
