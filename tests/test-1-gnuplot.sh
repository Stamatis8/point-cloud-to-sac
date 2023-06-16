set term pdfcairo font "Times-New-Roman,12" size 5.15in, 2in;

set xlabel "Longitudinal Position";
set ylabel "Sectional Area";
#set size ratio -1;
#set key at 0.74,0.205 maxrows 1;
#set key off;
set key outside;

set output "test-1-sac-comparison.pdf";
plot "test-1-analytic.dat" w l title 'analytic',\
     "test-1-approx.dat" w lp title 'approximate'