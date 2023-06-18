set term pdfcairo font "Times-New-Roman,12" size 5.15in, 2in;

set xlabel "Longitudinal Position";
set ylabel "Sectional Area";
#set size ratio -1;
#set key at 0.74,0.205 maxrows 1;
#set key off;
set key outside;

set output "test-4-sac-comparison.pdf";
plot "test-4-approx.dat" w l title '',\
     "test-4-approx-deriv-1-aft.dat" w l title '1-point tangent aft',\
     "test-4-approx-deriv-1-fwd.dat" w l title '1-point tangent fwd'