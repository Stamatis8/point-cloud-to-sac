set term pdfcairo font "Times-New-Roman,12" size 5.15in, 2in;

set xlabel "Longitudinal Position";
set ylabel "Sectional Area of KCSsim hull";
#set size ratio -1;
#set key at 0.74,0.205 maxrows 1;

set output "test-3-sac-random-uniform.pdf";
set key outside;
plot "test-3-approx-random.dat" w lp pt 0 title 'random',\
    "test-3-approx-uniform.dat" w lp pt 0 title 'uniform'

set output "test-3-sac.pdf";
plot "test-3-approx.dat" w lp pt 0 title ''

