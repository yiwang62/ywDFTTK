
eV = 96484.0
kJ = 1000.0

G(x) = a + b*x +c*x*log(x) + d*x*x + e*x*x*x + f/x + g*x**7 + h/x**9
S(x) = -b -c -c*log(x) - 2.0*d*x -3.0*e*x*x + f/x/x - 7.0*g*x**6 - 9.0*h/x**10
H(x) = G(x) + x*S(x)
Cp(x) = -c - 2.0*d*x -6.0*e*x*x - 2.0*f/x/x - 42.0*g*x**6 - 90.0*h/x**10

set key right top
set xlabel "T(K)"
set ylabel "Gibbs energy (eV/atom)"
plot [x=x_RANGE] "x_FILE" using 1:3 title "Calculated" w p pt 5, \
	"" using 1:(G($1)/eV) title "Fitted" w l lt 3

pause -1 

set key left top
set xlabel "T(K)"
set ylabel "Entropy (J/K/mol-atom)"
plot [x=x_RANGE] "x_FILE" using 1:4 title "Calculated" w p pt 5, \
        "" using 1:(S($1)) title "Fitted" w l lt 3

pause -1 

set key left top
set xlabel "T(K)"
set ylabel "Ethalpy (J/mol-atom)"
plot [x=x_RANGE] "x_FILE" using 1:5 title "Calculated" w p pt 5, \
        "" using 1:(H($1)) title "Fitted" w l lt 3

pause -1 

set key left top
set xlabel "T(K)"
set ylabel "C_p (J/K/mol-atom)"
plot [x=x_RANGE] "x_FILE" using 1:7 title "Calculated" w p pt 5, \
        "" using 1:(Cp($1)) title "Fitted" w l lt 3

pause -1 
