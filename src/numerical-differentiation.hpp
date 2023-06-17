#ifndef NUMERICAL_DIFFERENTIATION_HPP
#define NUMERICAL_DIFFERENTIATION_HPP

double Fwd1pointNumDiff(double f0, double f1, double h){
    /*
        Description: evaluate 1 point forward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (f1-f0)/h;
}

double Fwd2pointNumDiff(double f0, double f1, double f2, double h){
    /*
        Description: evaluate 2 point forward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (-3.0*f0 + 4.0*f1 - f2)/(2.0*h);
}

double Fwd3pointNumDiff(double f0, double f1, double f2, double f3, double h){
    /*
        Description: evaluate 3 point forward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (-11.0*f0+18.0*f1-9.0*f2+2.0*f3)/(6.0*h);
}

double Fwd4pointNumDiff(double f0, double f1, double f2, double f3, double f4, double h){
    /*
        Description: evaluate 4 point forward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (-25.0*f0 + 48.0*f1 -36.0*f2 + 16.0*f3 - 3.0*f4)/(12.0*h);
}

double Fwd5pointNumDiff(double f0, double f1, double f2, double f3, double f4, double f5, double h){
    /*
        Description: evaluate 5 point forward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (-137.0*f0 + 300.0*f1 -300.0*f2 + 200.0*f3 - 75.0*f4 + 12.0*f5)/(60.0*h);
}

double Bck1pointNumDiff(double f0, double fm1, double h){
    /*
        Description: evaluate 1 point backward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (f0-fm1)/h;
}

double Bck2pointNumDiff(double f0, double fm1, double fm2, double h){
    /*
        Description: evaluate 2 point backward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (3.0*f0 - 4.0*fm1 + fm2)/(2.0*h);
}

double Bck3pointNumDiff(double f0, double fm1, double fm2, double fm3, double h){
    /*
        Description: evaluate 3 point backward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (11.0*f0-18.0*fm1+9.0*fm2-2.0*fm3)/(6.0*h);
}

double Bck4pointNumDiff(double f0, double fm1, double fm2, double fm3, double fm4, double h){
    /*
        Description: evaluate 4 point backward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (25.0*f0 - 48.0*fm1 + 36.0*fm2 - 16.0*fm3 + 3.0*fm4)/(12.0*h);
}

double Bck5pointNumDiff(double f0, double fm1, double fm2, double fm3, double fm4, double fm5, double h){
    /*
        Description: evaluate 5 point backward numerical derivative
        
        assume h large enough to be a divisor
    */

    return (137.0*f0 - 300.0*fm1 + 300.0*fm2 - 200.0*fm3 + 75.0*fm4 - 12.0*fm5)/(60.0*h);
}

#endif //NUMERICAL_DIFFERENTIATION_HPP