#ifndef SIMPSON_HPP
#define SIMPSON_HPP

#include <functional>

double simpson(std::function<double(double)> f, std::vector<double> domain, long int n){
    /*
        Approximate the integral of f over domain.at(0) domain.at(1) with n subintervals using
            composite simpson's rule

        Notes: n must be even
    */

    double result = 0;

    double h = (domain.at(1)-domain.at(0))/double(n);

    for(int i = 1; i <= (n+1); i++){
        result += f(domain.at(0)+double(i-1)*h) * (2 + 2*(i%2))/(!(i==1||i==(n+1)) + (i==1||i==(n+1))*(2 + 2*(i%2)));
    }

    result *= h/3.0;

    return result;
}

double simpson(std::vector<double> f, std::vector<double> domain){
    /*
        Approximate the integral of f over domain.at(0) domain.at(1) with n subintervals using
            composite simpson's rule

        Notes: 
            - f.size() must have an odd number of elements greater or equal to 3
            - f.at(i) = f(domain.at(i))
            - f.size() == domain.size()

    */

    int n = f.size()-1;

    double result = 0;

    double h = (domain.back()-domain.at(0))/double(n);

    for(int i = 1; i <= (n+1); i++){
        result += f.at(i-1) * (2 + 2*(i%2))/(!(i==1||i==(n+1)) + (i==1||i==(n+1))*(2 + 2*(i%2)));
    }

    result *= h/3.0;

    return result;
}

double simpson(std::function<double(double,double)> f, std::vector<double> domain, long int n){
    /*
        Approximate the integral of f over [domain.at(0) domain.at(1)] for its first parameter
            and [domain.at(2) domain.at(3)] for its second parameter with n \times n subintervals using
            composite simpson's rule

        Notes: n must be even
    */

    double result = 0;

    double h1 = (domain.at(1)-domain.at(0))/double(n);
    
    double h2 = (domain.at(3)-domain.at(2))/double(n);

    for(int i = 1; i <= (n+1); i++){
        for(int j = 1; j <= (n+1); j++){
            result += f(domain.at(0)+double(i-1)*h1,domain.at(2)+double(j-1)*h2) 
                    * (2 + 2*(i%2))/(!(i==1||i==(n+1)) + (i==1||i==(n+1))*(2 + 2*(i%2)))
                    * (2 + 2*(j%2))/(!(j==1||j==(n+1)) + (j==1||j==(n+1))*(2 + 2*(j%2)));
    
        }
    }

    result *= h1*h2/9.0;

    return result;
}

#endif //SIMPSON_HPP