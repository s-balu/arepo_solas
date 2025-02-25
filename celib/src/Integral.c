#include "config.h"

/*!
 * This function numerically integrates the given function (the integrand) using the trapezoidal rule.
 */
double IntegralTrapezoidal(const double x_0, const double x_1, const int Steps, double (*func)(const double)){

    double h = (x_1-x_0)/(double)Steps;
    double Sum = 0.0;
    for(int i=1;i<Steps;i++){
        double x = i*h+x_0;
        Sum += func(x);
    }

    Sum += 0.5*(func(x_0)+func(x_1));

    return h*Sum;
}

/*!
 * This function numerically integrates the given function (the integrand) using
 * the trapezoidal rule.  This function accepts one extra parameter which is
 * necessary to evaluate the integrand.
 */
double IntegralTrapezoidalOneExtraParameter(const double x_0, const double x_1, const int Steps, double (*func)(const double, const double), const double Parameter){

    double h = (x_1-x_0)/(double)Steps;
    double Sum = 0.0;
    for(int i=1;i<Steps;i++){
        double x = i*h+x_0;
        Sum += func(x,Parameter);
    }

    Sum += 0.5*(func(x_0,Parameter)+func(x_1,Parameter));

    return h*Sum;
}

/*!
 * This function numerically integrates the given function (the integrand) using the Simpson rule.
 */
double IntegralSimpson(const double x_0, const double x_1, const int Steps, double (*func)(const double)){

    double h = (x_1-x_0)/(2.0*Steps);
    double Sum_odd = 0.0;
    double Sum_even = 0.0;
    /*
    for(int i=1;i<2*Steps-3;i+=2){
        Sum_odd += func(x_0+h*i);
        Sum_even += func(x_0+h*(i+1));
    }
    double Sum = (func(x_0)+func(x_1)+4*(Sum_odd+func(x_1-h))+2*Sum_even)*h/3.0;
    */
    for(int i=1;i<2*Steps;i+=2){
        Sum_odd += func(x_0+h*i);
    }
    for(int i=2;i<2*Steps;i+=2){
        Sum_even += func(x_0+h*i);
    }
    double Sum = (func(x_0)+func(x_1)+4*Sum_odd+2*Sum_even)*h/3.0;

    return Sum;
}

/*!
 * This function numerically integrates the given function (the integrand) using
 * the Simpson rule.  This function accepts one extra parameter which is
 * necessary to evaluate the integrand.
 */
double IntegralSimpsonOneExtraParameter(const double x_0, const double x_1, const int Steps, double (*func)(const double, const double), const double Parameter){

    double h = (x_1-x_0)/(2.0*Steps);
    double Sum_odd = 0.0;
    double Sum_even = 0.0;

    for(int i=1;i<2*Steps;i+=2){
        Sum_odd += func(x_0+h*i,Parameter);
    }
    for(int i=2;i<2*Steps;i+=2){
        Sum_even += func(x_0+h*i,Parameter);
    }
    double Sum = (func(x_0,Parameter)+func(x_1,Parameter)+4*Sum_odd+2*Sum_even)*h/3.0;

    return Sum;
}

/*!
 * This function numerically integrates the given function (the integrand) using
 * the Simpson rule.  This function accepts one extra control parameter which
 * can control the state of the integrand.
 */
double IntegralSimpsonOneExtraParameterType(const double x_0, const double x_1, const int Steps, double (*func)(const double, const int), const int Type){

    double h = (x_1-x_0)/(2.0*Steps);
    double Sum_odd = 0.0;
    double Sum_even = 0.0;

    for(int i=1;i<2*Steps;i+=2){
        Sum_odd += func(x_0+h*i,Type);
    }
    for(int i=2;i<2*Steps;i+=2){
        Sum_even += func(x_0+h*i,Type);
    }
    double Sum = (func(x_0,Type)+func(x_1,Type)+4*Sum_odd+2*Sum_even)*h/3.0;

    return Sum;
}
