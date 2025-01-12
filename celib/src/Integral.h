#pragma once 

double IntegralTrapezoidal(const double x_0, const double x_1, const int Steps, double (*func)(const double));
double IntegralTrapezoidalOneExtraParameter(const double x_0, const double x_1, const int Steps, double (*func)(const double, const double), const double Parameter);
double IntegralSimpson(const double x_0, const double x_1, const int Steps, double (*func)(const double));
double IntegralSimpsonOneExtraParameter(const double x_0, const double x_1, const int Steps, double (*func)(const double, const double), const double Parameter);
double IntegralSimpsonOneExtraParameterType(const double x_0, const double x_1, const int Steps, double (*func)(const double, const int), const int Type);
