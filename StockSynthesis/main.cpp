/* 
 * File:   main.cpp
 * Author: matthewsupernaw
 *
 * Created on May 27, 2014, 1:29 PM
 */

#include <cstdlib>
#include <iostream>
#include "test/CatchAtAge.hpp"
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    ss::CatchAtAgeData<double> data;

    ss::CatchAtAge<double, double> ca;
    ca.SetPhase(1);
    ca.SetData(data);

    ss::MortalityAndSurvivability<double, double> mortality;
    ca.AddFunctor(&mortality);

    ss::NumbersAtAgeFunctor<double, double> numbers;
    ca.AddFunctor(&numbers);

    ss::CatchFunctor<double, double> catch_functor;
    ca.AddFunctor(&catch_functor);

    ss::CostFunctor<double, double> cost;
    ca.AddFunctor(&cost);

    std::cout << ca.Evaluate();
    return 0;
}

