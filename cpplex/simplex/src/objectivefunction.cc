/*
Copyright (c) 2010-2013 Tommaso Urli

Tommaso Urli    tommaso.urli@uniud.it   University of Udine

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "objectivefunction.h"
#include "simplex.h"

// Using
using pilal::Matrix;
using pilal::AnonymousMatrix;

namespace optimization {

    /*
        ObjectiveFunction
        =================
        Class that represents an optimization objective function.
        
    */
    
    ObjectiveFunction::ObjectiveFunction() {}
    
    ObjectiveFunction::ObjectiveFunction ( ObjectiveFunctionType type, Matrix const & costs ) :
        type(type),
        costs(costs) {
        
    }
    
    ObjectiveFunction& ObjectiveFunction::operator=( ObjectiveFunction const & objective_function ) {
    
        type = objective_function.type;
        costs = objective_function.costs;
    
        return *this;
    }
    
    void ObjectiveFunction::log() const {
        
        if (type == OFT_MINIMIZE)
            std::cout << "min ( ";
        else
            std::cout << "max ( ";
        
        for (int i = 0; i < costs.dim().second; ++i)
            std::cout << costs(i) << "  ";
        std::cout << ") * x " << std::endl;
        
    }
    
    Matrix const & ObjectiveFunction::get_value( Matrix const & x) const {
        return costs * x;
    }
    
    void ObjectiveFunction::add_column(long double value) {
        AnonymousMatrix row(1, costs.dim().second+1);
        for (int i = 0; i < costs.dim().second; ++i)
            row(i) = costs(i);
        
        row(costs.dim().second) = value;
        costs = row;
    }

}

