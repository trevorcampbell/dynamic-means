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

#include "constraint.h"
#include "simplex.h"

// Using
using pilal::Matrix;
using pilal::AnonymousMatrix;

namespace optimization {

    /*
        Constraint
        ==========
        Class that represents a constraint: A_i * x_j = b_i.
        
    */
    
    Constraint::Constraint( Matrix const & coefficients, ConstraintType type, long double value ) {
        
        // Coefficients must be a row vector
        if ( coefficients.dim().first == 1 ) {
            this->coefficients = coefficients;
            this->type = type;
            this->value = value;
        } else {
            throw(DataMismatchException("Invalid coefficients vector."));
        }
    }
    
    Constraint::Constraint( Matrix const & coefficients, ConstraintType type, long double lower, long double upper ) {
        
        if ( type != CT_BOUNDS )
            throw(DataMismatchException("Invalid constraint type for provided data"));
            
        // Coefficients must be a row vector
        if ( coefficients.dim().first == 1 ) {
            this->coefficients = coefficients;
            this->type = type;
            this->lower = lower;
            this->upper = upper;
        } else {
            throw(DataMismatchException("Invalid coefficients vector."));
        }
    }
    
    int Constraint::size() const {
        return coefficients.dim().second;
    }
    
    void Constraint::log() const {
        for (int i = 0; i < coefficients.dim().second; ++i)
            std::cout << coefficients(i) << "\t";
        
        switch(type) {

            case CT_EQUAL:
            std::cout << "=\t";
            break;
            
            case CT_LESS_EQUAL:
            std::cout << "<=\t";
            break;
            
            case CT_MORE_EQUAL:
            std::cout << ">=\t";
            break;
            
            case CT_BOUNDS:
            std::cout << "bounded to ";
            break;
            
            case CT_NON_NEGATIVE:
            std::cout << "non-negative ";
        }
        
        if (type == CT_NON_NEGATIVE)
            std::cout << std::endl;
        else if ( type == CT_BOUNDS )
            std::cout << lower << " <= " << "value" << " <= " << upper << std::endl;
        else
            std::cout << value << std::endl;
    }
    
    void Constraint::add_column(long double value) {
        AnonymousMatrix row(1, coefficients.dim().second+1);
        for (int i = 0; i < coefficients.dim().second; ++i)
            row(i) = coefficients(i);
        
        row(coefficients.dim().second) = value;
        coefficients = row;
    }

}
