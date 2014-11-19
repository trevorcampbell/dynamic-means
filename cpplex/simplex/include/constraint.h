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

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "pilal.h"
using pilal::Matrix;

namespace optimization {

    enum ConstraintType {
        
        CT_LESS_EQUAL,
        CT_MORE_EQUAL,
        CT_EQUAL,
        CT_NON_NEGATIVE,
        CT_BOUNDS     
           
    };
        
            
    class Constraint {
        
        friend class Simplex;
        
        public:
                
            Constraint( Matrix const & coefficients, ConstraintType type, long double value );
            Constraint( Matrix const & coefficients, ConstraintType type, long double lower, long double upper );
            
            // Debug
            void log() const; 
            void add_column(long double value);
            int size() const; 
            
        private:
                
            ConstraintType type;
            Matrix coefficients;
            long double value;
            long double upper;
            long double lower;        
            
    };

}

#endif
