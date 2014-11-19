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

#ifndef VARIABLE_H
#define VARIABLE_H

#include "constraint.h"          
#include "simplex.h"

namespace optimization {

    class AuxiliaryVariable;

    class Variable {
        
        friend class Simplex;
        
        public:
            Variable( Simplex * creator, char const * name); 
            virtual ~Variable();
            virtual void process(Matrix& calculated_solution, Matrix& solution, int index);
            
        protected:
            std::string name;   
            Simplex * creator;         
            
    };
    
    class SplittedVariable : public Variable {
    
        friend class Simplex; 
        
        public:
            SplittedVariable( Simplex* creator, char const * name, AuxiliaryVariable* aux);
            ~SplittedVariable();
            void process(Matrix& calculated_solution, Matrix& solution, int index); 
        private:
            AuxiliaryVariable* aux;
    
    };
    
    class SlackVariable : public Variable {
    
        friend class Simplex;
    
        public:
            SlackVariable(Simplex * creator, char const * name);
            ~SlackVariable();
            void process(Matrix& calculated_solution, Matrix& solution, int index); 

        
    };
    
    class AuxiliaryVariable : public Variable {
        
        friend class Simplex;
        friend class SplittedVariable;
        
        public:
            AuxiliaryVariable(Simplex* creator, char const * name, int index);
            ~AuxiliaryVariable();   
            void process(Matrix& calculated_solution, Matrix& solution, int index); 
        private:
            int index;
       
    };    
}

#endif
