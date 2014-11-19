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

#include "variable.h"

namespace optimization {

    /*
        Variable
        ========
        Base variable class.
        
    */

    Variable::Variable( Simplex* creator,  char const * name) : 
        name(name),
        creator(creator) { 
    }
    
    Variable::~Variable() { }             
     
    void Variable::process(Matrix& calculated_solution, Matrix& solution, int index) {
        solution(index) = calculated_solution(index);
    }
    
    /*
        SplittedVariable
        ================
        Variables that are splitted in two (one) AuxiliaryVariables during
        the translation in standard form.
    
    */
    
    SplittedVariable::SplittedVariable( Simplex* creator, char const * name, AuxiliaryVariable* aux ) :
        Variable(creator, name),
        aux(aux) {
    }           
    
    SplittedVariable::~SplittedVariable() { }

    void SplittedVariable::process(Matrix& calculated_solution, Matrix& solution, int index) { 
        solution(index) = calculated_solution(index) - calculated_solution(aux->index);
    }

    /*
        SlackVariable
        =============
        Type of variable added when transforming a <= or >= constraint
        into a = constraint.
    
    */
    
    SlackVariable::SlackVariable( Simplex* creator, char const * name ) :
        Variable(creator, name) {
    }
   
    SlackVariable::~SlackVariable() { }

    void SlackVariable::process(Matrix& calculated_solution, Matrix& solution, int index) {}
   
    /*
        AuxiliaryVariable
        =================
        Variable created when transforming a variable in a splitted
        variable. The relation:
        
            x = x+ - x-
            
        holds between the original variable, the SplittedVariable and
        the AuxiliaryVariable.
    
    */
    
    AuxiliaryVariable::AuxiliaryVariable(  Simplex* creator, char const * name, int index ) :
        Variable(creator, name),
        index(index) {
    
    }
    
    AuxiliaryVariable::~AuxiliaryVariable() { }
    
    void AuxiliaryVariable::process(Matrix& calculated_solution, Matrix& solution, int index) {}
}

