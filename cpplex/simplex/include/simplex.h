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

#ifndef SIMPLEX_H
#define SIMPLEX_H

// PILAL
#include "pilal.h"

// Simplex classes
#include "simplexexceptions.h"
#include "constraint.h"
#include "objectivefunction.h"
#include "columnset.h"


// From the STL
#include <iostream>
#include <vector>

// Using
using pilal::Matrix;
using pilal::AnonymousMatrix;

namespace optimization {
    
    class ColumnSet;
    class Constraint;
    class ObjectiveFunction;
    class Variable;

    class Simplex {
        
        public:  
        
            // Constructor
            Simplex(char const * name);
            ~Simplex();
            
            // Settings                                
            void load_problem(char const * problem_name);
            void add_variable(Variable* variable);
            void add_constraint(Constraint const & constraint);
            void set_objective_function(ObjectiveFunction const & objective_function);
            
            // Solving procedures
            void solve();          

            // Print
            void print_solution() const;
            void log() const; 
            
            bool is_unlimited() const;
            bool has_solutions() const;
            bool must_be_fixed() const;
            Matrix const & get_dual_variables() const;
                                  
            
        protected:
                                                             
            std::string name;
            
            // Preprocessing
            void process_to_standard_form();
            void process_to_artificial_problem();
         
            // Solving
            void solve_with_base( ColumnSet const& base );
            
            // Column sets
            ColumnSet suggested_base;
            ColumnSet current_base;
            ColumnSet current_out_of_base;
            
            // Data
            ObjectiveFunction objective_function;
            std::vector< Constraint> constraints;
            std::vector< Constraint> nn_constraints;            
            std::vector< Variable*> variables;
            
            // Processed data
            Matrix costs;
            Matrix coefficients_matrix;
            Matrix constraints_vector;
            Matrix base_inverse;   
            Matrix dual_variables;
            Matrix column_p;
            int solution_dimension, old_column;
            
            // Results
            Matrix base_solution;  
            Matrix solution;
            Matrix reduced_cost;
            long double solution_value;
            
            bool    optimal, 
                    unlimited, 
                    overconstrained, 
                    has_to_be_fixed, 
                    changed_sign;
            
            int inverse_recalculation_rate;
          
    };      

}

#endif
