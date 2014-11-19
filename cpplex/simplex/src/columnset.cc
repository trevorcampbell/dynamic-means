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

#include "columnset.h"
#include "simplex.h"
#include <algorithm>
#include <vector>

using std::vector;

namespace optimization {
    
    /*
        ColumnSet
        =========
        Class that represents a set of columns.
        
    */
        
    void ColumnSet::insert( int column ) {
        columns.push_back(column);
    }
    
    void ColumnSet::substitute( int old_column, int new_column ) {
        if ( find( columns.begin(), columns.end(), old_column) != columns.end() )
            *(find( columns.begin(), columns.end(), old_column)) = new_column;
    }
    
    void ColumnSet::remove( int column ) {
        if ( find( columns.begin(), columns.end(), column) != columns.end() )
            columns.erase(find( columns.begin(), columns.end(), column));
    }
    
    int ColumnSet::index_of( int column ) {
        int pos = 0;
        for ( vector<int>::iterator it = columns.begin(); it != columns.end(); ++it)
            if ( *it != column )
                ++pos;
            else
                return pos;
        return -1;
    } 
    
    void ColumnSet::log(char const * prelude) const {
         
        
        std::cout << prelude;
        
        for (   vector<int>::const_iterator it = columns.begin(); 
                it != columns.end(); 
                ++it ) {
            std::cout << *it << " ";
        }
        
        std::cout << std::endl;
    }
    
    bool ColumnSet::contains(int column) const {
        return find( columns.begin(), columns.end(), column) != columns.end();
    }
    
    int& ColumnSet::column(int idx) {
        return columns.at(idx);
    }
    
    unsigned int ColumnSet::size() const {
        return columns.size();
    }
}
