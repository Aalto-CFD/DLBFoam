/*

*/

#pragma once

#include "Ostream.H"
#include "Istream.H"

namespace Foam{


///
///@brief A generic chemistryLoad object
///
///
struct chemistryLoad{

    chemistryLoad() = default;

    chemistryLoad(int rank_in, double val_in, int n_active_cells) : 
    rank(rank_in), 
    value(val_in), 
    number_of_active_cells(n_active_cells) 
    {} 


    int rank;     //rank of the process 
    double value; //an indicator how busy the processor is
    int    number_of_active_cells; //number of chemically active cells

    //overload comparisons for sorting
    bool operator==(const chemistryLoad& rhs) const {return value == rhs.value;}
    bool operator!=(const chemistryLoad& rhs) const {return !(value == rhs.value);}
    bool operator<=(const chemistryLoad& rhs) const {return value <= rhs.value;}
    bool operator>=(const chemistryLoad& rhs) const {return value >= rhs.value;}
    bool operator< (const chemistryLoad& rhs) const {return value <  rhs.value;}
    bool operator> (const chemistryLoad& rhs) const {return value >  rhs.value;}
};

Ostream& operator<<(Ostream& os, const chemistryLoad& l);
Istream& operator>>(Istream& is, chemistryLoad& l);





} //namespace Foam