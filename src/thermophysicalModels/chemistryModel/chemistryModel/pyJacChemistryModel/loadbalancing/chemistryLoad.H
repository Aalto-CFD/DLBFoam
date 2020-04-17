/*

*/

#pragma once

namespace Foam{


///
///@brief A generic chemistryLoad object
///
///
struct chemistryLoad{

    chemistryLoad() = default;

    chemistryLoad(int rank_in, double val_in) : rank(rank_in), value(val_in) {} 


    int rank;     //rank of the process 
    double value; //an indicator how busy the processor is


    //overchemistryLoad comparisons for sorting
    bool operator==(const chemistryLoad& rhs) const {return value == rhs.value;}
    bool operator!=(const chemistryLoad& rhs) const {return !(value == rhs.value);}
    bool operator<=(const chemistryLoad& rhs) const {return value <= rhs.value;}
    bool operator>=(const chemistryLoad& rhs) const {return value >= rhs.value;}
    bool operator< (const chemistryLoad& rhs) const {return value <  rhs.value;}
    bool operator> (const chemistryLoad& rhs) const {return value >  rhs.value;}
};

} //namespace Foam