#pragma once


///
///@brief A generic load object
///
///
struct Load{

    int rank;     //rank of the process 
    double value; //an indicator how busy the processor is


    //overload comparisons for sorting
    bool operator==(const Load& rhs) const {return value == rhs.value;}
    bool operator<=(const Load& rhs) const {return value <= rhs.value;}
    bool operator>=(const Load& rhs) const {return value >= rhs.value;}
    bool operator< (const Load& rhs) const {return value <  rhs.value;}
    bool operator> (const Load& rhs) const {return value >  rhs.value;}
};