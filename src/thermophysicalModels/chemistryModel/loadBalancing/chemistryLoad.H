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

    chemistryLoad(int _rank, scalar _value) : 
    rank(_rank), 
    value(_value) 
    {} 


    int rank;     //rank of the process 
    scalar value; //an indicator how busy the processor is

    //overload comparisons for sorting
    bool operator==(const chemistryLoad& rhs) const {return value == rhs.value;}
    bool operator!=(const chemistryLoad& rhs) const {return !(value == rhs.value);}
    bool operator<=(const chemistryLoad& rhs) const {return value <= rhs.value;}
    bool operator>=(const chemistryLoad& rhs) const {return value >= rhs.value;}
    bool operator< (const chemistryLoad& rhs) const {return value <  rhs.value;}
    bool operator> (const chemistryLoad& rhs) const {return value >  rhs.value;}


    bool operator==(scalar rhs) const {return value == rhs;}
    bool operator!=(scalar rhs) const {return !(*this == rhs);}
    bool operator<=(scalar rhs) const {return value <= rhs;}
    bool operator>=(scalar rhs) const {return value >= rhs;}
    bool operator<(scalar rhs) const {return value < rhs;}
    bool operator>(scalar rhs) const {return value > rhs;}


    chemistryLoad& operator+=(const chemistryLoad& rhs) {value += rhs.value; return *this;}    
    chemistryLoad& operator-=(const chemistryLoad& rhs) {value -= rhs.value; return *this;}    
    chemistryLoad& operator*=(const chemistryLoad& rhs) {value *= rhs.value; return *this;}    
    chemistryLoad& operator/=(const chemistryLoad& rhs) {value /= rhs.value; return *this;}

    chemistryLoad operator+(scalar rhs) {chemistryLoad l(*this); return l+=rhs;}    
    chemistryLoad operator-(scalar rhs) {chemistryLoad l(*this); return l-=rhs;}    
    chemistryLoad operator*(scalar rhs) {chemistryLoad l(*this); return l*=rhs;}    
    chemistryLoad operator/(scalar rhs) {chemistryLoad l(*this); return l/=rhs;}    

    chemistryLoad& operator+=(scalar rhs) {value += rhs; return *this;}    
    chemistryLoad& operator-=(scalar rhs) {value -= rhs; return *this;}    
    chemistryLoad& operator*=(scalar rhs) {value *= rhs; return *this;}    
    chemistryLoad& operator/=(scalar rhs) {value /= rhs; return *this;}



};

Ostream& operator<<(Ostream& os, const chemistryLoad& l);
Istream& operator>>(Istream& is, chemistryLoad& l);





} //namespace Foam