/*

*/

#pragma once

namespace Foam{


///
///@brief A generic chemistry problem. These are communicated through the load balancer.
///       essentially this should contain all parameters that the chemistrySolver.solve() requires
///
///
struct chemistryProblem{

    double value; //some value 
};

} //namespace Foam