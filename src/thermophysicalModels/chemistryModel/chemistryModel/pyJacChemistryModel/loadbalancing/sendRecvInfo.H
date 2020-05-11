#pragma once

#include <DynamicList.H>

namespace Foam{

struct sendRecvInfo{
    
    int my_rank;

    DynamicList<int>    source; //ranks which send to this process
    DynamicList<int>    destinations; //ranks to which this process sends to

    DynamicList<int>    number_of_problems; //number of problems which this rank sends/receivs

};

}

