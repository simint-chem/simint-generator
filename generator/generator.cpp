#include <iostream>

#include "Classes.h"
#include "Algorithms.h"
#include "Helpers.h"

using namespace std;

int main(void)
{
    // generate initial targets
    QuartetSet inittargets = GenerateInitialTargets({2,2,0,0});

    for(auto & it : inittargets)
        cout << it << "\n";

    return 0;
}
