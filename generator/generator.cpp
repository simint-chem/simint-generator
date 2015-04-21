#include <iostream>

#include "Classes.h"

using namespace std;

int main(void)
{
    cout << "Here\n";
    Gaussian g1{2,0,0};
    Gaussian g2{1,0,2};
    cout << "g1: " << g1.str() << "\n";
    cout << "g2: " << g2.str() << "\n";
    cout << " <: " << bool(g1 < g2) << "\n";
    cout << "==: " << bool(g1 == g2) << "\n";
//    cout << "!=: " << bool(g1 != g2) << "\n";

    Doublet d{DoubletType::BRA, g1, g2};
    cout << "d: " << d.str() << "\n"; 

    return 0;
}
