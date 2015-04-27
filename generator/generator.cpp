#include <iostream>
#include <memory>
#include <stdexcept>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/Create.hpp"

using namespace std;


int main(int argc, char ** argv)
{
    try {

    if(argc != 5)
    {
        printf("Need 4 arguments\n");
        return 1;
    }

    std::array<int, 4> amlist{
                               atoi(argv[1]),
                               atoi(argv[2]),
                               atoi(argv[3]),
                               atoi(argv[4])
                             };

    std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR);
    Create_Unrolled(amlist, hrralgo, cout);


    }
    catch(std::exception & ex)
    {
        cout << "\n\n";
        cout << "Caught exception\n";
        cout << "What = " << ex.what() << "\n\n";
        return 100;
    }
    return 0;
}
