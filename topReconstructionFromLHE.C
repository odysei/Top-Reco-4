#include "topReconstructionFromLHE.h"

using namespace std;

int main()
{
    topReconstructionFromLHE t;
    t.debug_verbosity = 3;
    t.Loop("output_files", 0, 1, 0, 6);

    //     t.Plot("plots");
    // t.Print();

    return 0;
}
