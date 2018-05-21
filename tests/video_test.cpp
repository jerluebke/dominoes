#include "../include/DominoChain.hpp"


int main()
{
    domino d;
    d.height = 2;
    d.width = 0.4;
    DominoChain dc( d );
    std::cerr << "making video ...\n";
    return dc.make_video( "./video_14.avi", 6.0, 1.9, 0.2 );
}
