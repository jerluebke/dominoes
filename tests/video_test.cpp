#include "../include/DominoChain.hpp"


int main()
{
    domino d;
    d.height = 2;
    d.width = 0.4;
    DominoChain dc( d );
    std::cerr << "makeing video ...\n";
    return dc.make_video( "./video_04.avi", 6.0, 0.8, 0.2 );
}
