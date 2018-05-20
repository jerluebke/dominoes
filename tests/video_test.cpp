#include "../include/DominoChainVideo.hpp"


int main()
{
    domino d;
    d.height = 2;
    d.width = 0.4;
    DominoChainVideo dc( d );
    return dc.make_video( "./video_01.avi", 6.0, 0.8, 0.2 );
}
