#include "HW_3.h"

#include <vector>

#include "../ImageProcessing/ImageProcessing.hpp"

void HW_3_1()
{
    dip::Image chessboard_distorted_256, cat_distorted_256, img;
    chessboard_distorted_256.openRAW("chessboard_distorted_256.raw", dip::Size(256, 256));
    cat_distorted_256.openRAW("cat_distorted_256.raw", dip::Size(256, 256));

    dip::Coordinate<int> xy[4], uv[4];
    xy[0].setCoord(1, 149);
    xy[1].setCoord(105, 1);
    xy[2].setCoord(254, 155);
    xy[3].setCoord(50, 254);
    uv[0].setCoord(0, 0);
    uv[1].setCoord(255, 0);
    uv[2].setCoord(255, 255);
    uv[3].setCoord(0, 255);

    dip::PerspectiveTransformation(chessboard_distorted_256, img, xy, uv);
    img.saveRAW("HW_3_1_chessboard.raw", true);
    dip::PerspectiveTransformation(cat_distorted_256, img, xy, uv);
    img.saveRAW("HW_3_1_cat.raw", true);
}
