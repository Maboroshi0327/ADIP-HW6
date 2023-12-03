#include "HW_1.h"

#include "../ImageProcessing/ImageProcessing.hpp"

void HW_1_1()
{
    dip::Image lena, ALNRF;

    lena.openRAW("lean_256_salt&pepper.raw", dip::Size(256, 256));
    dip::AdaptiveLocalNoiseReductionFiltering(lena, ALNRF, dip::Size(5, 5));
    ALNRF.saveRAW("lena_256_salt&pepper_ALNRF.raw", true);

    lena.openRAW("lena_gussuian_1_256.raw", dip::Size(256, 256));
    dip::AdaptiveLocalNoiseReductionFiltering(lena, ALNRF, dip::Size(5, 5));
    ALNRF.saveRAW("lena_gussuian_1_256_ALNRF.raw", true);

    lena.openRAW("lena_gussuian_2_256.raw", dip::Size(256, 256));
    dip::AdaptiveLocalNoiseReductionFiltering(lena, ALNRF, dip::Size(5, 5));
    ALNRF.saveRAW("lena_gussuian_2_256_ALNRF.raw", true);
}

void HW_1_2()
{
    dip::Image lena, ATMF;

    lena.openRAW("lean_256_salt&pepper.raw", dip::Size(256, 256));
    dip::AlphaTrimmedMeanFiltering(lena, ATMF, dip::Size(5, 5), 16);
    ATMF.saveRAW("lena_256_salt&pepper_ATMF.raw", true);

    lena.openRAW("lena_gussuian_1_256.raw", dip::Size(256, 256));
    dip::AlphaTrimmedMeanFiltering(lena, ATMF, dip::Size(5, 5), 16);
    ATMF.saveRAW("lena_gussuian_1_256_ATMF.raw", true);

    lena.openRAW("lena_gussuian_2_256.raw", dip::Size(256, 256));
    dip::AlphaTrimmedMeanFiltering(lena, ATMF, dip::Size(5, 5), 16);
    ATMF.saveRAW("lena_gussuian_2_256_ATMF.raw", true);
}
