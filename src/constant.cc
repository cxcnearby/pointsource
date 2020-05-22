#include "constant.h"

const double D2R = 0.017453292519943296;  /*   pi/180  */
const double R2D = 57.295779513082322;    /*   180/pi  */
const double PI = 3.14159265358979323846; /* pi */

/********** Location of Tibet **********/

const double observer_la = 29.357656306;  /*  latitude  */
const double observer_lo = 100.138794639; /* longtitude */

// simulation parameters
const double kArea = 2000. * 2000.;
std::vector<int> kTypeCode = {0, 14, 402, 1407, 2513, 5626};
std::vector<double> kTypePowerIndex = {-2.62, -2.71, -2.64,
                                       -2.68, -2.66, -2.59};
std::vector<double> kEnergyBounds = {0.01, 0.1, 1., 1.e1, 1.e2, 1.e3};
const int kNEbin = kEnergyBounds.size() - 1;
const int kNESmallbin = 4;
const int kNEnergyBin = kNESmallbin * kNEbin;

// histogram parameters
const double kZenBinWidth = 0.5;
const double kZenRange = 60.;
const int kNZenBin = kZenRange / kZenBinWidth;

const int kPmtBinWidth = 10;
const int kPmtRange = 900;
const int kNPmtBin = kPmtRange / kPmtBinWidth;
