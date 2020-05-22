#ifndef POINTSOURCE_INCLUDE_CONSTANT_H_
#define POINTSOURCE_INCLUDE_CONSTANT_H

#include <vector>

/********** Constant **********/

extern const double D2R; /*   pi/180  */
extern const double R2D; /*   180/pi  */
extern const double PI;  /* pi */

/********** Location of Tibet **********/

extern const double observer_la; /*  latitude  */
extern const double observer_lo; /* longtitude */

/********** simulation ***********/
extern const double kArea;
extern std::vector<int> kTypeCode;
extern std::vector<double> kTypePowerIndex;
extern std::vector<double> kEnergyBounds;
extern const int kNEbin;
extern const int kNESmallbin;
extern const int kNEnergyBin;

extern const double kZenBinWidth;
extern const double kZenRange;
extern const int kNZenBin;

extern const int kPmtBinWidth;
extern const int kPmtRange;
extern const int kNPmtBin;

#endif