#ifndef POINTSOURCE_INCLUDE_CONSTANT_H_
#define POINTSOURCE_INCLUDE_CONSTANT_H
/********** Constant **********/

#define D2R 0.017453292519943296  /*   pi/180  */
#define R2D 57.295779513082322    /*   180/pi  */
#define PI 3.14159265358979323846 /* pi */

/********** Location of Tibet **********/

#define observer_la 29.357656306  /*  latitude  */
#define observer_lo 100.138794639 /* longtitude */

const double zen_bin_width = 0.5;
const int n_zen_bin = 60. / zen_bin_width;
double area = 2000. * 2000.;

#endif