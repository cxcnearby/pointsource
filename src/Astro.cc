#include "Astro.h"
#include "constant.h"

#include <math.h>
#include <stdio.h>

/********** GMST **********/

double gmst(double mjd) {

  double ut1, Tu, alpham;

  ut1 = mjd - floor(mjd);
  Tu = (mjd - 51544.5) / 36525.0;
  alpham = 280.460618370 + 36000.770053608 * Tu + 3.8793333331e-4 * Tu * Tu -
           2.5833333331e-8 * Tu * Tu * Tu;

  return (180.0 + ut1 * 360.0 + alpham);
}

double azirange(double azi) {
  while (azi > 180. || azi < -180.) {
    if (azi > 180.) {
      azi -= 360.;
    } else {
      azi += 360.;
    }
  }
  return azi;
}

/********** transfomation of (zenith & azimuth)
            from (decrination & right asension)  **********/

void equator_horizon(double mjd, double ras, double dec, double *zen,
                     double *azi)
// double mjd, ras, dec, *zen, *azi;
{

  double lst, H, g;
  double sH, sdec, sphi;
  double cH, cdec, cphi;
  double zenith, azimuth;
  double sazi, cazi;

  g = gmst(mjd) + observer_lo;
  lst = (g - floor(g / 360.0) * 360.0);
  H = lst - ras;

  sH = sin(H * D2R);
  sdec = sin(dec * D2R);
  sphi = sin(observer_la * D2R);
  cH = cos(H * D2R);
  cdec = cos(dec * D2R);
  cphi = cos(observer_la * D2R);

  zenith = acos(sphi * sdec + cphi * cdec * cH) * R2D;

  sazi = (-cdec * sH) / sin(zenith * D2R);
  cazi = (cphi * sdec - sphi * cdec * cH) / sin(zenith * D2R);

  azimuth = atan2(sazi, cazi) * R2D;
  if (azimuth > 180.0)
    azimuth -= 360.0;
  if (azimuth < -180.0)
    azimuth += 360.0;

  *zen = zenith;
  *azi = azimuth;
}

/********** transfomation of (decrination & right asension)
                                    from (zenith & azimuth) **********/

void horizon_equator(double mjd, double zen, double azi, double *ras,
                     double *dec)
// double mjd, *ras, *dec, zen, azi;
{

  double lst, H, g;
  double sazi, szen, sphi;
  double cazi, czen, cphi;
  double DEC, RAS;
  double sH, cH;

  g = gmst(mjd) + observer_lo;
  lst = (g - floor(g / 360.0) * 360.0);

  sazi = sin(azi * D2R);
  szen = sin(zen * D2R);
  sphi = sin(observer_la * D2R);
  cazi = cos(azi * D2R);
  czen = cos(zen * D2R);
  cphi = cos(observer_la * D2R);

  DEC = asin(sphi * czen + cphi * szen * cazi) * R2D;

  sH = (-szen * sazi) / cos(DEC * D2R);
  cH = (cphi * czen - sphi * szen * cazi) / cos(DEC * D2R);

  H = atan2(sH, cH) * R2D;

  RAS = lst - H;
  if (RAS > 360.0)
    RAS -= 360.0;
  if (RAS < 0.0)
    RAS += 360.0;

  *ras = RAS;
  *dec = DEC;
}

/********** caluculation of distance
 *                          by (altitude & longitude)
 *                             altitude=0 at equator!!
 *                                                   *********/

double distance_equatorial(double ra1, double dec1, double ra2, double dec2) {
  return (acos(cos((ra2 - ra1) * D2R) * cos(dec1 * D2R) * cos(dec2 * D2R) +
               sin(dec1 * D2R) * sin(dec2 * D2R)) *
          R2D);
}

/********** caluculation of distance
 *                          by (zenith & azimuth)
 *                             zenith=0 at north pole!!
 *                                                   *********/

double distance_horizontal(double zen1, double azi1, double zen2, double azi2) {
  return (acos(cos((azi1 - azi2) * D2R) * sin(zen1 * D2R) * sin(zen2 * D2R) +
               cos(zen1 * D2R) * cos(zen2 * D2R)) *
          R2D);
}

/********** caluculation of position angle
 *
 *		the direction is point(ra2,dec2) relative to point(ra1,dec1)
 *		angle=0 when the direction is to the north pole, and clockwise.
 *
 * **********/

double direction_equatorial(double ra1, double dec1, double ra2, double dec2) {
  return (
      atan2(cos(dec2 * D2R) * sin((ra2 - ra1) * D2R),
            cos(dec1 * D2R) * sin(dec2 * D2R) -
                sin(dec1 * D2R) * cos(dec2 * D2R) * cos((ra2 - ra1) * D2R)) *
      R2D);
}

/********** caluculation of position angle
 *  *
 *   *              the direction is point(ra2,dec2) relative to point(ra1,dec1)
 *    *              angle=0 when the direction is to the north pole, and
 * clockwise.
 *     *
 *      * **********/

double direction_horizontal(double zen1, double azi1, double zen2,
                            double azi2) {
  return (
      atan2(sin(zen2 * D2R) * sin((azi1 - azi2) * D2R),
            sin(zen1 * D2R) * cos(zen2 * D2R) -
                cos(zen1 * D2R) * sin(zen2 * D2R) * cos((azi1 - azi2) * D2R)) *
      R2D);
}
