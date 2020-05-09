#ifndef POINTSOURCE_INCLUDE_ASTRO_H_
#define POINTSOURCE_INCLUDE_ASTRO_H_

void moon_orbit(double mjd, double *MRAS, double *MDEC);
void sun_orbit(double mjd, double *SRAS, double *SDEC);
double gmst(double mjd);
double azirange(double azi);
void equator_horizon(double mjd, double ras, double dec, double *zen,
                     double *azi);
void horizon_equator(double mjd, double zen, double azi, double *ras,
                     double *dec);
double distance_equatorial(double ra1, double dec1, double ra2, double dec2);
double distance_horizontal(double zen1, double azi1, double zen2, double azi2);
double direction_equatorial(double ra1, double dec1, double ra2, double dec2);
double direction_horizontal(double zen1, double azi1, double zen2, double azi2);

#endif