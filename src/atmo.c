/* ============================================================================

   Copyright (C) 1990, 1997, 1998, ..., 2010  Konrad Bernloehr

   This file is part of the IACT/atmo package for CORSIKA.

   The IACT/atmo package is free software; you can redistribute it 
   and/or modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This package is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this package. If not, see <http://www.gnu.org/licenses/>.

============================================================================ */

/* ================================================================== */
/**
 *  @file atmo.c
 *  @short Use of tabulated atmospheric profiles and atmospheric refraction.
 *
 *  @author  Konrad Bernloehr 
 *  @date    @verbatim CVS $Date: 2015/01/17 16:22:03 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.3 $ @endverbatim
 *
 *  --------------------------------------------------------------------
 *
 *  This file provides code for use of external atmospheric models
 *  (in the form of text-format tables) with the CORSIKA program.
 *  Six atmospheric models as implemented in the MODTRAN program
 *  and as tabulated in MODTRAN documentation (F.X. Kneizys et al. 1996,
 *  'The MODTRAN 2/3 Report and LOWTRAN 7 Model', Phillips Laboratory,
 *  Hanscom AFB, MA 01731-3010, U.S.A.) are provided as separate files
 *  (atmprof1.dat ... atmprof6.dat). User-provided atmospheric
 *  models should be given model numbers above 6.
 *
 *  Note that for the Cherenkov part and the hadronic (and muon) part
 *  of CORSIKA the table values are directly interpolated but the
 *  electron/positron/gamma part (derived from EGS) uses special
 *  layers (at present 4 with exponential density decrease and the
 *  most upper layer with constant density). Parameters of these
 *  layers are fitted to tabulated values but not every possible
 *  atmospheric model fits very well with an exponential profile.
 *  You are adviced to check that the fit matches tabulated values to
 *  sufficient precision in the altitude ranges of interest to you.
 *  Try to adjust layer boundary altitudes in case of problems.
 *  The propagation of light without refraction (as implemented in
 *  CORSIKA, unless using the CURVED option) and with refraction (as
 *  implemented by this software) assumes a plane-parallel atmosphere.
*/
/* ==================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../inc/atmo.h"
#include "../inc/fileopen.h"

/*
 * The fast interpolation method sidesteps calls to LOG and EXP functions
 * by a narrower internal step size but can result in little artifacts
 * around the sampling points of the atmospheric profile.
 * To compile with fast interpolation disabled add '-DNO_FAST_INTERPOLATION'
 * to your compilation flags.
 */

#ifndef NO_FAST_INTERPOLATION
#define FAST_INTERPOLATION 1
#endif

/* Local C called functions (parameter types are always checked) */
static void interp(double x, double *v, int n, int *ipl, double *rpl);
static void init_refraction_tables(void);
static void init_fast_interpolation(void);
// static void init_corsika_atmosphere(void);
static void init_atmosphere(void);
static double sum_log_dev_sq(double a, double b, double c, int np, 
   double *h, double *t, double *rho);
static double atm_exp_fit(double h1, double h2, double *ap, double *bp, 
   double *cp, double *s0, int *npp);

/* Variables used for atmospheric profiles */

#define MAX_PROFILE 50
int atmosphere;  ///< The atmospheric profile number, 0 for built-in.
static int num_prof;
static double p_alt[MAX_PROFILE], p_log_alt[MAX_PROFILE];
static double p_log_rho[MAX_PROFILE], p_rho[MAX_PROFILE];
static double p_log_thick[MAX_PROFILE];
static double p_log_n1[MAX_PROFILE];
static double p_bend_ray_hori_a[MAX_PROFILE];
static double p_bend_ray_time0[MAX_PROFILE];
static double p_bend_ray_time_a[MAX_PROFILE];

static double top_of_atmosphere = 112.83e5;
static double bottom_of_atmosphere = 0.;

#ifdef FAST_INTERPOLATION
#define MAX_FAST_PROFILE 10000
static double fast_p_alt[MAX_FAST_PROFILE];
static double fast_p_log_rho[MAX_FAST_PROFILE];
static double fast_p_log_thick[MAX_FAST_PROFILE];
static double fast_p_log_n1[MAX_FAST_PROFILE];
static double fast_h_fac;
#endif

/* ================================================================== */
/*
   Linear interpolation functions from an older program of K.B.
   A binary search algorithm is used for fast interpolation.
*/

/* --------------------------- interp ------------------------------- */
/**
 *  @short Linear interpolation with binary search algorithm.
 *
 *  Linear interpolation between data point in sorted (i.e. monotonic
 *  ascending or descending) order. This function determines between
 *  which two data points the requested coordinate is and where between
 *  them. If the given coordinate is outside the covered range, the
 *  value for the corresponding edge is returned.
 *
 *  A binary search algorithm is used for fast interpolation.
 *
 *  @param  x Input: the requested coordinate
 *  @param  v Input: tabulated coordinates at data points
 *  @param  n Input: number of data points
 *  @param  ipl Output: the number of the data point following the requested
 *	    coordinate in the given sorting (1 <= ipl <= n-1)
 *  @param  rpl Output: the fraction (x-v[ipl-1])/(v[ipl]-v[ipl-1])
 *	    with 0 <= rpl <= 1
*/      

static void interp ( double x, double *v, int n, int *ipl, double *rpl )
{
   int i, l, m, j, lm;

#ifdef DEBUG_TEST_ALL
   if ( v == NULL || n <= 2 )
   {
      fprintf(stderr,"Invalid parameters for interpolation.\n");
      *ipl = 1;
      *rpl = 0.;
      return;
   }
#endif

   if ( v[0] < v[n-1] )
   {
      if (x <= v[0])
      {
         *ipl = 1;
         *rpl = 0.;
         return;
      }
      else if (x >= v[n-1])
      {
         *ipl = n-1;
         *rpl = 1.;
         return;
      }
      lm = 0;
   }
   else
   {
      if (x >= v[0])
      {
         *ipl = 1;
         *rpl = 0.;
         return;
      }
      else if (x <= v[n-1])
      {
         *ipl = n-1;
         *rpl = 1.;
         return;
      }
      lm = 1;
   }

   l = (n+1)/2-1;
   m = (n+1)/2;
   for (i=1; i<=30; i++ )
   {
      j = l;
      if (j < 1) j=1;
      if (j > n-1) j=n-1;
      if (x >= v[j+lm-1] && x <= v[j-lm])
      {
         *ipl = j;
         if ( v[j] != v[j-1] )
            *rpl = (x-v[j-1])/(v[j]-v[j-1]);
         else
            *rpl = 0.5;
         return;
      }
      m = (m+1)/2;
      if (x > v[j-1])
         l = l + (1-2*lm)*m;
      else
         l = l - (1-2*lm)*m;
   }
   fprintf(stderr,"Interpolation error.\n");
}

/* ----------------------------- rpol ------------------------------- */
/**
 *  @short Linear interpolation with binary search algorithm.
 *
 *  Linear interpolation between data point in sorted (i.e. monotonic
 *  ascending or descending) order. The resulting interpolated value
 *  is returned as a return value.
 *
 *  This function calls interp() to find out where to interpolate.
 *  
 *  @param   x  Input: Coordinates for data table
 *  @param   y  Input: Corresponding values for data table
 *  @param   n  Input: Number of data points
 *  @param   xp Input: Coordinate of requested value
 *
 *  @return  Interpolated value
 *
*/

double rpol ( double *x, double *y, int n, double xp )
{
   int ipl = 1;
   double rpl = 0.;

   interp ( xp, x, n, &ipl, &rpl );
   return y[ipl-1]*(1.-rpl) + y[ipl]*rpl;
}

/* ======================================================================= */

static double etadsn;  /**< About the same as in CORSIKA Cherenkov function */
                       /**< (but doesn't need to be the same). */
static double observation_level;  /**< Altitude [cm] of observation level */
static double obs_level_refidx;
static double obs_level_thick;

/* ------------------- init_refraction_tables ---------------------- */
/**
 *  @short Initialize tables needed for atmospheric refraction.
 *
 *  Initialize the correction tables used for the refraction bending
 *  of the light paths. It is called once after the atmospheric
 *  profile has been defined.
*/

static void init_refraction_tables()
{
   int ialt;
   
   /* Etadsn is the parameter used in CORSIKA for the scaling */
   /* between density and index of refraction minus one (n-1). */
   etadsn = 0.000283 * 994186.38 / 1222.656; /* CORSIKA default */
   /* Look for better approximation above the observation level. */
   for (ialt=0; ialt<num_prof; ialt++)
      if ( p_alt[ialt] > observation_level + 1.5e5 )
      {
         etadsn = exp(p_log_n1[ialt]) / exp(p_log_rho[ialt]);
         break;
      }
      
#ifdef DEBUG_TEST_ALL
   if ( p_alt[0] != p_alt[4] && p_log_rho[0] != p_log_rho[4] )
   {
      double dscl_p, dscl_t;
      dscl_p = (p_alt[4]-p_alt[0]) / (p_log_rho[0]-p_log_rho[4]);
      dscl_t = (p_alt[4]-p_alt[0]) / (p_log_thick[0]-p_log_thick[4]);
      printf(" Pressure scale height near ground: %7.2f m\n",dscl_p*0.01);
      printf(" Thickness scale height near ground: %7.2f m\n",dscl_t*0.01);
      printf(" Etadsn=(n-1)/density parameter: %g\n",etadsn);
   }
#endif
   
   /* Initialize tables by numerical integration for vertical and slanted */
   /* (45 degrees) paths, taking into account known angular dependencies. */
   
   for (ialt=0; ialt<num_prof; ialt++)
   {
      double t0_vt, t0_45, x0_45, n_sint0;
      double t_vt, t_45, x_45, dt_vt, dt_45, dx_45, dz, z, zm, ds;
      double vc = 29.9792458; /* velocity of light [cm/ns] */
      double theta0 = 45. * (M_PI/180.), theta1, theta2;
      int nz, iz;
      double c, s;
      
      x0_45 = (p_alt[ialt]-observation_level) * tan(theta0);
      t0_vt = 1./vc * (p_alt[ialt] - observation_level +
           etadsn*(thickx_(&observation_level)-thickx_(&p_alt[ialt])));
      t0_45 = t0_vt / cos(theta0);
      nz = 1000; /* Number of steps for numerical ray tracing */
      dz = (observation_level-p_alt[ialt]) / (double)nz;
      n_sint0 = refidx_(&p_alt[ialt]) * sin(theta0);
      for (iz=0, z=p_alt[ialt], theta2=theta0, t_vt=t_45=x_45=0.; iz<nz; iz++)
      {
         z += dz;
         zm = z-0.5*dz;
         theta1 = theta2;
         theta2 = asin(n_sint0/refidx_(&z));
         ds = fabs(dz) / cos(0.5*(theta1+theta2));
         dt_vt = fabs(dz) * refidx_(&zm) / vc;
         dt_45 = ds * refidx_(&zm) / vc;
         dx_45 = fabs(dz) * tan(0.5*(theta1+theta2));
         t_vt += dt_vt;
         t_45 += dt_45;
         x_45 += dx_45;
      }
      if ( p_alt[ialt] < observation_level )
      {
         t_vt *= -1.;
         t_45 *= -1.;
         x_45 *= -1.;
      }
      theta1 = asin(n_sint0/refidx_(&observation_level));
      c = cos(theta0+0.28*(theta1-theta0));
      s = sin(theta0+0.28*(theta1-theta0));
      if ( x_45 < x0_45 ) /* Offset is normally less than for straight line */
         p_bend_ray_hori_a[ialt] = sqrt((x0_45 - x_45) * (c*c*c)/s);
      else
         p_bend_ray_hori_a[ialt] = 0.;
      p_bend_ray_time0[ialt]  = t_vt - t0_vt;
      if ( ((t_45 - t0_45) - (t_vt - t0_vt)) < 0. )
         p_bend_ray_time_a[ialt] = 
          sqrt(((t0_45 - t_45) - (t0_vt - t_vt)) * (c*c*c)/(s*s));
      else
         p_bend_ray_time_a[ialt] = 0.;
   }
}

/* ------------------- init_fast_interpolation ---------------------- */

#ifdef FAST_INTERPOLATION
/** 
 *  @short An alternate interpolation method (which requires that the
 *         table is sufficiently fine-grained and equidistant) has to
 *         be initialized first.
 */
static void init_fast_interpolation()
{
   int i;
   for ( i=0; i<MAX_FAST_PROFILE; i++)
   {
      if ( i<MAX_FAST_PROFILE-1 )
         fast_p_alt[i] = bottom_of_atmosphere + (double) i /
            (double)(MAX_FAST_PROFILE-1) * 
            (top_of_atmosphere - bottom_of_atmosphere);
      else /* avoid rounding errors */
         fast_p_alt[i] = top_of_atmosphere;
      fast_p_log_rho[i]   = rpol(p_alt,p_log_rho,num_prof,fast_p_alt[i]);
      fast_p_log_thick[i] = rpol(p_alt,p_log_thick,num_prof,fast_p_alt[i]);
      fast_p_log_n1[i]    = rpol(p_alt,p_log_n1,num_prof,fast_p_alt[i]);
   }
   
   fast_h_fac = (double)(MAX_FAST_PROFILE-1) / 
        (top_of_atmosphere - bottom_of_atmosphere);
}
#endif

/* ------------------- init_corsika_atmosphere -------------------- */
/**
 *  @short Take the atmospheric profile from CORSIKA built-in functions.
 *
 *  For use of the refraction bending corrections together with the
 *  CORSIKA built-in atmospheres, the atmosphere tables are constructed
 *  from the CORSIKA RHOF and THICK functions.
 *  Note that the refraction index in this case is without taking the
 *  effect of wator vapour into account.
*/

/*

 removed by GM

*/


/* ----------------------- init_atmosphere ------------------------ */
/**
 *  @short Initialize atmospheric profiles.
 *
 *  Internal function for initialising both external and CORSIKA
 *  built-in atmospheric profiles. If any CORSIKA built-in profile
 *  should be used, it simply calls init_corsika_atmosphere().
 *
 *  Otherwise, atmospheric models are read in from text-format tables.
 *  The supplied models 1-6 are based on output of the MODTRAN program.
 *  For the interpolation of relevant parameters (density, thickness,
 *  index of refraction, ...) all parameters are transformed such
 *  that linear interpolation can be easily used.
 *
*/

static void init_atmosphere ()
{
   char fname[128];
   FILE *f;
   char line[1024];
   int count;
   double alt,rho,thick,n_1;
#ifdef LONG_ATMPROF
   double p,t,N,O3,H2O;
#endif
   
   /* CORSIKA built-in atmospheres have atmosphere numbers <= 0 */
   if ( atmosphere <=0 )
   {
//      init_corsika_atmosphere();
      return;
   }
   
   /* There are two different versions of data files. */
#ifndef LONG_ATMPROF
   sprintf(fname,"data/atmprof%d.dat",atmosphere);
#else
   sprintf(fname,"atm_profile_model_%d.dat",atmosphere);
#endif
   if ( (f=fileopen(fname,"r")) == NULL )
   {
      perror(fname);
#ifdef LONG_ATMPROF /* Try the other variant before giving up. */
      sprintf(fname,"data/atmprof%d.dat",atmosphere);
#else
      sprintf(fname,"atm_profile_model_%d.dat",atmosphere);
#endif
      fprintf(stderr,"Trying file %s instead.\n", fname);
      if ( (f=fileopen(fname,"r")) == NULL )
      {
         perror(fname);
         exit(1);
      }
   }

   count = num_prof = 0;
   while ( fgets(line,sizeof(line)-1,f) != NULL && num_prof < MAX_PROFILE )
   {
      char *s;
      count++;

      for (s=line;*s==' ';s++)
         ;
      if ( *s=='#' ) /* Comment line */
         continue;
#ifndef LONG_ATMPROF
      /* The short files contain only data relevant for CORSIKA. */
      if ( sscanf(s,"%lf %lf %lf %lf",
        &alt,&rho,&thick,&n_1) != 4 )
#else
      /* The long files contain other data as well. */
      if ( sscanf(s,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &alt,&p,&t,&N,&rho,&thick,&O3,&H2O,&n_1) != 9 )
#endif
      {
         fprintf(stderr,"Syntax error in %s line %d.\n",fname,count);
         exit(1);
      }
      
      p_alt[num_prof] = alt*1e5; /* Altitude in file was in km */
      p_log_alt[num_prof] = (alt>0.)?log(alt*1e5):0.;
      p_log_rho[num_prof] = (rho>0.)?log(rho):-1000.;
      p_rho[num_prof] = rho;
      p_log_thick[num_prof] = (thick>0.)?log(thick):-1000.;
      p_log_n1[num_prof] = (n_1>0.)?log(n_1):-1000.;
      num_prof++;
   }
   
   fclose(f);
   fflush(stdout);
//   printf("\n Atmospheric profile %d with %d levels read from file %s\n\n",
//      atmosphere,num_prof,fname);

   if ( num_prof < 5 )
   {
      fprintf(stderr,
         "There are definitely too few atmospheric levels in this file.\n");
      fprintf(stderr,
         "Normally this kind of file should have 50 levels.\n");
      exit(1);
   }
   
   bottom_of_atmosphere = p_alt[0];
   top_of_atmosphere    = p_alt[num_prof-1];

#ifdef FAST_INTERPOLATION
   /* Initialize faster tables for most frequent lookups */
   init_fast_interpolation();
#endif

   /* Initialize the tables for the refraction bending */
   init_refraction_tables();
}

/* -------------------------- atmset_ ---------------------------- */
/**
 *  @short Set number of atmospheric model profile to be used.
 *
 *  The atmospheric model is initialized first before the
 *  interpolating functions can be used. For efficiency reasons,
 *  the functions rhofx_(), thickx_(), ... don't check if the
 *  initialisation was done.
 *
 *  This function is called if the 'ATMOSPHERE' keyword is
 *  present in the CORSIKA input file.
 *
 *  The function may be called from CORSIKA to initialize
 *  the atmospheric model via 'CALL ATMSET(IATMO,OBSLEV)' or such.
 *
 *  @param iatmo   (pointer to) atmospheric profile number;
 *   	           negative for CORSIKA built-in profiles.
 *  @param obslev  (pointer to) altitude of observation level [cm]
 *
 *  @return (none)
*/

void atmset_ (int *iatmo, double *obslev)
{
   atmosphere = *iatmo;
   observation_level = *obslev;
   init_atmosphere();
   obs_level_refidx = refidx_(obslev);
   obs_level_thick = thickx_(obslev);
}

/* ---------------------------- rhofx_ ----------------------------- */
/**
 *
 *  @short Density of the atmosphere as a function of altitude.
 *  This function can be called from Fortran code as RHOFX(HEIGHT).
 *
 *  @param  height (pointer to) altitude [cm]
 *
 *  @return density [g/cm**3]
*/

double rhofx_ (double *height)
{
#ifdef FAST_INTERPOLATION
   int i;
   double r;
   if ( (*height) < bottom_of_atmosphere )
      return p_rho[0];
   else if ( (*height) >= top_of_atmosphere )
      return 0.;
   i = (int) ( fast_h_fac * ((*height)-bottom_of_atmosphere) );
   if ( i >= MAX_FAST_PROFILE-1 )
      return 0.;
   r = fast_h_fac * ((*height)-fast_p_alt[i]);
   return exp((1.-r)*fast_p_log_rho[i] + r*fast_p_log_rho[i+1]);
#else
   return exp(rpol(p_alt,p_log_rho,num_prof,*height));
#endif
}

/* ---------------------------- thickx_ ----------------------------- */
/**
 *
 *  @short Atmospheric thickness [g/cm**2] as a function of altitude.
 *  This function can be called from Fortran code as THICKX(HEIGHT).
 *
 *  @param  height (pointer to) altitude [cm]
 *
 *  @return thickness [g/cm**2]
 *
*/

double thickx_ (double *height)
{
#ifdef FAST_INTERPOLATION
   int i;
   double r;
   if ( (*height) < bottom_of_atmosphere )
      return exp(fast_p_log_thick[0]);
   else if ( (*height) >= top_of_atmosphere )
      return 0.;
   i = (int) ( fast_h_fac * ((*height)-bottom_of_atmosphere) );
   if ( i >= MAX_FAST_PROFILE-1 )
      return 0.;
   r = fast_h_fac * ((*height)-fast_p_alt[i]);
   return exp((1.-r)*fast_p_log_thick[i] + r*fast_p_log_thick[i+1]);
#else
   return exp(rpol(p_alt,p_log_thick,num_prof,*height));
#endif
}

/* ---------------------------- refidx_ ----------------------------- */
/**
 *
 *  @short Index of refraction as a function of altitude [cm].
 *  This function can be called from Fortran code as REFIDX(HEIGHT).
 *
 *  @param height (pointer to) altitude [cm]
 *
 *  @return index of refraction
 *
*/

double refidx_ (double *height)
{
#ifdef FAST_INTERPOLATION
   int i;
   double r;
   if ( (*height) < bottom_of_atmosphere )
      return 1.+exp(fast_p_log_n1[0]);
   else if ( (*height) >= top_of_atmosphere )
      return 1.;
   i = (int) ( fast_h_fac * ((*height)-bottom_of_atmosphere) );
   if ( i >= MAX_FAST_PROFILE-1 )
      return 1.;
   r = fast_h_fac * ((*height)-fast_p_alt[i]);
   return 1.+exp((1.-r)*fast_p_log_n1[i] + r*fast_p_log_n1[i+1]);
#else
   return 1.+exp(rpol(p_alt,p_log_n1,num_prof,*height));
#endif
}

/* ---------------------------- heighx_ ----------------------------- */
/**
 *
 *  @short Altitude [cm] as a function of atmospheric thickness [g/cm**2].
 *  This function can be called from Fortran code as HEIGHX(THICK).
 *
 *  @param   thick  (pointer to) atmospheric thickness [g/cm**2]
 *
 *  @return   altitude [cm]
*/

double heighx_ (double *thick)
{
   double h;
   if ( (*thick) <= 0. )
      return top_of_atmosphere;
   h = rpol(p_log_thick,p_alt,num_prof,log(*thick));
   if ( h < top_of_atmosphere )
      return h;
   else
      return top_of_atmosphere;
}

/* ---------------------------- raybnd_ ---------------------------- */
/**
 *  @short Calculate the bending of light due to atmospheric refraction.
 *
 *  Path of light through the atmosphere including the bending by refraction.
 *  This function assumes a plane-parallel atmosphere.
 *  Coefficients for corrections from straight-line propagation to
 *  refraction-bent path are numerically evaluated when the atmospheric 
 *  model is defined.
 *  Note that while the former mix of double/float data types may appear odd,
 *  it was determined by the variables present in older CORSIKA to save
 *  conversions. With CORSIKA 6.0 all parameters are of double type.
 *
 *  This function may be called from FORTRAN as
 *    CALL RAYBND(ZEM,U,V,W,DX,DY,DT)
 *
 *  @param zem    Altitude of emission above sea level [cm]
 *  @param u      Initial/Final direction cosine along X axis (updated)
 *  @param v      Initial/Final direction cosine along Y axis (updated)
 *  @param w      Initial/Final direction cosine along Z axis (updated)
 *  @param dx     Position in CORSIKA detection plane [cm] (updated)
 *  @param dy     Position in CORSIKA detection plane [cm] (updated)
 *  @param dt     Time of photon [ns]. Input: emission time.
 *                Output: time of arrival in CORSIKA detection plane.
*/

void raybnd_(double *zem, cors_real_now_t *u, cors_real_now_t *v, 
   double *w, cors_real_now_t *dx, cors_real_now_t *dy, cors_real_now_t *dt)
{
   double sin_t_em, sin_t_obs, theta_em, theta_obs;
   double c, s, h, t, rho;
   double hori_off, travel_time;
   double vc = 29.9792458; /* velocity of light [cm/ns] */
   
   /* (Sine of) emission zenith angle */
   sin_t_em = sqrt((double)((*u)*(*u)+(*v)*(*v)));
   if ( sin_t_em <= 0. )
   {   /* Exactly vertical: no bending; just calulate travel time. */
      *dt += (((*zem) - observation_level) + 
         etadsn*(obs_level_thick-thickx_(zem))) / (*w) / vc +
         rpol(p_alt,p_bend_ray_time0,num_prof,*zem);
      return;
   }
   if ( sin_t_em > 1. || (*w) <= 0. )
      return;
   theta_em = asin(sin_t_em);
   
   /* (Sine of) observed zenith angle */
   sin_t_obs = sin_t_em*refidx_(zem) / obs_level_refidx;
   if ( sin_t_obs > 1. )
      return;
   theta_obs = asin(sin_t_obs);
   
#ifdef TEST_RAYBND
fflush(NULL);
printf(" raybnd: theta = %5.3f, %5.3f\n",(double)(theta_em*180./M_PI),
   (double)(theta_obs*180./M_PI));
#endif

   /* Calculate horizontal displacement with respect to straight line */
   /* and total light travel time from emission to observation level. */
   
   c = cos(theta_em+0.28*(theta_obs-theta_em));
   s = sin(theta_em+0.28*(theta_obs-theta_em));
   
   rho = rhofx_(zem);
   h = rpol(p_rho,p_bend_ray_hori_a,num_prof,rho);
   hori_off = -(h*h) * s/(c*c*c);
   t = rpol(p_rho,p_bend_ray_time_a,num_prof,rho);
#ifdef TEST_RAYBND
printf(" raybnd: horizontal displacement = %5.2f\n",hori_off);
printf(" raybdn: time = %5.3f + %5.3f +%5.3f + %5.3f\n",
   *dt,(((*zem) - observation_level) + 
         etadsn*(obs_level_thick-thickx_(zem))) / (*w) / vc,
   rpol(p_alt,p_bend_ray_time0,num_prof,*zem),
      -(t*t) * (s*s)/(c*c*c));
#endif
   travel_time = rpol(p_alt,p_bend_ray_time0,num_prof,*zem) -
      (t*t) * (s*s)/(c*c*c);
   travel_time += (((*zem) - observation_level) + 
         etadsn*(obs_level_thick-thickx_(zem))) / (*w) / vc;
   
   /* Update arguments: */
   /* Emission direction replaced by observed direction. */
   *u *= sin_t_obs/sin_t_em;
   *v *= sin_t_obs/sin_t_em;
   *w = sqrt(1.-sin_t_obs*sin_t_obs);
   
   /* Position in observation level corrected for displacement. */
   *dx += hori_off * (*u)/sin_t_obs;
   *dy += hori_off * (*v)/sin_t_obs;
   
   /* Light travel time added to emission time. */
   *dt += travel_time;
}

/* ============================================================== */
/*
 *  Functions for fitting the tabulated density profile for CORSIKA EGS part.
*/

/* ------------------------ sum_log_dev_sq ------------------------- */
/**
 *  Measure of deviation of model layers from tables.
*/

static double sum_log_dev_sq(double a, double b, double c, int np,
   double *h, double *t, double *rho)
{
   int ip;
   double s = 0., al;
   for (ip=0; ip<np; ip++ )
   {
#if 0
      /* Fit minimizes relative deviations (better fits at high altitude) */
      if ( a+b*exp(-h[ip]/c) > 0. && t[ip] > 0. )
         al = log((a+b*exp(-h[ip]/c)) / t[ip]);
      else
         al = 0.;
      s += al*al;
#elif 1
      /* Compromise between relative and absolute deviations */
      if ( a+b*exp(-h[ip]/c) > 0. && t[ip] > 0. )
         al = log((a+b*exp(-h[ip]/c)) / t[ip]);
      else
         al = 0.;
      s += fabs(al)*fabs(a+b*exp(-h[ip]/c) - t[ip]);
#elif 0
      /* Fit minimizes absolute deviations (better fits at low altitude) */
      al = a+b*exp(-h[ip]/c) - t[ip];
      s += al*al;
#endif
   }
   return s;
}

/* ------------------------ atm_exp_fit ------------------------- */
/**
 *  Fit one atmosphere layer by an expontential density model.
*/

static double atm_exp_fit ( double h1, double h2, double *ap, 
   double *bp, double *cp, double *s0, int *npp )
{
   int ip, np, iter;
   double h[MAX_PROFILE], t[MAX_PROFILE], rho[MAX_PROFILE], t1, t2;
   double a = *ap, b = *bp, c = *cp;
   double s, dc;
   
   for (ip=np=0; ip<num_prof; ip++)
      if ( p_alt[ip] >= h1 && p_alt[ip] <= h2 && p_alt[ip] < 86e5 )
      {
         h[np] = p_alt[ip];
         t[np] = exp(p_log_thick[ip]);
         rho[np] = exp(p_log_rho[ip]);
         np++;
      }
   t1 = thickx_(&h1);
   t2 = thickx_(&h2);
   
   *s0 = sum_log_dev_sq(a,b,c,np,h,t,rho);
   *npp = np;
   if ( np <= 0 )
      return 0.;
   else if ( h1 == h2 || t1 == t2 )
      return sum_log_dev_sq(a,b,c,np,h,t,rho);
   else if ( np <= 2 || ( np == 3 && (h1==h[0] || h2==h[1])) )
   {
      *cp = c = *cp;
      *bp = b = (t2-t1) / (exp(-h2/c)-exp(-h1/c));
      *ap = a = t1 - b*exp(-h1/c);
      return sum_log_dev_sq(a,b,c,np,h,t,rho);
   }
   
   c = *cp;
   b = (t2-t1) / (exp(-h2/c)-exp(-h1/c));
   a = t1 - b*exp(-h1/c);
   s = sum_log_dev_sq(a,b,c,np,h,t,rho);
   for ( iter=0, dc=2000.e2; iter<30; iter++ )
   {
      double a1, a2, a3, an, b1, b2, b3, bn, c1, c2, c3, cn, s1, s2, s3, sn;
      c2 = c;
      c1 = c - dc;
      if ( c1 <= 0. )
         c1 = c / 2.;
      c3 = c + dc;
      
      b1 = (t2-t1) / (exp(-h2/c1)-exp(-h1/c1));
      a1 = t1 - b1*exp(-h1/c1);
      b2 = (t2-t1) / (exp(-h2/c2)-exp(-h1/c2));
      a2 = t1 - b2*exp(-h1/c2);
      b3 = (t2-t1) / (exp(-h2/c3)-exp(-h1/c3));
      a3 = t1 - b3*exp(-h1/c3);
      s1 = sum_log_dev_sq(a1,b1,c1,np,h,t,rho);
      s2 = sum_log_dev_sq(a2,b2,c2,np,h,t,rho);
      s3 = sum_log_dev_sq(a3,b3,c3,np,h,t,rho);
      /* This works only for a roughly parabolic minimum */
      cn = ((c1*c1-c2*c2)*(s3-s2)-(c3*c3-c2*c2)*(s1-s2)) /
          (2.*(c1-c2)*(s3-s2)-2.*(c3-c2)*(s1-s2));
      if ( cn <= 0. )
         cn = 0.5*c;
      bn = (t2-t1) / (exp(-h2/cn)-exp(-h1/cn));
      an = t1 - bn*exp(-h1/cn);
      sn = sum_log_dev_sq(an,bn,cn,np,h,t,rho);
      if ( sn > s )
         break;
      if ( sn < s*1.00001 )
      {
         a = an;
         b = bn;
         c = cn;
         s = sn;
         break;
      }
      if ( cn < 3000.e2 )
      {
         cn = 3000.e2;
         bn = (t2-t1) / (exp(-h2/cn)-exp(-h1/cn));
         an = t1 - bn*exp(-h1/cn);
         dc = 1000.e2;
      }
      else if ( cn > c + 1.5*dc )
      {
         cn = c + 1.5*dc;
         bn = (t2-t1) / (exp(-h2/cn)-exp(-h1/cn));
         an = t1 - bn*exp(-h1/cn);
         dc *= 2.;
      }
      else if ( cn < c -1.5*dc )
      {
         cn = c -1.5*dc;
         bn = (t2-t1) / (exp(-h2/cn)-exp(-h1/cn));
         an = t1 - bn*exp(-h1/cn);
         dc *= 2.;
      }
      else if ( fabs(c-cn) < 0.25*dc )
         dc = 2.*fabs(c-cn) + 0.2*dc;
      else
         dc = 1.5*fabs(c-cn);
      a = an;
      b = bn;
      c = cn;
      s = sn;
   }
   *cp = c;
   *bp = b;
   *ap = a;
   return s;
}

/** Corresponding to CORSIKA built-in function THICK; 
    only used to show fit results. */

static double fn_thick(double h, int nl, double *hl, 
   double *a, double *b, double *c)
{
   int i;
   if ( h > top_of_atmosphere || nl < 2 )
      return 0;
   for ( i=0; i<nl-1; i++)
      if ( h < hl[i+1] )
         return a[i] + b[i] * exp((-h)/c[i]);
   return a[i] - h/c[i];
}

/** Corresponding to CORSIKA built-in function RHOF; 
    only used to show fit results. */

static double fn_rhof(double h, int nl, double *hl, 
   double *a, double *b, double *c)
{
   int i;
   if ( h > top_of_atmosphere || nl < 2 )
      return 0;
   for ( i=0; i<nl-1; i++)
      if ( h < hl[i+1] )
         return b[i]/c[i] * exp((-h)/c[i]);
   return 1./c[i];
}

#if 0
/** Corresponding to CORSIKA built-in function HEIGHT. */

static double fn_height(double t, int nl, double *hl, 
   double *a, double *b, double *c)
{
   int i;
   if ( t < 0. )
      return top_of_atmosphere;
   for ( i=0; i<nl-1; i++)
      if ( t > fn_thick(hl[i+1],nl,hl,a,b,c) )
         return c[i] * log(b[i] / (t-a[i]));
   return (a[i]-t) * c[i];
}
#endif

/* ---------------------------- atmfit_ ------------------------------ */
/**
 *  @short Fit the tabulated density profile for CORSIKA EGS part.
 *
 *  Fitting of the tabulated atmospheric density profile by
 *  piecewise exponential parts as used in CORSIKA.
 *  The fits are constrained by fixing the atmospheric thicknesses
 *  at the boundaries to the values obtained from the table.
 *  Note that not every atmospheric profile can be fitted well
 *  by the CORSIKA piecewise models (4*exponential + 1*constant
 *  density). In particular, the tropical model is known to 
 *  be a problem. Setting the boundary heights manually might help.
 *  The user is advised to check at least once that the fitted
 *  layers represent the tabulated atmosphere sufficiently well,
 *  at least at the altitudes most critical for the observations
 *  (usually at observation level and near shower maximum but
 *  depending on the user's emphasis, this may vary).
 *
 *  Fit all layers (except the uppermost) by exponentials and (if *nlp > 0)
 *  try to improve fits by adjusting layer boundaries.
 *  The uppermost layer has constant density up to the 'edge' of the atmosphere.
 *  
 *  This function may be called from CORSIKA.
 * 
 *  Parameters (all pointers since function is called from Fortran):
 *  @param   nlp    Number of layers (or negative of that if boundaries set manually)
 *  @param   hlay   Vector of layer (lower) boundaries.
 *  @param   aatm,batm,catm    Parameters as used in CORSIKA.
*/

void atmfit_(int *nlp, double *hlay, double *aatm, double *batm, double *catm)
{
   int il, np, k;
   double *a, *b, *c, *h, *s, *s0;
   double atmp[2], btmp[2], ctmp[2], htmp[3];
   int nl = (*nlp<0) ? -*nlp : *nlp;     /* Actual number of layers */
   double factmx = (*nlp<0) ? 1.0 : 1.4; /* Max. scale for boundary adjustment */
   int show_fit;

   if ( nl < 2 || atmosphere == 0 )
      return;

#ifdef DEBUG_ATM_FIT
   printf("\n Parameters of atmospheric layers as suggested by CORSIKA:\n");
   for ( il=0; il<nl; il++)
      printf(" Layer %d: %6.2f km < h < %6.2f km: a = %13.6g, b = %13.6g, c = %13.6g\n",
         il+1, hlay[il]/1e5,
         (il+1<(*nlp)) ? hlay[il+1]/1e5 : aatm[il]*catm[il]/1e5,
         aatm[il], batm[il], catm[il]);
   printf("\n");
#endif

   /* The lowest layer boundary must not be below the lowest table entry */
   /* because values can only be interpolated, not extrapolated. */
   if ( hlay[0] < bottom_of_atmosphere )
      hlay[0] = bottom_of_atmosphere;
   /* The default layers are known to be a rather bad choice for */
   /* the tropical atmosphere. Replace them with better starting values. */
   if ( *nlp > 0 && (atmosphere == 1 || (atmosphere >= 10 && atmosphere <30) ) )
   {
      hlay[1] = 9.25e5;
      hlay[2] = 19.0e5;
      hlay[3] = 37.5e5;
      hlay[4] = 105.0e5;
   }

   if ( (h = calloc((size_t)(nl+1),sizeof(double))) == NULL ||
        (s = calloc((size_t)(nl),sizeof(double))) == NULL ||
        (s0 = calloc((size_t)(nl),sizeof(double))) == NULL ||
        (a = calloc((size_t)(nl),sizeof(double))) == NULL ||
        (b = calloc((size_t)(nl),sizeof(double))) == NULL ||
        (c = calloc((size_t)(nl),sizeof(double))) == NULL )
   {
      return;
   }
   
   for ( il=0; il<nl; il++ )
   {
      a[il] = aatm[il];
      b[il] = batm[il];
      c[il] = catm[il];
      h[il] = hlay[il];
   }
   // h[nl] = a[nl-1]*c[nl-1];
   h[nl] = top_of_atmosphere;
   
   fflush(NULL);

#ifdef DEBUG_ATM_FIT
   printf("\n Parameters of atmospheric layers before the fit, after adjusting layers:\n");
   for ( il=0; il<nl; il++)
      printf(" Layer %d: %6.2f km < h < %6.2f km: a = %13.6g, b = %13.6g, c = %13.6g\n",
         il+1,h[il]/1e5,h[il+1]/1e5,a[il],b[il],c[il]);
   printf("\n");

   /* Functions and variables for cut-and-paste into gnuplot or C: */
   
   puts("thicka(h)=(h<h2a?a1a+b1a*exp(-h*1e5/c1a):(h<h3a?a2a+b2a*exp(-h*1e5/c2a):(h<h4a?a3a+b3a*exp(-h*1e5/c3a):(h<h5a?a4a+b4a*exp(-h*1e5/c4a):(h<h6a?a5a-h*1e5/c5a:0)))))");
   puts("rhoa(h)=(h<h2a?b1a/c1a*exp(-h*1e5/c1a):(h<h3a?b2a/c2a*exp(-h*1e5/c2a):(h<h4a?b3a/c3a*exp(-h*1e5/c3a):(h<h5a?b4a/c4a*exp(-h*1e5/c4a):(h<h6a?1./c5a:0)))))");
   for ( il=0; il<=nl;il++)
     printf("h%da=%6.2f;",il+1,h[il]/1.e5);
   printf("\n");
   for ( il=0; il<nl;il++)
     printf("a%da=%13.6g;",il+1,a[il]);
   printf("\n");
   for ( il=0; il<nl;il++)
     printf("b%da=%13.6g;",il+1,b[il]);
   printf("\n");
   for ( il=0; il<nl;il++)
     printf("c%da=%13.6g;",il+1,c[il]);
   printf("\n");
#endif

   for ( il=0; il<nl-1; il++ )
      s[il] = atm_exp_fit(h[il],h[il+1],&a[il],&b[il],&c[il],&s0[il],&np);

   c[nl-1] = 2./rhofx_(&h[nl-1]);
   a[nl-1] = thickx_(&h[nl-1]) + h[nl-1]/c[nl-1];
   b[nl-1] = 1.;
   h[nl] = a[nl-1]*c[nl-1];

   if ( factmx > 1.0 ) /* Try to improve by adjusting boundaries */
   {
#ifdef DEBUG_ATM_FIT
      printf("\n Intermediate results of atmosphere fit:\n");
      for ( il=0; il<nl; il++)
         printf(" Layer %d: %6.2f km < h < %6.2f km: a = %13.6g, b = %13.6g, c = %13.6g, s = %10.4e -> %10.4e\n",
            il+1,h[il]/1e5,h[il+1]/1e5,a[il],b[il],c[il],s0[il],s[il]);
      printf("\n");
      for ( il=0; il<=nl;il++)
        printf("h%db=%6.2f;",il+1,h[il]/1.e5);
      printf("\n");
      for ( il=0; il<nl;il++)
        printf("a%db=%13.6g;",il+1,a[il]);
      printf("\n");
      for ( il=0; il<nl;il++)
        printf("b%db=%13.6g;",il+1,b[il]);
      printf("\n");
      for ( il=0; il<nl;il++)
        printf("c%db=%13.6g;",il+1,c[il]);
      printf("\n");
#endif

      for (k=0; k<2; ++k)
      {
       for ( il=nl-2; il>=0; il-- )
       {
         double smin[2], stmp[2], s0tmp[2], hvar;
         int nptmp[2], npmin[2];
         atmp[0] = a[il]; atmp[1] = a[il+1];
         btmp[0] = b[il]; btmp[1] = b[il+1];
         ctmp[0] = c[il]; ctmp[1] = c[il+1];
         htmp[0] = h[il]; htmp[1] = h[il+1]; htmp[2] = h[il+2];
         smin[0] = atm_exp_fit(htmp[0],htmp[1],&atmp[0],&btmp[0],&ctmp[0],&s0tmp[0],&npmin[0]);
         if ( il < nl-2 )
            smin[1] = atm_exp_fit(htmp[1],htmp[2],&atmp[1],&btmp[1],&ctmp[1],&s0tmp[1],&npmin[1]);
         else
         {
            smin[1] = 0.;
            npmin[1] = 4;
         }
#ifdef DEBUG_ATM_FIT
         printf("In layers %d/%d: smin=%f/%f, npmin=%d/%d\n",
            il, il+1, smin[0], smin[1], npmin[0], npmin[1]);
#endif
         if ( npmin[0] <= 3 || npmin[1] <= 3 )
            continue;

#ifdef DEBUG_ATM_FIT
         printf("In layers %d/%d try all hvar from %f km to %f km with factmx=%f\n",
            il, il+1, (htmp[0]+3.e5)/1e5, (htmp[2]-3.e5)/1e5, factmx);
#endif
         for ( hvar=(il<nl-2?htmp[0]+3.e5:htmp[1]); hvar<=htmp[2]-3.e5; hvar+=0.25e5 )
         {
            if ( hvar < h[il+1]/factmx || hvar > h[il+1]*factmx || hvar > h[il+1]+10.e5 )
               continue;
            atmp[0] = a[il]; atmp[1] = a[il+1];
            btmp[0] = b[il]; btmp[1] = b[il+1];
            ctmp[0] = c[il]; ctmp[1] = c[il+1];
            stmp[0] = atm_exp_fit(htmp[0],hvar,&atmp[0],&btmp[0],&ctmp[0],&s0tmp[0],&nptmp[0]);
            if ( il < nl-2 )
               stmp[1] = atm_exp_fit(hvar,htmp[2],&atmp[1],&btmp[1],&ctmp[1],&s0tmp[1],&nptmp[1]);
            else
            {
               stmp[1] = 0.;
               nptmp[1] = 4;
            }
#ifdef DEBUG_ATM_FIT
            printf("Try hvar=%f: stmp=%f/%f, nptmp=%d/%d\n", hvar/1e5, 
               stmp[0], stmp[1], nptmp[0], nptmp[1]);
            printf("             smin=%f/%f, npmin=%d/%d\n",             
               smin[0], smin[1], npmin[0], npmin[1]);
#endif
            if ( nptmp[0] >= 3 && nptmp[1] >= 3 &&
                 stmp[0]/(nptmp[0]-2.9)+stmp[1]/(nptmp[1]-2.5) <
                 smin[0]/(npmin[0]-2.9)+smin[1]/(npmin[1]-2.5) )
            {
#ifdef DEBUG_ATM_FIT
               printf("Looks like we got an improvement\n");
#endif
               htmp[1] = hvar;
               smin[0] = stmp[0];
               smin[1] = stmp[1];
               npmin[0] = nptmp[0];
               npmin[1] = nptmp[1];
            }
         }
         h[il+1] = htmp[1];
       }
      }

      for ( il=0; il<nl-1; il++ )
         s[il] = atm_exp_fit(h[il],h[il+1],&a[il],&b[il],&c[il],&s0[il],&np);

      c[nl-1] = 2./rhofx_(&h[nl-1]);
      a[nl-1] = thickx_(&h[nl-1]) + h[nl-1]/c[nl-1];
      b[nl-1] = 1.;
      h[nl] = a[nl-1]*c[nl-1];
   }
#ifdef DEBUG_ATM_FIT
   else
      printf(" No layer adjustment enabled\n");
#endif

   printf("\n Results of the atmosphere fit:\n");
   for ( il=0; il<nl; il++)
      printf(" Layer %d: %6.2f km < h < %6.2f km: a =%13.6g, b =%13.6g, c =%13.6g\n",
         il+1,h[il]/1e5,h[il+1]/1e5,a[il],b[il],c[il]);
   printf("\n");
#ifdef DEBUG_ATM_FIT
   puts("thickb(h)=(h<h2b?a1b+b1b*exp(-h*1e5/c1b):(h<h3b?a2b+b2b*exp(-h*1e5/c2b):(h<h4b?a3b+b3b*exp(-h*1e5/c3b):(h<h5b?a4b+b4b*exp(-h*1e5/c4b):(h<h6b?a5b-h*1e5/c5b:0)))))");
   puts("rhob(h)=(h<h2b?b1b/c1b*exp(-h*1e5/c1b):(h<h3b?b2b/c2b*exp(-h*1e5/c2b):(h<h4b?b3b/c3b*exp(-h*1e5/c3b):(h<h5b?b4b/c4b*exp(-h*1e5/c4b):(h<h6b?1./c5b:0)))))");
   for ( il=0; il<=nl;il++)
     printf("h%db=%6.2f;",il+1,h[il]/1.e5);
   printf("\n");
   for ( il=0; il<nl;il++)
     printf("a%db=%13.6g;",il+1,a[il]);
   printf("\n");
   for ( il=0; il<nl;il++)
     printf("b%db=%13.6g;",il+1,b[il]);
   printf("\n");
   for ( il=0; il<nl;il++)
     printf("c%db=%13.6g;",il+1,c[il]);
   printf("\n");
#endif

#ifdef RESPECT_EGS_TOP_OF_ATMOSPHERE
   if ( h[nl] > 113e5 ) /* There was a fixed limit of 113 km in EGS part */
   {
      fflush(NULL);
      if ( h[nl-1] >= 112.99e5 )
      {
         fprintf(stderr,
         "\n Atmospheric boundaries cannot satisfy requirements for CORSIKA\n");
         fprintf(stderr," EGS part (upper end at %6.2f instead of 113.00 km)\n\n",
            h[nl]/1e5);
      }
      else
      {
         double rho;
         fprintf(stderr,"\n Upper atmospheric boundary forced from %6.2 to 113.00 km\n",
             h[nl]/1e5);
         fprintf(stderr," to satisfy requirements for CORSIKA EGS part\n");
         h[nl] = 113e5;
         rho = thickx_(&h[nl-1]) / (h[nl]-h[nl-1]);
         c[nl-1] = 1./rho;
         a[nl-1] = thickx_(&h[nl-1]) + h[nl-1]/c[nl-1];
      }
   }
#endif

   for ( il=0; il<nl; il++ )
   {
      aatm[il] = a[il];
      batm[il] = b[il];
      catm[il] = c[il];
      hlay[il] = h[il];
   }

   /* Finally, we amend our table such that the density and thickness */
   /* at the level corresponding to the top of the atmosphere are zero */
   /* and any further levels are ignored. */
   /* Note that results from the interpolation and the CORSIKA formulae */
   /* will still not agree in the upper linear density gradient component. */
   top_of_atmosphere = h[nl];
   for (il=0; il<num_prof; il++)
      if ( p_alt[il] >= top_of_atmosphere )
      {
         p_alt[il] = top_of_atmosphere;
         p_rho[il] = p_rho[il-1];
         p_log_rho[il] = p_log_rho[il-1];
         p_log_thick[il] = -100;
         p_log_n1[il] = p_log_n1[il-1];
         num_prof = il+1;
         break;
      }

   /* An environment variable can be used to show fit compared to table  */
   if ( getenv("SHOWFIT") != NULL )
      show_fit = atoi(getenv("SHOWFIT"));
   else
      show_fit = 1;

   if ( show_fit )
   {
      int ip;
      printf("\n Altitude [km]    rho(table)     rho(fit)       thick(table)  thick(fit)\n");
      for (ip=0; ip<num_prof; ip++ )
      {
      	 if ( p_alt[ip] >= bottom_of_atmosphere && 
	      p_alt[ip] <= top_of_atmosphere )
	 {
	    printf("      %5.1f     %12.5e  %12.5e     %12.5e  %12.5e\n",
	       p_alt[ip]/1e5, p_rho[ip], /* rhofx_(&p_alt[ip]) */ 
                  fn_rhof(p_alt[ip],nl,h,a,b,c),
	       exp(p_log_thick[ip]), /* thickx_(&p_alt[ip]) */ 
                  fn_thick(p_alt[ip],nl,h,a,b,c));
	 }
      }
      printf("\n");
      if ( show_fit == 999 )
      {
      	 printf("Terminating now as requested by the SHOWFIT environment variable.\n");
	 exit(0);
      }
   }

   free(h); free(s); free(s0); free(a); free(b); free(c);

   fflush(stdout);

   return;
}


#ifdef TEST_ATMO

double rhof_(double *h) { exit(1); return 0.; }
double heigh_(double *t) { exit(2); return 0.; }
double thick_(double *h) { exit(3); return 0.; }

int main (int argc, char **argv)
{
   int iatmo = 0;
   double olev = 1800e2;
   double h, t, h2;
   if ( argc > 1 )
      iatmo = atoi(argv[1]);
   if ( argc > 2 )
      olev = atof(argv[2]);
   atmset_(&iatmo,&olev);
   
   for ( h=olev; h<30e5; h+=0.001e5 )
   {
      t  = thickx_(&h);
      h2 = heighx_(&t);
      printf("%f %f %f\n",h,t,h2-h);
   }
   return 0;
}

#endif
