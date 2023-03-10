/* ============================================================================

   Copyright (C) 2001, 2005, 2009  Konrad Bernloehr

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
 *  @file atmo.h
 *  @short Use of tabulated atmospheric profiles and atmospheric refraction.
 *
 *  @author  Konrad Bernloehr
 *  @date    @verbatim CVS $Date: 2012/12/03 08:23:07 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.2 $ @endverbatim
 *
 */
/* ================================================================== */

#ifndef ATMO_H__LOADED

#define ATMO_H__LOADED 1

/* The CORSIKA version against which this software should match. */
/* If your CORSIKA is somewhat newer than 5.901 there is probably */
/* no reason to worry; incompatible changes should not happen */
/* all too often. */
#ifndef CORSIKA_VERSION
# define CORSIKA_VERSION 6000
#endif

#if (CORSIKA_VERSION < 5901)
typedef float cors_real_now_t;
#else
typedef double cors_real_now_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Function prototypes for functions implemented in this file */

/* FORTRAN called functions (beware changes of parameter types !!) */
void atmset_(int *iatmo, double *obslev);
double rhofx_(double *height);
double thickx_(double *height);
double refidx_(double *height);
double heighx_(double *thick);
void raybnd_(double *zem, cors_real_now_t *u, cors_real_now_t *v, double *w, 
   cors_real_now_t *dx, cors_real_now_t *dy, cors_real_now_t *dt);
void atmfit_(int *nlp, double *hlay, double *aatm, double *batm, double *catm);

/* C called functions (parameter types are always checked) */
double rpol(double *x, double *y, int n, double xp);

/* FORTRAN functions called from C */
/// The CORSIKA built-in density lookup function.
double rhof_(double *height);
/// The CORSIKA built-in function for vertical atmospheric thickness (overburden).
double thick_(double *height);
/// The CORSIKA built-in function for the height as a function of overburden.
double heigh_(double *thick);

#ifdef __cplusplus
}
#endif

#endif
