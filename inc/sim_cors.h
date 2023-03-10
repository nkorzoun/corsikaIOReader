#ifndef SIM_CORS_H
#define SIM_CORS_H

/** Basic parameters of the CORSIKA run */

#define MAX_TEL 1200            /**< The largest no. of telescopes/array. */
#define MAX_ARRAY 100         /**< The largest no. of arrays to be handled */

struct mc_run
{
   int data;
   double version;
   double runnumber;
   double height;          /**< Height of observation level [m] */
   double e_min;           /**< Lower limit of simulated energies [TeV] */
   double e_max;           /**< Upper limit of simulated energies [TeV] */
   double slope;           /**< Spectral index of power-law spectrum */
   double radius;          /**< Radius within which cores are thrown at random. [m] */
   int num_arrays;         /**< Number of arrays simulated. */
   double theta_min;       /**< Lower limit of zenith angle [degrees] */
   double theta_max;       /**< Upper limit of zenith angle [degrees] */
   double phi_min;         /**< Lower limit of azimuth angle [degrees] */
   double phi_max;         /**< Upper limit of azimuth angle [degrees] */
   double wlen_min;        /**< Lower limit of Cherenkov wavelength range [nm] */
   double wlen_max;        /**< Upper limit of Cherenkov wavelength range [nm] */
   double bunchsize;       /**< Cherenkov bunch size. */
};

/** Basic parameters of a simulated shower */

struct simulated_shower_parameters
{
   double energy;                /**< Shower energy [TeV] */
   double azimuth;               /**< Shower direction azimuth [deg] */
   double altitude;              /**< Shower direction altitude above horizon */
   double xcore, ycore, zcore;   /**< Shower core position [m] */
   double core_dist_3d;          /**< Distance of core from reference point */
   double tel_core_dist_3d[MAX_TEL];  //*< Offset of telescopes from shower axis */
   int particle;                 /**< Primary particle type [CORSIKA code] */
   double xmax;            /**< Depth of shower maximum from all particles [g/cm**2] */
   double emax;            /**< Depth of shower maximum from positrons and electrons */
   double cmax;            /**< Depth of maximum of Cherenkov light emission [g/cm**2] */
   double hmax;            /**< Height of shower maximum (from xmax above) [m] a.s.l. */
   double firstint;	//height of first interaction in m
   int shower_id;	//corsika event id
};

/**  Description of telescope position, array offets and shower parameters. */

struct telescope_array
{
   int ntel;               /**< Number of telescopes simulated per array */
   int max_tel;            /**< Maximum number of telescopes acceptable (MAX_TEL) */
   int narray;             /**< Number of arrays with random shifts per shower */
   double refpos[3];       /**< Reference position with respect to obs. level [cm] */
   double obs_height;      /**< Height of observation level [cm] */
   double xtel[MAX_TEL];   /**< X positions of telescopes ([cm] -> north) */
   double ytel[MAX_TEL];   /**< Y positions of telescopes ([cm] -> west) */
   double ztel[MAX_TEL];   /**< Z positions of telescopes ([cm] -> up) */
   double rtel[MAX_TEL];   /**< Radius of spheres enclosing telescopes [cm] */
   double toff;            /**< Time offset from first interaction to the moment */
                           /**< when the extrapolated primary flying with the vacuum */
                           /**< speed of light would be at the observation level. */
   double xoff[MAX_ARRAY]; /**< X offsets of the randomly shifted arrays [cm] */
   double yoff[MAX_ARRAY]; /**< Y offsets of the randomly shifted arrays [cm] */
   double azimuth;         /**< Nominal azimuth angle of telescope system [deg]. */
   double altitude;        /**< Nominal altitude angle of telescope system [deg]. */
   double source_azimuth;  /**< Azimuth of assumed source. */
   double source_altitude; /**< Altitude of assumed source. */
   struct simulated_shower_parameters shower_sim;
   struct mc_run mc_run;
   /* ... */
};

#endif
