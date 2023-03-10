//! VCORSIKARunheader
/*
=============================================================================
    corsikaIOreader is a tool to read CORSIKA eventio files
    Copyright (C) 2004, 2013, 2019 Gernot Maier and Henrike Fleischhack

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
=============================================================================
*/

#ifndef VCORSIKARunheader_H
#define VCORSIKARunheader_H

#include "TNamed.h"

#include <bitset>
#include <iostream>
#include <ostream>

using namespace std;

class VCORSIKARunheader : public TNamed
{
    private:
    
        void reset();
        
    public:
    
        unsigned int runnumber;                    // run number
        unsigned int production_date;              // date of simulations
        float corsika_version;                     // corsika version
        
        float observation_level_m;                 // observation level [m]
        
        unsigned int particleID;                   // CORSIKA particle ID
        
        float startingaltitude_gcm2;               // starting altitude [g/cm2]
        float tstart;                              // tstart
        
        float energy_slope;                        // slope of energy spectrum
        float energy_min_GeV;                      // lower limit in energy range [GeV]
        float energy_max_GeV;                      // upper limit in energy range [GeV]
        
        float zenith_min_deg;
        float zenith_max_deg;
        float azimuth_min_deg;
        float azimuth_max_deg;
        
        float viewcone_min_deg;
        float viewcone_max_deg;
        
        unsigned int nscatt;                       // number of uses of each Cherenkov event
        float xscatt_m;                            // scatter range in x direction
        float yscatt_m;                            // scatter range in y direction
        
        /* cherenkov flag
           1 CERENKOV option compiled in
           2 IACT option compiled in
           3 CEFFIC option compiled in
           4 ATMEXT option compiled in
           5 ATMEXT option used with refraction enabled
           6 VOLUMEDET option compiled in
           7 CURVED option compiled in
           9 SLANT option compiled in
           11-21 table number for external atmosphere table */
        unsigned long int cherenkov_flag;
        
        float cherenkov_bunchsize;                 // Cherenkov photon bunch size
        float cherenkov_bandwidth_min_nm;          // Cherenkov bandwidth lower end
        float cherenkov_bandwidth_max_nm;          // Cherenkov bandwidth upper end
        
        float geomagneticfield_arrang_deg;         // magnetic field rotation
        float geomagneticfield_x_muT;              // x component of Earths magnetic field in muT
        float geomagneticfield_z_muT;              // z component of Earths magnetic field in muT
        
        unsigned int hadronic_model_low;           // low-energy hadr. model flag (1.=GHEISHA, 2.=UrQMD, 3.=FLUKA
        unsigned int hadronic_model_high;          // 0.=HDPM,1.=VENUS, 2.=SIBYLL,3.=QGSJET, 4.=DPMJET, 5.=NEXUS, 6.=EPOS)
        float hadronic_model_transition_energy_GeV; // transition energy high-energy/low-energy model (in GeV)
        
        VCORSIKARunheader();
        ~VCORSIKARunheader() {}
        void printHeader( std::ostream& );
        // void printHeader( ofstream );
        
        ClassDef( VCORSIKARunheader, 1 );
};

#endif

