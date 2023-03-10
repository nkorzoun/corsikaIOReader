//! VGrisu  CORSIKA to GrIsu writer
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

#ifndef VGRISU_H
#define VGRISU_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

#include "VCORSIKARunheader.h"

#include "mc_tel.h"
#include "sim_cors.h"
#include "atmo.h"
#include "fileopen.h"

using namespace std;

class VGrisu
{
    private:
        bool bSTDOUT;                        //!< write output to stdout
        ofstream of_file;                    //!< output file
        map<int, int> particles;             //!< particle ID transformation matrix: first: CORSIKA ID, second: kascade ID
        float degrad;                        //!< convertion deg->rad
        int primID;                          //!< primary particle ID
        double xoff;                         //!< offset in x-coordinate for scattered core
        double yoff;                         //!< offset in y-coordinate for scattered core
        double qeff;                         //!< global quantum efficiency
        double observation_height;           //!< observation height [m]
        
        int atm_id;			  //!< corsika atmprof number (for depth calculation)
        
        
        string fVersion;
        
        void transformCoord( float&, float&, float& );     //!< transform from CORSIKA to GrIsu coordinates
        void makeParticleMap();              //!< make map with  particle ID transformation matrix
        float redang( float );             //! reduce large angle to intervall 0, 2*pi
        
    public:
        VGrisu( string fVersion = "", int id = -1 );
        ~VGrisu() {}
        void setOutputfile( string );       //!< create grisu readable output file
        void setObservationHeight( double ih )
        {
            observation_height = ih;    //!< set observation height
        }
        void setQeff( double iq )
        {
            qeff = iq;    //!< set global quantum efficiency
        }
        void setQueff( double iq )
        {
            qeff = iq;    //!< set global quantum efficiency
        }
        void writeRunHeader( float*, VCORSIKARunheader* );      //!<  write some information about CORSIKA run intot the runheader
        void writeEvent( telescope_array, bool );  //!< write MC information
        void writePhotons( bunch, int );         //!< write next photon to grisu file ("P" line)
};

#endif
