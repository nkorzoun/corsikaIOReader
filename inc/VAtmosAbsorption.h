//! VAtmosAbsorption atmospheric absorption of Cherenkov photons
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

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "TMath.h"
#include "TRandom3.h"

// #include "TCanvas.h"  // only needed for visualisation
// #include "TGraph.h"   // only needed for visualisation
// #include "TH1D.h" // only needed for visualisation
// #include "TLegend.h" // only needed for visualisation
// #include "TStyle.h"  // only needed for visualisation

#ifndef VATMOSABSORPTION_H
#define VATMOSABSORPTION_H

using namespace std;

class VAtmosAbsorption
{
    private:
        string fModel;                   //!< used model (CORSIKA,kascade)
        string fSourceFile;
        
        double fObservationLevel;        //!< observation height in [m]
        
        double fminWave;                 //!< minimum wavelength in [nm]
        double fmaxWave;                 //!< maximum wavelength in [nm]
        TRandom3* fRandom;                //!< random generator
        
        map< int, vector<double> > fCoeff;   //!<   vector with atmospheric data (CORSIKA)
        map< int, double > fCoeffObs;        //!<   vector with atmospheric data at observation level (CORSIKA)
        
        vector< vector<double> > extint;     //!<   vector with atmospheric data (kascade)
        
        double getLinearInterpolate( double x, double x0, double x1, double y0, double y1 );  //!< linear interpolation
        void   readCorsikaAtmabs();                  //!< read CORSIKA atmospheric extinction file (atmabs.dat)
        void   read_extint( int );                        //!< read kascade atmospheric extinction file (kextint.dat)
        void   read_extint_F2( int );                        //!< read kascade atmospheric extinction file (kextint.dat)
        void   read_extint_M5( );                         //!< read henrikes atmospheric extinction file (modtran 5)
        
    public:
        VAtmosAbsorption( string, int, string iSourceFile = "" );
        ~VAtmosAbsorption() {}
        void setWavelengthintervall( double iminwavelength, double imaxwavelength );   //!< in [nm]
        void setObservationlevel( double iobslevel );  //!< in [m]
        double probAtmAbsorbed( double wavelength, double emissionheigth, double emissionangle );   //!< calculates survival probability for photon
        double probAtmAbsorbed( double wavelength, double emissionheigth, double emissionangle, double& obsdepth );   //!< calculates survival probability for photon
        double getWavelength( double emissionheigth, double emissionangle );  //!< get random wavelength
};

#endif


