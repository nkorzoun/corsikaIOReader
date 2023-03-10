//!  VIOHistosgrams  histogram filling class
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

#ifndef VHISTO_H
#define VHISTO_H

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <utility>//NK

#include "mc_tel.h"
#include "sim_cors.h"

using namespace std;

class VIOHistograms
{
    private:
        TList* hisList;
        
        TH1D* hBunch;
        TH1D* hT0;
        TH1D* hZem;
        
        TH2D* hGXY;
        TH2D* hGCXCY; //NK
        //TH2D* hGZeAz; NK
        TH1D* hGLambda;
        TH1D* hGProb;
        TH1D* hGZem;
        TH2D* hSXY;
        TH2D* hSCXCY; //NK
        //TH2D* hSZeAz; NK
        TH1D* hSLambda;
        TH1D* hSProb;
        TH1D* hSZem;
        
        TClonesArray* hCXYZ;
        
        TFile* fout;
        TTree* fTree;
        
        vector< double > fXYZlevelsHeight;
        vector< double > fXYZlevelsThickness;

        
        int nevent;
        double degrad;
        
        // event block
        
        int eventNumber;
        int particleID;
        double energy;
        double startAlt;
        double firstInt;
        double ze;
        double az;
        int date;
        int runNumber;
        double corsVersion;
        double obsLevel;
        double eslope;
        double e0min;
        double e0max;
        double cutHad;
        double cutMuon;
        double cutEM;
        double cutPhot;
        double magX;
        double magZ;
        int bCher;
        double cherBunch;
        double cherLambdaMin;
        double cherLambdaMax;
        int    telNumber;
        int    arrayNumber;
        int    coreNumber;
        double telXpos[MAX_TEL];
        double telYpos[MAX_TEL];
        double telZpos[MAX_TEL];
        double telR[MAX_TEL];
        double xCore;
        double yCore;
        double rCore;
        double xmax;
        double emax;
        double cmax;
        double hmax;
        double NCp[MAX_TEL];
        
        std::vector <double > CX;  //NK
        std::vector <double > CY;  //NK
        std::vector <int > telID; //NK
        //TH2D* CXCY[MAX_TEL];

        double toff;
        bool bShort;
        
        bool bCORSIKA_coordinates;   // fill photons in corsika coordinates
        
        bool bSmallFile; // reduce output
        bool bMuon;      // adjust histograms for muon input
        
        double redang( double );
        void transformCoord( double&, double&, double& );
        
    public:
        VIOHistograms();
        ~VIOHistograms() {}
        void init( string, bool );
        void newEvent( float*, telescope_array, int );
        void initXYZhistograms();
        void fillBunch( bunch, double );
        void fillGenerated( bunch, double );
        void fillNPhotons( int iTel, double iphotons );
        void fillSurvived( bunch, double, float*, int );
        void setCORSIKAcoordinates()
        {
            bCORSIKA_coordinates = true;
        }
        void set_small_file()
        {
            bSmallFile = true;
        }
        void setMuonSettings()
        {
            bMuon = true;
        }
        void setXYZlevelsHeight( vector<double> iL )
        {
            fXYZlevelsHeight = iL;
        }
        void setXYZlevelsThickness( double istart, double idiff );
        void terminate();
};

#endif
