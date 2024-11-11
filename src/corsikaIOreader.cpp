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

    This code is based on a skeleton program provided by Konrad Bernloehr (sim_skeleton.c)
    as part of the CORSIKA IACT package

    see README file for documentation

    options: try corsikaIOreader -help


*/

#include "initial.h"      /* This file includes others as required. */
#include "io_basic.h"     /* This file includes others as required. */
#include "mc_tel.h"
#include "atmo.h"
#include "sim_cors.h"

#include <cmath>
#include <bitset>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "VAtmosAbsorption.h"        // atmospheric extinction class
#include "VCORSIKARunheader.h"
#include "VIOHistograms.h"           // histogramming class (only needed for test histograms)
#include "VGrisu.h"                  // writing of grisu format

#include "TRandom3.h"                 // if you don't like root -> use your own random generator
// + delete all VIOHistograms lines

#define MAX_BUNCHES 50000000   // (GM) why this limitation? (original 50000)

static double airlightspeed = 29.9792458 / 1.0002256; /* [cm/ns] at H=2200 m */

/*! Refraction index of air as a function of height in km (0km<=h<=8km) */
#define Nair(hkm) (1.+0.0002814*exp(-0.0947982*(hkm)-0.00134614*(hkm)*(hkm)))

#ifndef Nint
#define Nint(x) ((x)>0?(int)((x)+0.5):(int)((x)-0.5))
#endif

using namespace std;

bool bDebug = false;

string fVersion = "corsikaIOreader (v 2.0.0)";


struct linked_string corsika_inputs;

/* ========================== Utility functions ======================= */

/* ------------------- line_point_distance --------------------- */
/**
 *  Distance between a straight line and a point in space
 *
 *  @param  x1,y1,z1  reference point on the line
 *  @param  cx,cy,cz  direction cosines of the line
 *  @param  x,y,z     point in space
 *
 *  @return distance
 *
*/

double line_point_distance( double x1, double y1, double z1,
                            double cx, double cy, double cz,
                            double x, double y, double z )
{
    double a, a1, a2, a3, b;
    
    a1 = ( y - y1 ) * cz - ( z - z1 ) * cy;
    a2 = ( z - z1 ) * cx - ( x - x1 ) * cz;
    a3 = ( x - x1 ) * cy - ( y - y1 ) * cx;
    a  = a1 * a1 + a2 * a2 + a3 * a3;
    b = cx * cx + cy * cy + cz * cz;
    if( a < 0. || b <= 0. )
    {
        return -1;
    }
    return sqrt( a / b );
}

/*
     read telescopes matrix to convert corsika telescope numbers to a subset for grisudet

     (example: corsika simulated for 5 telescopes, but grisudet should be for 4)

*/
vector< int > readTelescopeMatrix( string iCFGFile, int ntel, double* xtel, double* ytel )
{
    vector< int > m( ntel, 0 );
    for( int i = 0; i < ntel; i++ )
    {
        m[i] = i;
    }
    
    if( iCFGFile.size() == 0 )
    {
        return m;
    }
    for( int i = 0; i < ntel; i++ )
    {
        m[i] = -1;
    }
    
    ifstream is;
    is.open( iCFGFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "readTelescopeMatrix error opening grisudet cfg file " << iCFGFile << endl;
        cout << "...exiting" << endl;
        exit( -1 );
    }
    unsigned int iTelID;
    double x = 0.;
    double y = 0.;
    string iTemp;
    string is_line;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            is_stream >> iTemp;
            if( iTemp != "*" )
            {
                continue;
            }
            is_stream >> iTemp;
            if( iTemp != "TLLOC" )
            {
                continue;
            }
            
            is_stream >> iTemp;
            iTelID = atoi( iTemp.c_str() );
            is_stream >> iTemp;
            x = atof( iTemp.c_str() );
            is_stream >> iTemp;
            y = atof( iTemp.c_str() );
            
            // check which telescope position is consistent here
            // observe: corsika in cm
            // observe: corsika and grisudet have different coordinate systems
            // observe: ignore z coordinate
            for( int i = 0; i < ntel; i++ )
            {
                if( sqrt( ( x + ytel[i] / 1.e2 ) * ( x + ytel[i] / 1.e2 ) + ( y - xtel[i] / 1.e2 ) * ( y - xtel[i] / 1.e2 ) ) < 0.5 )
                {
                    m[i] = iTelID - 1;
                }
            }
        }
    }
    is.close();
    
    return m;
}

/*!
    get XYZ levels for histograms in [m]

    end at 1005 gms

*/
vector< double > getXYZlevels( double istart, double idiff )
{
    int iatmo = 6;
    double obs_height = 0.;
    double ih;
    
    try
    {
        atmset_( &iatmo, &obs_height );
    }
    catch(...)
    {
        cout << "error initialising atmospheres" << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    
    vector< double > iL;
    if( idiff <= 0. )
    {
        return iL;
    }
    
    int nlevel = ( int )( ( 1005. - istart ) / idiff );
    
    if( bDebug )
    {
        cout << "filling XYZ levels: " << endl;
    }
    for( int i = 0; i <= nlevel; i++ )
    {
        ih = istart + i * idiff;
        iL.push_back( heighx_( &ih ) / 100. );
        if( bDebug )
        {
            cout << "\t XYZ level " << i << "\t" << ih << " [g/cm2] " << iL.back() << " [m]" << endl;
        }
    }
    return iL;
}


/**********************************************************************************

**********************************************************************************/
int main( int argc, char** argv )
{
    IO_BUFFER* iobuf;
    IO_ITEM_HEADER item_header, sub_item_header, block_header;
    real runh[273], rune[273], evth[273], evte[273];
    
    static struct telescope_array array;
    
    bool bGRISU = false;    // if true, nothing is printed to stdout except from VGrisu (and error messages)
    bool bstdout = false;
    bool bHisto = false;    // if true, tree and histograms are filled
    bool bPrintHeaders = false;
    int nbunches;
    int itc, iarray, jarray, ibunch;
    double lambda;
    double photons;
    double wl_bunch, airmass, cx, cy, cz, prob;
    double tel_dist, tel_delay;
    FILE* data_file;
    // static double elow;
    static struct bunch bunches[MAX_BUNCHES];
    static int particle_type;
    static double primary_energy;
    static double wl_lower_limit, wl_upper_limit;
    int run = -1;
    // int event = -1;
    double alt = 0.;
    double az = 0.;
    static int have_atm_profile = 0;
    // double toffset = 0.;
    double corstime;
    bunch Chphoton;
    string fCorsikaIO = "";                              // corsika io file
    string fAtmosModel = "noExtinction";
    string fAtmosFile  = "data/us76.50km.ext";
    double queff = 1.;
    bitset<32> EVTH76;
    bool bCEFFICWARNING = true;
    bool bPrintMoreInfo = false;
    
    int readNevent = 0;               // event counter
    int nevents = -1;                 // number of events to be read
    int readNarray = 0;               // array per event counter
    int narray = -1;                  // number of arrays per event to be read
    int nTel = -1;                    // telescope number to be read (-1: all in file)
    int fSeed = 0;
    int atmid = -1;
    
    // this is the grisu format output class
    vector< VGrisu* > fGrisu;
    string fGrisuOutputFile = "";
    // histogramming class (only filled with switch -histo/shorthisto)
    VIOHistograms* fHisto = new VIOHistograms();
    // matrix of telescope numbering: needed if telescope numbers in grisudet and corsika disagree
    // (example corsika: 400 telescopes, grisudet: 49 telescopes)
    // grisu cfg file
    string fGrisuConfigurationFile;
    vector< int > fTelescopeMatrix;
    // corsika run header class
    VCORSIKARunheader* fRunHeader = new VCORSIKARunheader();
    
    // reading of command line arguments
    int i = 0;
    while( i++ < argc )
    {
        string iTemp = argv[i - 1];
        string iTemp2 = "";
        if( i < argc )
        {
            iTemp2 = argv[i];
        }
        // some command line arguments are wrong or need of help
        if( argc == 1 || iTemp.find( "-help" ) < iTemp.size() )
        {
            cout << endl;
            cout << fVersion << endl;
            cout << "=======================" << endl << endl;
            cout << "corsikaIOreader is a tool to read CORSIKA eventio files and" << endl;
            cout << "  - convert photon information into grisu format" << endl;
            cout << "  - fill root histograms with particle and Cherenkov photon information" << endl;
            cout << endl;
            cout << "Based on skeleton.c in the CORSIKA IACT package. Use of Konrad Bernloehr's IACT routines, which are part of the CORSIKA distribution.";
            cout << endl << endl;
            cout << "Command line options: " << endl << endl;
            cout << "\t -cors  IOFILENAME     CORSIKA io-style particle file" << endl;
            cout << "\t -ioread FILENAME      write eventio file contents in Grisu style into FILENAME (stdout if output to stdout is wanted)" << endl;
            cout << "\t -histo FILE.root      fill eventio file contents into histograms" << endl;
            cout << "\t -xyz FILE.root        fill  eventio file contents into histograms (with photon xy positions for different heights)" << endl;
            cout << "\t -shorthisto FILE.root      fill eventio file contents into histograms (compact version)" << endl;
            cout << "\t -smallfile            try to reduce histogram outputfile (reduction > factor of 3)" << endl;
            cout << "\t -muon                 set histogram limits for muons" << endl;
            cout << "\t -absfile              use atmospheric absorption routines from this extinction file (full path and file; default: ./data/us76.50km.ext)" << endl;
            cout << "\t                       (use '-absfile noExtinction' to ignore atmospheric extinction)" << endl;
            cout << "\t -queff FLOAT[0,1]     apply global quantum efficiency" << endl;
            cout << "\t -nevents INT          read only nevents events" << endl;
            cout << "\t -narray INT           read only narray arrays per event" << endl;
            cout << "\t -tel INT              telescope number to be processed (<0: process all telescopes, -1: output into one file; -2: one file per telescope" << endl;
            cout << "\t -seed INT             set seed for random generators" << endl;
            cout << "\t -COCO                 fill photon impact coordinates in CORSIKA coordinates" << endl;
            cout << "\t -verbose              print parameters for each event (default off)" << endl;
            cout << "\t -cfg FILENAME         grisu configuration file (only needed when telescope numbering in corsika and grisudet is different)" << endl;
            cout << "\t -printmoreinfo INT    print additional information (corsika event number & depth) into grisudet photon file, needs atmosphere number" << endl;
            /*         cout << endl << "(unsigned int " << numeric_limits<unsigned int>::max() << ")" << endl;
                     cout << endl << "(unsigned short int " << numeric_limits<unsigned short int>::max() << ")" << endl;
                     cout << endl << "(uint8_t " << numeric_limits<uint8_t>::max() << ")" << endl;
                     cout << endl << "(short " << numeric_limits<short>::max() << ")" << endl;
                     cout << endl << "(int " << numeric_limits<int >::max() << ")" << endl;
                     cout << endl << "(long " << numeric_limits<long>::max() << ")" << endl; */
            cout << endl;
            exit( 0 );
        }
        // grisu output file
        else if( ( iTemp.find( "-grisu" ) < iTemp.size() || iTemp.find( "-ioread" ) < iTemp.size() ) && iTemp2.size() > 0 )
        {
            bGRISU = true;
            fGrisuOutputFile = iTemp2;
            if( iTemp2 == "stdout" )
            {
                bstdout = true;
            }
            i++;
        }
        else if( iTemp.find( "-verbo" ) < iTemp.size() )
        {
            bDebug = true;
        }
        else if( iTemp.find( "-histo" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            if( bHisto )
            {
                cout << "error: histo already set (use -histo or -xyz, but not both)" << endl;
                exit( 0 );
            }
            bHisto = true;
            fHisto->init( iTemp2, false );
            i++;
        }
        else if( iTemp.find( "-xyz" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            if( bHisto )
            {
                cout << "error: histo already set (use -histo or -xyz, but not both)" << endl;
                exit( 0 );
            }
            bHisto = true;
            fHisto->init( iTemp2, false );
            // calculate levels for XYZ histograms
            fHisto->setXYZlevelsHeight( getXYZlevels( 5., 40. ) );
            fHisto->setXYZlevelsThickness( 5., 40. );
            fHisto->initXYZhistograms();
            i++;
        }
        else if( iTemp.find( "-shorthisto" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            bHisto = true;
            fHisto->init( iTemp2, true );
            i++;
        }
        else if( iTemp.find( "-COCO" ) < iTemp.size() )
        {
            fHisto->setCORSIKAcoordinates();
        }
        else if( iTemp.find( "-smallfile" ) < iTemp.size() )
        {
            //(preli)
            if( bHisto )
            {
                cout << "set small file command before -histo" << endl;
                exit( 0 );
            }
            fHisto->set_small_file();
        }
        else if( iTemp.find( "-muon" ) < iTemp.size() )
        {
            fHisto->setMuonSettings();
        }
        else if( iTemp.find( "-abs" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            fAtmosFile = iTemp2;
            if( fAtmosFile == "kascade" )
            {
                fAtmosFile = "data/kextint.dat";
                fAtmosModel = "kascade";
            }
            else if( fAtmosFile == "us76_new" )
            {
                fAtmosFile = "data/us76.50km.ext";
                fAtmosModel = "modtran4";
            }
            else if( fAtmosFile == "corsika" || fAtmosFile == "CORSIKA" )
            {
                fAtmosFile = "data/atmabs.dat";
                fAtmosModel = "corsika";
            }
            else if( fAtmosFile.find( "M5" )  < fAtmosFile.size() )
            {
                fAtmosModel = "modtran5";
            }
            else if( fAtmosFile == "noExtinction" )
            {
                fAtmosModel = "noExtinction";
            }
            else
            {
                fAtmosModel = "modtran4";
            }
            i++;
        }
        // global quantum efficiency should be between 0. and 1.
        else if( iTemp.find( "-queff" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            queff = atof( iTemp2.c_str() );
            i++;
            if( queff < 0. || queff > 1. )
            {
                cout << "invalid global quantum efficiency (0.<qeff<1.): " << queff << endl;
                exit( -1 );
            }
        }
        else if( iTemp.find( "-nevents" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            nevents = atoi( iTemp2.c_str() );
            i++;
        }
        else if( iTemp.find( "-narray" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            narray = atoi( iTemp2.c_str() );
            i++;
        }
        else if( iTemp.find( "-tel" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            nTel = atoi( iTemp2.c_str() );
            i++;
        }
        else if( iTemp.find( "-seed" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            fSeed = atoi( iTemp2.c_str() );
            i++;
        }
        else if( iTemp.find( "-cors" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            fCorsikaIO = iTemp2;
            i++;
            if( fCorsikaIO.size() <= 0 )
            {
                cout << "no CORSIKA input file, exiting... " << endl;
                exit( 0 );
            }
        }
        else if( iTemp.find( "-cfg" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            fGrisuConfigurationFile = iTemp2;
            i++;
        }
        else if( iTemp.find( "-printmoreinfo" ) < iTemp.size() && iTemp2.size() > 0 )
        {
            atmid = atoi( iTemp2.c_str() );
            i++;
            bPrintMoreInfo = true;
        }
        else if( i > 1 )
        {
            cout << "unknown run parameter: " << iTemp << endl;
            exit( -1 );
        }
    }
    if( !bstdout )
    {
        cout << fVersion << endl;
        cout << "=======================" << endl << endl;
        cout << "Based on skeleton.c in the CORSIKA IACT package. Use of Konrad Bernloehr's IACT routines, which are part of the CORSIKA distribution.";
        cout << endl << endl << endl;
    }
    
    TRandom3 fRandom( fSeed );
    if( !bstdout )
    {
        cout << "SEED (for Cherenkov photon wavelengths): " << ( int )fRandom.GetSeed() << endl;
    }
    if( !bstdout )
    {
        cout << "ntel mode " << nTel << endl;
    }
    
    // set the atmospheric absorption model
    VAtmosAbsorption fAtabso( fAtmosModel, fSeed, fAtmosFile );
    
    /* I/O buffer for input needed */
    if( ( iobuf = allocate_io_buffer( 0 ) ) == NULL )
    {
        fprintf( stderr, "Input I/O buffer not allocated\n" );
        exit( 1 );
    }
    iobuf->max_length = numeric_limits<long>::max();
    //   iobuf->max_length = 10000000000;
    if( !bstdout )
    {
        printf( "Maximum buffer length: %ld\n", iobuf->max_length );
    }
    if( !bstdout )
    {
        cout << "Photon output file " << fGrisuOutputFile << endl;
    }
    
    // try to open Corsika file
    data_file = 0;
    if( !bstdout )
    {
        printf( "Input file: %s\n", fCorsikaIO.c_str() );
    }
    if( ( data_file = fopen( fCorsikaIO.c_str(), "r" ) ) == NULL )
    {
        perror( fCorsikaIO.c_str() );
        exit( 1 );
    }
    
    iobuf->input_file = data_file;
    
    for( ;; ) /* Loop over all data in the input file */
    {
        /* Find and read the next block of data. */
        /* In case of problems with the data, print an error message and skip the rest of the file. */
        
        //possible return values: 0 (O.k.),  -1 (error),  or  -2 (end-of-file)
        int i_find = find_io_block( iobuf, &block_header );
        
        if( i_find != 0 )
        {
            switch( i_find )
            {
                case -1:
                    cerr << "corsikaIOreader: There was an error finding the next IO block; will skip the rest of the input." << endl;
                    break;	//break switch statement.
                    
                case -2:
                    break;	//break switch statement.
                    
                default:
                    cerr << "corsikaIOreader: There was an undefined error finding the next IO block; find_io_block returned " << i_find << "." << endl;
            }
            break;		//break for loop.
        }
        
        //possible return values:  0 (O.k.), -1 (error), -2 (end-of-file), -3 (block skipped because it is too large).
        int i_block = read_io_block( iobuf, &block_header );
        
        if( i_block != 0 )
        {
            switch( i_block )
            {
                case -1:
                    cerr << "corsikaIOreader: There was an error reading the next IO block; will skip the rest of the input." << endl;
                    break;	//break switch statement.
                    
                case -2:
                    break;	//break switch statement.
                    
                case -3:
                    cerr << "corsikaIOreader: There was an error reading the next IO block (block skipper because it is too large); will skip the rest of the input." << endl;
                    break;	//break switch statement.
                    
                default:
                    cerr << "corsikaIOreader: There was an undefined error finding the next IO block; find_io_block returned " << i_find << "." << endl;
            }
            break;		//break for loop.
        }
        
        
        /* What did we actually get? */
        switch( block_header.type )
        {
            /* CORSIKA run header */
            case IO_TYPE_MC_RUNH:
                if( bPrintHeaders )
                {
                    cout << "IO_TYPE_MC_RUNH" << endl;
                }
                read_tel_block( iobuf, IO_TYPE_MC_RUNH, runh, 273 );
                {
                    int nht = ( int )runh[4];
                    if( nht > 0 && nht <= 10 )
                    {
                        array.obs_height = runh[4 + nht];
                    }
                    else
                    {
                        array.obs_height = -100;
                    }
                }
                run = ( int ) runh[1];
                if( bDebug )
                {
                    printf( "Run %d: observation level is at %6.1f m\n", run, 0.01 * array.obs_height );
                }
                airlightspeed = 29.9792458 / Nair( 1e-5 * array.obs_height );
                // elow = runh[16];
                if( bDebug )
                {
                    printf( "Events created between %5.3f and %5.3f TeV\n", ( double )runh[16] / 1e3, ( double )runh[17] / 1e3 );
                }
                /* CORSIKA run information in run header */
                array.mc_run.height = array.obs_height * 0.01;
                array.mc_run.e_min = runh[16] * 0.001;
                array.mc_run.e_max = runh[17] * 0.001;
                array.mc_run.slope = runh[15];
                /* Further information has to wait for event header */
                array.mc_run.radius = 0.;
                array.mc_run.num_arrays = 0;
                array.mc_run.theta_min = array.mc_run.theta_max = -1.;
                array.mc_run.phi_min = array.mc_run.phi_max = -1.;
                array.mc_run.wlen_min = array.mc_run.wlen_max = 0.;
                
                // observation level from Corsika
                if( array.obs_height > 0. )
                {
                    fAtabso.setObservationlevel( array.obs_height * 0.01 );
                    for( unsigned int p = 0; p < fGrisu.size(); p++ )
                    {
                        fGrisu[p]->setObservationHeight( array.obs_height * 0.01 );
                    }
                }
                
                // fill corsika run header class
                fRunHeader->runnumber = ( unsigned int )runh[1];
                fRunHeader->production_date = ( unsigned int )runh[2];
                fRunHeader->corsika_version = ( float )runh[3];
                fRunHeader->observation_level_m = ( float )array.mc_run.height;
                fRunHeader->energy_slope = ( float )array.mc_run.slope;
                fRunHeader->energy_min_GeV = ( float )array.mc_run.e_min;
                fRunHeader->energy_max_GeV = ( float )array.mc_run.e_max;
                fRunHeader->xscatt_m = ( float )runh[247] * 0.01;
                fRunHeader->yscatt_m = ( float )runh[248] * 0.01;
                
                break;
                
            /* CORSIKA inputs */
            case IO_TYPE_MC_INPUTCFG:
                if( bPrintHeaders )
                {
                    cout << "IO_TYPE_MC_INPUTCFG" << endl;
                }
                read_input_lines( iobuf, &corsika_inputs );
                if( corsika_inputs.text != NULL )
                {
                    struct linked_string* xl, *xln;
                    for( xl = &corsika_inputs; xl != NULL; xl = xln )
                    {
                        free( xl->text );
                        xl->text = NULL;
                        xln = xl->next;
                        xl->next = NULL;
                        if( xl != &corsika_inputs )
                        {
                            free( xl );
                        }
                    }
                    fflush( stdout );
                }
                break;
                
            /* Telescope positions (relative positions in array) */
            case IO_TYPE_MC_TELPOS:
                if( bPrintHeaders )
                {
                    cout << "IO_TYPE_MC_TELPOS" << endl;
                }
                read_tel_pos( iobuf, MAX_TEL, &array.ntel, array.xtel, array.ytel, array.ztel, array.rtel );
                fTelescopeMatrix = readTelescopeMatrix( fGrisuConfigurationFile, array.ntel, array.xtel, array.ytel );
                //	    cout << "telescope matrix: " << endl;
                //	    for( unsigned int i = 0; i < fTelescopeMatrix.size(); i++ ) if( fTelescopeMatrix[i] >= 0 ) cout << "Tel: " << i+1 << "\t" << fTelescopeMatrix[i]+1 << endl;
                
                // initialize grisu readers:
                if( fGrisu.size() == 0 )
                {
                    if( nTel > -2 )
                    {
                        fGrisu.push_back( new VGrisu( fVersion, atmid ) );
                        if( fGrisuOutputFile.size() > 0 )
                        {
                            fGrisu.back()->setOutputfile( fGrisuOutputFile );
                        }
                    }
                    else if( nTel == -2 )
                    {
                        char hO[2000];
                        for( int pt = 0; pt < array.ntel; pt++ )
                        {
                            fGrisu.push_back( new VGrisu( fVersion, atmid ) );
                            if( fGrisuOutputFile.size() > 0 )
                            {
                                sprintf( hO, "%s_%d", fGrisuOutputFile.c_str(), pt + 1 );
                                fGrisu.back()->setOutputfile( hO );
                            }
                        }
                    }
                    for( unsigned int pt = 0; pt < fGrisu.size(); pt++ )
                    {
                        fGrisu[pt]->setQueff( queff );
                    }
                }
                break;
                
            /* CORSIKA event header */
            case IO_TYPE_MC_EVTH:
                if( bPrintHeaders )
                {
                    cout << "IO_TYPE_MC_EVTH" << endl;
                }
                if( bDebug )
                {
                    cout << "reading IO_TYPE_MC_EVTH" << endl;
                }
                read_tel_block( iobuf, IO_TYPE_MC_EVTH, evth, 273 );
                // event = ( int )evth[1];
                wl_lower_limit = evth[95];
                wl_upper_limit = evth[96];
                primary_energy = evth[3];
                EVTH76 = ( unsigned long int )evth[76];
                if( EVTH76.test( 2 ) && bCEFFICWARNING )
                {
                    cout << endl;
                    cout << "WARNING: ignoring any efficiencies applied in CORSIKA (CEFFIC options)" << endl;
                    cout << endl;
                    bCEFFICWARNING = false;
                }
                fRunHeader->particleID = ( unsigned int )evth[2];
                fRunHeader->startingaltitude_gcm2 = ( float )evth[4];
                fRunHeader->tstart = ( int )evth[6];
                fRunHeader->nscatt = ( unsigned int )evth[97];
                fRunHeader->zenith_min_deg = ( float )evth[80];
                fRunHeader->zenith_max_deg = ( float )evth[81];
                fRunHeader->azimuth_min_deg = ( float )evth[82];
                fRunHeader->azimuth_max_deg = ( float )evth[83];
                fRunHeader->viewcone_min_deg = ( float )evth[152];
                fRunHeader->viewcone_max_deg = ( float )evth[153];
                fRunHeader->cherenkov_flag = ( unsigned long int )( evth[76] + 0.5 );
                fRunHeader->cherenkov_bunchsize = ( float )evth[84];
                fRunHeader->cherenkov_bandwidth_min_nm = ( float )evth[95];
                fRunHeader->cherenkov_bandwidth_max_nm = ( float )evth[96];
                fRunHeader->geomagneticfield_arrang_deg = ( float )evth[92];
                fRunHeader->geomagneticfield_x_muT = ( float )evth[70];
                fRunHeader->geomagneticfield_z_muT = ( float )evth[71];
                fRunHeader->hadronic_model_low = ( unsigned int )evth[74];
                fRunHeader->hadronic_model_high = ( unsigned int )evth[75];
                fRunHeader->hadronic_model_transition_energy_GeV = ( float )evth[154];
                if( !have_atm_profile )
                {
                    // read atmospheric profile from CHERENKOV OPTIONS
                    EVTH76 >>= 10;
                    int iatmo = ( int )EVTH76.to_ulong();
                    if( iatmo <= 0 )
                    {
                        iatmo = 6;    /* US standard atmosphere */
                    }
                    try
                    {
                        atmset_( &iatmo, &array.obs_height );
                    }
                    catch(...)
                    {
                        cout << "error initialising atmospheres" << endl;
                        cout << "...exiting" << endl;
                        exit( EXIT_FAILURE );
                    }
                    have_atm_profile = 1;
                }
                
                array.shower_sim.energy = 0.001 * primary_energy; /* in TeV */
                array.shower_sim.xmax = array.shower_sim.emax =
                                            array.shower_sim.cmax = 0.;
                array.shower_sim.hmax = 0.;
                particle_type = Nint( evth[2] );
                // toffset = ( evth[6] - array.obs_height ) / cos( evth[10] ) / 29.9792458;
                if( bDebug )
                {
                    printf( "Event %d: particle type %d with energy %5.2f TeV \n", ( int )evth[1], particle_type, 0.001 * primary_energy );
                }
                alt = 90. - ( 180. / M_PI ) * evth[10];
                //(GM) keep Corsika coordinate system            az  = 180. - (180./M_PI)*(evth[11]-evth[92]);
                az = ( 180. / M_PI ) * ( evth[11] - evth[92] );
                az -= floor( az / 360. ) * 360.;
                if( bDebug )
                {
                    printf( "   zenith angle %4.2f deg, azimuth %4.2f deg\n", 90. - alt, az );
                }
                
                array.mc_run.theta_min = evth[80];
                array.mc_run.theta_max = evth[81];
                array.mc_run.phi_max = 180. - ( evth[82] - evth[92] );
                array.mc_run.phi_min = 180. - ( evth[83] - evth[92] );
                array.mc_run.bunchsize = evth[84];
                array.mc_run.wlen_min = evth[95];
                array.mc_run.wlen_max = evth[96];
                
                array.shower_sim.azimuth = az;
                array.shower_sim.altitude = alt;
                array.shower_sim.firstint = abs( evth[6] * 0.01 );
                array.shower_sim.shower_id = evth[1];
                array.shower_sim.particle = particle_type;
                
                if( readNevent == 0 && bGRISU )
                {
                    for( unsigned int p = 0; p < fGrisu.size(); p++ )
                    {
                        fGrisu[p]->setObservationHeight( array.obs_height * 0.01 );
                        fGrisu[p]->writeRunHeader( evth, fRunHeader );
                    }
                }
                readNarray = 0;
                break;
                
            /* Offsets of telescope array instances for the following event */
            case IO_TYPE_MC_TELOFF:
                if( bPrintHeaders )
                {
                    cout << "IO_TYPE_MC_TELOFF" << endl;
                }
                read_tel_offset( iobuf, MAX_ARRAY, &array.narray, &array.toff, array.xoff, array.yoff );
                // toffset = array.toff; /* Should be about the same again as above */
                array.mc_run.num_arrays = array.narray;
                if( bDebug )
                {
                    cout << "\t total number of arrays simulated: " << array.mc_run.num_arrays << endl;
                }
                break;
                
            /* Photon data for a complete array (one of perhaps many instances) */
            case IO_TYPE_MC_TELARRAY:
                if( bPrintHeaders )
                {
                    cout << "IO_TYPE_MC_TELARRAY" << endl;
                }
                if( readNarray >= narray && narray >= 0 )
                {
                    break;
                }
                
                begin_read_tel_array( iobuf, &item_header, &iarray );
                
                array.shower_sim.xcore = -0.01 * array.xoff[iarray]; /* in meters */
                array.shower_sim.ycore = -0.01 * array.yoff[iarray]; /* in meters */
                if( bDebug )
                {
                    cout << "\t shower core for array " << readNarray << " at (x/m) [m]: " << array.shower_sim.xcore << "\t" << array.shower_sim.ycore << endl;
                }
                /* Observation level is now defining z=0.: */
                array.shower_sim.zcore = 0.; /* Note: this is below lowest telescope */
                
                array.shower_sim.core_dist_3d =
                    line_point_distance( array.shower_sim.xcore,
                                         array.shower_sim.ycore, array.shower_sim.zcore,
                                         cx = -1.*cos( alt * ( M_PI / 180. ) ) * cos( az * ( M_PI / 180. ) ),
                                         cy = -1.*cos( alt * ( M_PI / 180. ) ) * sin( az * ( M_PI / 180. ) ),
                                         cz = sin( alt * ( M_PI / 180. ) ), 0.01 * array.refpos[0],
                                         0.01 * array.refpos[1], 0.01 * array.refpos[2] );
                /* Distances of telescopes from shower axis */
                for( int itel = 0; itel < array.max_tel; itel++ )
                {
                    array.shower_sim.tel_core_dist_3d[itel] = line_point_distance( array.shower_sim.xcore, array.shower_sim.ycore, array.shower_sim.zcore, cx, cy, cz, 0.01 * array.xtel[itel], 0.01 * array.ytel[itel], 0.01 * array.ztel[itel] );
                }
                
                if( bHisto )
                {
                    fHisto->newEvent( evth, array, iarray );    // start new event for each array
                }
                if( bGRISU )
                {
                    for( unsigned int p = 0; p < fGrisu.size(); p++ )
                    {
                        fGrisu[p]->writeEvent( array, bPrintMoreInfo );
                    }
                }
                
                for( itc = 0; itc < array.ntel; itc++ )
                {
                    sub_item_header.type = IO_TYPE_MC_PHOTONS;
                    if( search_sub_item( iobuf, &item_header, &sub_item_header ) < 0 )
                    {
                        break;
                    }
                    /* Read the photon bunches for one telescope */
                    int itel = 0;
                    if( read_tel_photons( iobuf, MAX_BUNCHES, &jarray, &itel, &photons, bunches, &nbunches ) < 0 )
                    {
                        fprintf( stderr, "Error reading %d photon bunches\n", nbunches );
                        continue;
                    }
                    if( nTel >= 0 && itel != nTel )
                    {
                        continue;
                    }
                    // fill number of photons per telescope (ignore telescope matrix, this is for tcors!)
                    //if( bHisto )
                    //{
                    //    fHisto->fillNPhotons( itel, photons );
                    //}
                    
                    // check if this telescope should be analysed
                    if( itel < ( int )fTelescopeMatrix.size() && fTelescopeMatrix[itel] < 0 )
                    {
                        continue;
                    }
                    
                    // loop over all bunches
                    for( ibunch = 0; ibunch < nbunches; ibunch++ ) // loop over all bunches for this telescope
                    {
                        wl_bunch = bunches[ibunch].lambda;
                        cx = bunches[ibunch].cx;
                        cy = bunches[ibunch].cy;
                        cz = -1.*sqrt( 1. - cx * cx - cy * cy ); /* direction is downwards */
                        /* Use secans(zenith angle) for airmass, */
                        /* i.e. assume a plane atmosphere. */
                        if( cz != 0. )
                        {
                            airmass = -1. / cz;
                        }
                        else
                        {
                            airmass = 1.e16;
                        }
                        /* Distance between point of emission and */
                        /* the CORSIKA observation level */
                        // distance = ( bunches[ibunch].zem - array.obs_height ) * airmass;  // bunches[ibunch].zem is above sea level
                        /* Distance between CORSIKA observation level and */
                        /* telescope fixed position. */
                        tel_dist = array.ztel[itel] * airmass;
                        /* Note that, although tracing starts at the CORSIKA */
                        /* level, the bunch time corresponds to the crossing */
                        /* of the telescope level. */
                        tel_delay = tel_dist / airlightspeed;
                        /* Note also that the photon bunch might be created */
                        /* behind the telescope mirror. Check in raytracing. */
                        
                        // (GM) restore arrival time at ground:
                        // add travel time from telescope plane to ground plane
                        corstime = bunches[ibunch].ctime + tel_delay;
                        
                        // fill all bunch specific stuff into histograms
                        if( bHisto )
                        {
                            fHisto->fillBunch( bunches[ibunch], corstime );
                        }
                        // now loop over bunch
                        for( ; bunches[ibunch].photons > 0; bunches[ibunch].photons -= 1. )
                        {
                            // photon wavelength
                            if( wl_bunch == 0. )
                            {
                                /* get photon wavelength according to 1./lambda^2 distribution */
                                lambda = 1. / ( 1. / wl_lower_limit - fRandom.Uniform( 1. ) * ( 1. / wl_lower_limit - 1. / wl_upper_limit ) );
                            }
                            else if( wl_bunch < 0. )
                                /* This indicates that quantum efficiency, mirror */
                                /* reflectivity, and atmospheric transmission */
                                /* have already been applied in CORSIKA (which */
                                /* was CMZ extracted then with the CEFFIC option). */
                            {
                                // IGNORE CEFFIC!!!!
                                //		        if( cherenkov_flag == 6175 )
                                {
                                    lambda = 1. / ( 1. / wl_lower_limit - fRandom.Uniform( 1. ) * ( 1. / wl_lower_limit - 1. / wl_upper_limit ) );
                                }
                                /*			else
                                			{
                                			   lambda = -1.;
                                                        } */
                            }
                            else
                                /* Wavelength already generated in Corsika */ // (GM) for non-standard CORSIKA
                            {
                                lambda = wl_bunch;
                            }
                            
                            // atmospheric extinction
                            if( lambda >= 1000 )
                            {
                                continue;
                            }
                            else if( lambda >= 0 )
                            {
                                prob = fAtabso.probAtmAbsorbed( lambda, ( double )bunches[ibunch].zem * 0.01, -1. * cz );
                            }
                            else
                            {
                                prob = 1.;
                            }
                            // fill photon structure
                            Chphoton.photons = 1.;
                            Chphoton.x = bunches[ibunch].x * 0.01 + array.xtel[itel] * 0.01;
                            Chphoton.y = bunches[ibunch].y * 0.01 + array.ytel[itel] * 0.01;
                            Chphoton.cx = bunches[ibunch].cx;
                            Chphoton.cy = bunches[ibunch].cy;
                            Chphoton.ctime = corstime;
                            Chphoton.zem = bunches[ibunch].zem * 0.01;
                            Chphoton.lambda = lambda;

                            // fill generated photons into histograms
                            if( bHisto )
                            {
                                fHisto->fillGenerated( Chphoton, prob );
                            }
                            // extinction + efficiencies
                            if( bunches[ibunch].photons < 1. )
                            {
                                prob *= bunches[ibunch].photons;
                            }
                            if( prob <= 1. )
                            {

                                double iRand = fRandom.Uniform( 1. );
                                if( iRand > prob )
                                {
                                    continue;
                                }
                                //// fill histograms without quantum efficiency applied (only atm.extinction)
                                // moved to after QE and lens transmission by NK
                                //if( bHisto )
                                //{
                                //    fHisto->fillSurvived( Chphoton, prob );
                                //}
                                //apply global quantum efficiency
                                prob *= queff;
                                if( iRand > prob )
                                {
                                    continue;
                                }
                                /* NK
                                * PANOSETI quantum efficiency
                                */
                                // apply PANOSETI quantum efficiency
                                prob *= 0.9189 / (1. + ( exp(-0.2046*(lambda-384.2)) ) );
                                if( iRand > prob )
                                {
                                    continue;
                                }
                                // apply PANOSETI lens transmission
                                prob *= ( (-3.244e-11 * pow(lambda,4) ) + (9.376e-8 * pow(lambda,3) ) + (-9.880e-5 * pow(lambda,2) ) + (4.402e-2 * lambda) - 6.623 );
                                if( iRand > prob )
                                {
                                    continue;
                                }
                                // fill number of photons per telescope (after extinction)
                                if( bHisto )
                                {
                                    fHisto->fillNPhotons( itel, 1.0 );
                                    fHisto->fillSurvived( Chphoton, prob, evth, itel );
                                }
                                // write photons to iotxt output file (after quantum efficiency)
                                if( bGRISU )
                                {
                                    if( nTel > -2 )
                                    {
                                        if( fGrisu.size() == 1 )
                                        {
                                            fGrisu[0]->writePhotons( Chphoton, fTelescopeMatrix[itel] );
                                        }
                                    }
                                    else if( nTel == -2 )
                                    {
                                        // move all photons around coordinates centre
                                        Chphoton.x -= array.xtel[itel] / 1.e2;
                                        Chphoton.y -= array.ytel[itel] / 1.e2;
                                        // telescope ID is always 0
                                        if( itel < ( int )fGrisu.size() )
                                        {
                                            fGrisu[itel]->writePhotons( Chphoton, 0 );
                                        }
                                    }
                                }
                            }
                        }
                    }
                } /* End of loop over telescopes */
                end_read_tel_array( iobuf, &item_header );
                
                readNarray++;
                break;
                
            /* CORSIKA event trailer */
            case IO_TYPE_MC_EVTE:
                if( bPrintHeaders )
                {
                    cout << "IO_TYPE_MC_EVTE" << endl;
                }
                read_tel_block( iobuf, IO_TYPE_MC_EVTE, evte, 273 );
                readNevent++;
                /* All array instances for this shower are finished */
                break;
                
            /* CORSIKA run trailer */
            case IO_TYPE_MC_RUNE:
                if( bPrintHeaders )
                {
                    cout << "IO_TYPE_MC_RUNE" << endl;
                }
                read_tel_block( iobuf, IO_TYPE_MC_RUNE, rune, 273 );
                break;
                
            /* Unknown / any other material */
            default:
                if( bPrintHeaders )
                {
                    cout << "Read block header of type " << block_header.type << ", will ignore." << endl;
                }
                break;
        } /* End of switch over all input data types */
        if( readNevent >= nevents && nevents > 0 )
        {
            break;
        }
    } /* End of loop over all data in the input file */
    fclose( iobuf->input_file );
    iobuf->input_file = NULL;
    if( bHisto && fHisto )
    {
        fHisto->terminate();
    }
    if( !bstdout )
    {
        if( fRunHeader )
        {
            fRunHeader->printHeader( cout );
        }
        cout << endl;
        cout << "END OF RUN ( " << readNevent << " showers )" << endl;
    }
    
    return 0;
}
