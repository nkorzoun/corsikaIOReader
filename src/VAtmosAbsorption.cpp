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

    \class VAtmosAbsorption
    \brief atmospheric absorption of Cherenkov photons

     equiv. to CORSIKA function ATABSO or kascade function atm_pass
          (routines are copied from these programes)

     sourcefiles for extinction values are:
        - atmabs.dat   default CORSIKA data
	- kextint.dat  default cherenkf dat
     ( S.Valley (ed.), Handbook of Geophysics and Space Enviornments, Cambridge MA, McGraw-Hill. 1965
            Table 7.4 (270-900nm) )
	- us76.50km.ext MODTRAN4 calculation of Michael Daniel (from 200nm)

     assume extinction values are constant above 50 km

     all units for input parameters are in meter and nanometer

     ROOT is only needed for random generator (and histogramming classes)

    \attention
      finetuned to tables in extinction values files - do not change

    \section example example

    \code
    #include "VAtmosAbsorption.h"
    ....

    // use extinction values from MODTRAN4
     VAtmosAbsorption *atabso = new VAtmosAbsorption( "us76_new", 0 );
     // set observation height to 1500m
     atabso->setObservationlevel( 1500. );
     // set wavelength intervall in [nm]
     atabso->setWavelengthintervall( 185., 700. );

     // get photon wavelength from emission height 5km and with direction cosinus 0.95 (at ground)
     cout <<  fatmabs->getWavelength( 5000., 0.95 ) << endl;

     \endcode

    \observe
     very sensitive to input file format!

    \date
      08/05/2004

    \author
      Gernot Maier

*/


#include "VAtmosAbsorption.h"

/*!
    \param model atmospheric absorption from "CORSIKA", "kascade", us76_new
*/
VAtmosAbsorption::VAtmosAbsorption( string model, int fSeed, string iSourceFile )
{
    fModel = model;
    fObservationLevel = 1000.;       // default observation level in [m]
    fRandom = new TRandom3( fSeed );
    fminWave = 300.;                 // default CORSIKA values
    fmaxWave = 450.;                 // default CORSIKA values
    
    cerr << "Atmospheric extinction model : " << model << endl;
    
    // CORSIKA extinction values
    
    if( fModel == "modtran5" )
    {
        fSourceFile = iSourceFile;
        read_extint_M5( );
    }
    
    else if( fModel == "CORSIKA" || fModel == "corsika" )
    {
        fModel = "corsika";
        
        if( iSourceFile.size() > 0 )
        {
            fSourceFile = iSourceFile;
        }
        else
        {
            fSourceFile = "data/atmabs.dat";
        }
        
        readCorsikaAtmabs();
        
    }
    // kascade extinction values
    else if( fModel == "kascade" )
    {
        if( iSourceFile.size() > 0 )
        {
            fSourceFile = iSourceFile;
        }
        else
        {
            fSourceFile = "data/kextint.dat";
        }
        // for the original kextint, this should be 180.
        read_extint( 180 );
        
    }
    // MODTRAN4 US76 extinction values
    else if( fModel == "us76_new" || fModel == "us76.50km" || fModel == "us76.50" || fModel == "modtran4" )
    {
        fModel = "modtran4";
        
        if( iSourceFile.size() > 0 )
        {
            fSourceFile = iSourceFile;
        }
        else
        {
            fSourceFile = "data/us76.50km.ext";
        }
        read_extint( 200 );
    }
    else if( fModel == "us76.23km" || fModel == "us76.23" || fModel == "modtran4_2" )
    {
        fModel = "modtran4_2";
        
        if( iSourceFile.size() > 0 )
        {
            fSourceFile = iSourceFile;
        }
        else
        {
            fSourceFile = "data/us76.23km.ext";
        }
        read_extint_F2( 200 );
    }
    else if( fModel == "artemis" )
    {
        if( iSourceFile.size() > 0 )
        {
            fSourceFile = iSourceFile;
        }
        else
        {
            fSourceFile = "data/extinction_uv.dat";
        }
        read_extint( 180 );
    }
    else if( fModel == "noExtinction" )
    {
        fSourceFile = "noExtinction";
        cerr << "VAtmosAbsorption: no atmospheric extinction applied" << endl;
    }
    else
    {
        cout << "VAtmosAbsorption::VAtmosAbsorption: error, unknown model: " << model << endl;
        exit( -1 );
    }
}

/*!
   \param obslevel observation level in [m]
*/
void VAtmosAbsorption::setObservationlevel( double obslevel )
{
    fObservationLevel = obslevel;
    
    if( fModel == "corsika" )
    {
        int xobs = ( int )( obslevel / 1000 );
        map< int, vector<double> >::const_iterator m_iter;
        for( m_iter = fCoeff.begin(); m_iter != fCoeff.end(); ++m_iter )
        {
            fCoeffObs[m_iter->first] = getLinearInterpolate( obslevel / 1000., xobs, xobs + 1, m_iter->second[xobs], m_iter->second[xobs + 1] );
        }
    }
}

/*!
   calulates survival probability for photons

   \param  wl wavelength in nm
   \param  zemis emission height in m
   \param  wemis cos of emission angle (cos theta)
   \return survival probability for photon
*/

double VAtmosAbsorption::probAtmAbsorbed( double wl, double zemis, double wemis )
{
    double a = 0;
    return probAtmAbsorbed( wl, zemis, wemis, a );
}

double VAtmosAbsorption::probAtmAbsorbed( double wl, double zemis, double wemis, double& optdepth )
{
    // fixed altitude steps [m]
    double fAltitudeStep = 1000.;
    // optical depth
    optdepth = 0.;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // copy from CORSIKA (translated to C++)
    if( fModel == "corsika" )
    {
        const double fWlMin = 180.;
        const double fWlStep = 5.;
        int riwl;
        int wli0;
        int wli1;
        int hti0;
        int hti1;
        double htkm;
        double coatex;
        double fx0;
        double fx1;
        double phi0 = 0.;
        double phi1 = 0.;
        double probs;
        
        //  CALCULATE THE REFERENCE WL AND INDEX OF WL FOR THE INTERPOLATIONS
        riwl = 1 + ( int )( ( wl - fWlMin ) / fWlStep );
        wli0 = riwl * ( int )fWlStep + ( int )( fWlMin - fWlStep );
        wli1 = riwl * ( int )fWlStep + ( int )( fWlMin );
        
        if( fCoeff.find( wli0 ) !=  fCoeff.end() && fCoeff.find( wli1 ) !=  fCoeff.end() )
        {
            // CONSIDER ATMOSPHERIC EXTINCITION
            htkm = zemis / 1000;   // m -> km
            hti0 = ( int )htkm;
            hti1 = ( int )htkm + 1;
            
            if( hti0 < 0 )
            {
                phi0 = fCoeff[wli0][0];
                phi1 = fCoeff[wli1][0];
            }
            else if( hti1 > 50 )
            {
                if( 50 < fCoeff[wli0].size() && 50 < fCoeff[wli1].size() )
                {
                    phi0 = fCoeff[wli0][50];
                    phi1 = fCoeff[wli1][50];
                }
            }
            else
            {
                // INTERPOLATION IN HEIGHT
                fx0 = fCoeff[wli0][hti0];
                fx1 = fCoeff[wli0][hti1];
                phi0 = getLinearInterpolate( htkm, ( double )hti0, ( double )hti1, fx0, fx1 );
                phi0 -= fCoeffObs[wli0];
                fx0 = fCoeff[wli1][hti0];
                fx1 = fCoeff[wli1][hti1];
                phi1 = getLinearInterpolate( htkm, ( double )hti0, ( double )hti1, fx0, fx1 );
                phi1 -= fCoeffObs[wli1];
            }
            coatex = getLinearInterpolate( wl, ( double )wli0, ( double )wli1, phi0, phi1 );
            optdepth = coatex / wemis;
            probs = exp( -1. * optdepth );
            return probs;
        }
        else
        {
            cout << "VAtmosAbsorption::isAbsorped: coeff. matrix not valid" << endl;
            exit( -1 );
        }
    }
    else if( fModel == "noExtinction" )
    {
        return 1.;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    else
    {
        // copy from GrIsu-code cherenk7.c
        /* int atm_pass(double z, double dn, double wave )
          Returns 1 if the photon survives the atmosphere, zero if it does not.
          The probability of passing is determined by the optical depth, which is
          calculated from the altitude of emission, the observatory altitude and
          the photon wavelength.                                               */
        
        int hobs = ( int )fObservationLevel;
        int z = ( int )zemis;
        int wave = ( int )wl;
        
        int    iwave1, iwave2, ihobs, ihgt;
        double tlow, thigh, atmprob;
        double p1, p2;
        unsigned int iext_index;
        
        /* Convert from m to nm. */
        //(GM) already in nm wave /= NM;
        
        int minwave = 180;
        // MODTRAN4 data is from 200nm only
        if( fModel == "modtran4" || fModel == "modtran4_2" )
        {
            minwave = 200;
        }
        if( fModel == "modtran5" )
        {
            minwave = 205;
        }
        
        // if wavelength is below minimal wavelength return 0
        if( wave < minwave )
        {
            return 0.;
        }
        
        // step size fixed to 5 nm (if this is not in the data file -> interpolation between values)
        iwave1 = ( wave - minwave ) / 5;
        iwave2 = ( wave - minwave ) / 5 + 1;
        
        ihobs = hobs / ( int )fAltitudeStep;
        
        
        p1 = extint[ihobs][iwave1] + ( extint[ihobs + 1][iwave1] - extint[ihobs][iwave1] ) * ( hobs / fAltitudeStep - ihobs );
        p2 = extint[ihobs][iwave2] + ( extint[ihobs + 1][iwave2] - extint[ihobs][iwave2] ) * ( hobs / fAltitudeStep - ihobs );
        tlow = getLinearInterpolate( wave, ( double )( iwave1 * 5 + minwave ), ( double )( ( iwave2 ) * 5 + minwave ), p1, p2 );
        
        ihgt = z / ( int )fAltitudeStep;
        iext_index = ( int )ihgt;
        
        // extinction calculated up to 50 km only (for modtran4 input format)
        if( iext_index > extint.size() - 2 )
        
        {
            // extinction coefficinents don't change much above 50km, use values of 50km
            iext_index = extint.size() - 2 ;
            //	   cout << "warning: photon production height above 50km, skipping" << endl;
            //	   return 0.;
        }
        p1 = extint[iext_index][iwave1] + ( extint[iext_index + 1][iwave1] - extint[iext_index][iwave1] ) * ( z / fAltitudeStep - ihgt );
        p2 = extint[iext_index][iwave2] + ( extint[iext_index + 1][iwave2] - extint[iext_index][iwave2] ) * ( z / fAltitudeStep - ihgt );
        thigh = getLinearInterpolate( wave, ( double )( iwave1 * 5 + minwave ), ( double )( ( iwave2 ) * 5 + minwave ), p1, p2 );
        
        /* tlow and thigh are, respectively, the optical depths at the observation
           and emission altitudes.                                              */
        if( wemis != 0.0 )
        {
            optdepth = -1.*( tlow - thigh ) / wemis;
            atmprob = TMath::Exp( -1. * optdepth );
        }
        else
        {
            atmprob = 0.;
        }
        // check validity of results
        if( !isnormal( atmprob ) )
        {
            cerr << "VAtmosAbsorption::probAtmAbsorbed not normal " << wl << "\t" << atmprob << endl;
            cerr << "\t tlow " << tlow << "\t thigh " << thigh << "\t wemis " << wemis << "\t optdepth " << optdepth << endl;
            cerr << "\t p1 " << p1 << "\t " << p2 << "\t iext_index " << iext_index << "\t iwave2 " << iwave2 << "\t z " << z << "\t fAltitudeStep " << fAltitudeStep << "\t ihgt " << ihgt << endl;
            cerr << "\t ve " << extint.size() << "\t" << extint[iext_index].size() << "\t" << extint[iext_index + 1].size() << endl;
            cerr << "\t extint[iext_index][iwave2] " << extint[iext_index][iwave2] << "\t extint[iext_index+1][iwave2] " << extint[iext_index + 1][iwave2] << endl;
            cerr << "\t inter " << wave << "\t" << ( double )( iwave1 * 5 + minwave ) << "\t" << ( double )( ( iwave2 ) * 5 + minwave ) << "\t" << p1 << "\t" << p2 << endl;
            return 0.;
        }
        return atmprob;
    }
    return 0.;
}

/*!
   \attention
      fine tuned to the CORSIKA V6.031 atmospheric extinction file atmabs.dat
*/
void VAtmosAbsorption::readCorsikaAtmabs()
{
    ifstream is( fSourceFile.c_str() );
    if( !is )
    {
        cout << "Atmospheric extinction file not found: " << fSourceFile << endl;
        exit( -1 );
    }
    
    cerr << "VAtmosAbsorption: reading atmospheric extinction file (from CORSIKA): " << fSourceFile << endl;
    
    string is_line;
    string is_Temp;
    
    int wl = 0;
    vector< double > i_coeff;
    int nStep = 0;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            if( is_line.size() < 5 )   // wavelength line
            {
                if( nStep > 0 )
                {
                    fCoeff[wl] = i_coeff;
                }
                istringstream is_stream( is_line );
                is_stream >> is_Temp;
                wl = atoi( is_Temp.c_str() );
                i_coeff.clear();
                nStep++;
            }
            else if( is_line.find( "EXTINCTION" ) >= is_line.size() )
            {
                istringstream is_stream( is_line );
                while( !is_stream.eof() )
                {
                    is_stream >> is_Temp;
                    i_coeff.push_back( atof( is_Temp.c_str() ) );
                }
            }
        }
    }
    if( nStep > 0 )
    {
        fCoeff[wl] = i_coeff;
    }
}

/*!
   from CORSIKA
*/
double VAtmosAbsorption::getLinearInterpolate( double x, double x0, double x1, double y0, double y1 )
{
    if( x0 == x1 )
    {
        return 0.;
    }
    return ( y0 + ( y1 - y0 ) * ( x - x0 ) / ( x1 - x0 ) );
}

/*!

*/
void VAtmosAbsorption::setWavelengthintervall( double imin, double imax )
{
    fminWave = imin;
    fmaxWave = imax;
}

/*!
   get a random wavelength from Jelly formula and apply atmospheric absorption

   \param emissionheigth    emission height of photon [m]
   \param emissionangle     emission angle of photon
   \return photon wavelength if photon survives atmosphere, otherwise -1
*/
double VAtmosAbsorption::getWavelength( double emissionheigth, double emissionangle )
{
    double lambda = 0.;
    double prob = 0.;
    
    lambda = 1. / ( 1. / fminWave - fRandom->Uniform( 1. ) * ( 1. / fminWave - 1. / fmaxWave ) );
    prob = probAtmAbsorbed( lambda, emissionheigth, emissionangle );
    if( fRandom->Uniform( 1. ) > prob )
    {
        return -1;
    }
    
    return lambda;
}

/*!
     copy from GrIsu-code cherenk7.c

     51 alitude steps of 1km, from
     wavelength step size is 5nm

     \param LAMBDA_MIN minimum of Cherenkov photon wavelength intervall

*/
void VAtmosAbsorption::read_extint( int LAMBDA_MIN )
{
    /* Program to read in extinction coefficients. Very finely tuned to file
       format in "kextint.dat"                                            */
    
    //   cout << "VAtmosAbsorption: reading atmospheric extinction file (lambda > " << LAMBDA_MIN << " nm): " << fSourceFile << endl;
    
    const int LAMBDA_MAX = 900;
    const int STEP_SIZE_NEEDED = 5;
    const int LSTEPS = ( 1 + ( LAMBDA_MAX - LAMBDA_MIN ) / STEP_SIZE_NEEDED );
    const int ALT_STEPS = 51;
    
    vector< double > i_yext( LSTEPS, 0. );
    for( int i = 0; i < ALT_STEPS; i++ )
    {
        extint.push_back( i_yext );
    }
    
    int    i, j, bin, mid_step, step, lambda, lambda_in;
    FILE*   finp;
    
    finp = fopen( fSourceFile.c_str(), "r" );
    if( finp != NULL )
    {
        for( lambda = LAMBDA_MIN, step = STEP_SIZE_NEEDED; lambda <= LAMBDA_MAX; lambda += step )
        {
            /* Read in wavelength. */
            if( fModel != "artemis" )
            {
                fscanf( finp, "%d %*[^\n]\n", &lambda_in );
            }
            else
            {
                fscanf( finp, "%d \n", &lambda_in );
            }
            
            if( lambda != lambda_in && fModel != "artemis" )
            {
                cerr << "VAtmosAbsorption::read_extint: Error in reading " << fSourceFile;
                cerr << " Lambda expected: " << lambda << "\t Lambda read " << lambda_in << endl;
            }
            
            bin = ( lambda - LAMBDA_MIN ) / STEP_SIZE_NEEDED;
            
            /* For given wavelength, read in  extinction coefficients.      */
            for( i = 0 ; i < ALT_STEPS ; i++ )
            {
                if( fModel != "artemis" )
                {
                    fscanf( finp, " %lf %*c", &extint[i][bin] );
                }
                // artemis file without ","
                else
                {
                    fscanf( finp, " %lf ", &extint[i][bin] );
                }
                if( !isnormal( extint[i][bin] ) && TMath::Abs( extint[i][bin] ) > 1.e-5 )
                {
                    cerr << "Warning: value not normal in " << fSourceFile << " for lambda = " << lambda << " and altitude " << i << " (setting value to 200.)" << endl;
                    extint[i][bin] = 200.;
                }
            }
            
            /* Fill in blanks in coefficients array, if any */
            if( step != STEP_SIZE_NEEDED )
            {
                mid_step = step / STEP_SIZE_NEEDED;
                
                for( i = 1 ; i < mid_step ; i++ )
                    for( j = 0 ; j < ALT_STEPS ; j++ )
                        extint[j][bin - mid_step + i] = extint[j][bin - mid_step] +
                                                        ( ( double ) i / mid_step ) *
                                                        ( extint[j][bin] -  extint[j][bin - mid_step] );
            }
            
            if( fModel != "artemis" )
            {
                switch( lambda )     /* These step sizes are specific to file
				     format in kextint.dat              */
                {
                    case 270:
                        step = 10;
                        break;
                        
                    case 280:
                        step = 20;
                        break;
                        
                    case 400:
                        step = 50;
                        break;
                        
                    case 700:
                        step = 100;
                        break;
                }
            }
        }
    }
    
    else
    {
        cout << "Unable to open input file: ";
        cout << fSourceFile << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
}

/*

    read extinction values from new files (provided by M.Daniel)
*/
void VAtmosAbsorption::read_extint_F2( int LAMBDA_MIN )
{
    /* Program to read in extinction coefficients. Very finely tuned to file
       format in "kextint.dat"                                            */
    
    cerr << "\t reading atmospheric extinction file (from kascade): " << fSourceFile << endl;
    
    const int LAMBDA_MAX = 900;
    const int STEP_SIZE_NEEDED = 5;
    const int LSTEPS = ( 1 + ( LAMBDA_MAX - LAMBDA_MIN ) / STEP_SIZE_NEEDED );
    const int ALT_STEPS = 51;
    
    vector< double > i_yext( LSTEPS, 0. );
    for( int i = 0; i < ALT_STEPS; i++ )
    {
        extint.push_back( i_yext );
    }
    
    ifstream is;
    is.open( fSourceFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VAtmosAbsorption::read_extint_F2, source file not found: " << fSourceFile << endl;
        exit( -1 );
    }
    
    string is_line;
    string temp;
    string temp2;
    
    int lambda, bin;
    
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            is_stream >> temp;
            lambda = atoi( temp.c_str() );
            bin = ( lambda - LAMBDA_MIN ) / STEP_SIZE_NEEDED;
            
            for( int i = 0; i < ALT_STEPS; i++ )
            {
                is_stream >> temp;
                extint[i][bin] = atof( temp.c_str() );
                if( !isnormal( extint[i][bin] ) )
                {
                    extint[i][bin] = 1.e10;
                }
            }
        }
    }
}


/*

    read extinction values from new files (Henrike)
*/
void VAtmosAbsorption::read_extint_M5( )
{
    /* Program to read in extinction coefficients. File format has one row per lamba. First column is lambda in nm, assumed to be in steps of 5. next is transmission from 1 - 2 km, 2 - 3 km, etc. Max height is not fixed.
    */
    
    const int LAMBDA_MAX = 900;
    const int LAMBDA_MIN = 205;
    const int STEP_SIZE_NEEDED = 5;
    const int LSTEPS = ( 1 + ( LAMBDA_MAX - LAMBDA_MIN ) / STEP_SIZE_NEEDED );
    //const int ALT_STEPS = 51;
    
    //vector< double > i_yext( LSTEPS, 0. );
    //for( int i = 0; i < ALT_STEPS; i++ ) extint.push_back( i_yext );
    
    ifstream is;
    is.open( fSourceFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VAtmosAbsorption::read_extint_M5, source file not found: " << fSourceFile << endl;
        exit( -1 );
    }
    else
    {
        cerr << "VAtmosAbsorption::read_extint_M5, reading source file: " << fSourceFile << endl;
    }
    
    string is_line;
    string temp;
    string temp2;
    
    int lambda;
    unsigned int bin;
    
    vector<double> bla( LSTEPS, 0 );
    extint.push_back( bla ); //extinction from 0 to 0 km is 0
    extint.push_back( bla ); //extinction from 0 to 1 km is unknown->set 0
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            is_stream >> temp;
            lambda = atoi( temp.c_str() );
            bin = ( lambda - LAMBDA_MIN ) / STEP_SIZE_NEEDED;
            unsigned int i = 2;
            
            while( is_stream >> temp )
            {
                double ext = atof( temp.c_str() );
                double minus_ln = - log( ext );
                if( extint.size() < i + 1 )
                {
                    vector<double> a;
                    extint.push_back( a );
                }
                if( extint[i].size() < bin + 1 )
                {
                    extint[i].push_back( minus_ln + extint[i - 1][bin] );
                }
                else
                {
                    cerr << "VAtmosAbsorption::read_extint_M5: Overwriting extint[" << i << "][" << bin << "]" << endl;
                    extint[i][bin] = minus_ln + extint[i - 1][bin];
                }
                if( !isnormal( extint[i][bin] ) )
                {
                    extint[i][bin] = 1.e10;
                }
                i++;
            }
        }
    }
    
    if( false ) //debugging
    {
        for( unsigned int bin = 0; bin < extint[0].size(); bin++ )
        {
            printf( "%d,   !Nanometers\n", LAMBDA_MIN + STEP_SIZE_NEEDED * bin );
            for( unsigned int i = 0; i < extint.size(); i++ )
            {
                printf( "%.8f, ", extint[i][bin] );
                if( i % 10 == 9 )
                {
                    printf( "\n" );
                }
            }
            printf( "\n" );
        }
    }
    
}
