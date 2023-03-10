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
============================================================================= */
/*! \class VCORSIKARunheader
    \brief data class for CORSIKA run options

    \author Gernot Maier
*/

#include "VCORSIKARunheader.h"

// ClassImp(VCORSIKARunheader)

VCORSIKARunheader::VCORSIKARunheader()
{
    reset();
}

void VCORSIKARunheader::reset()
{
    runnumber = 0;
    production_date = 0;
    corsika_version = 0;
    
    particleID = 0;
    
    startingaltitude_gcm2 = 0.;
    tstart = 0.;
    
    observation_level_m = 0.;
    
    energy_slope = 0.;
    energy_min_GeV = 0.;
    energy_max_GeV = 0.;
    
    zenith_min_deg = 0.;
    zenith_max_deg = 0.;
    azimuth_min_deg = 0.;
    azimuth_max_deg = 0.;
    
    viewcone_min_deg = 0.;
    viewcone_max_deg = 0.;
    
    nscatt = 0;
    xscatt_m = 0;
    yscatt_m = 0;
    
    cherenkov_flag = 0;
    
    cherenkov_bunchsize = 0.;
    cherenkov_bandwidth_min_nm = 0.;
    cherenkov_bandwidth_max_nm = 0.;
    
    geomagneticfield_arrang_deg = 0.;
    geomagneticfield_x_muT = 0.;
    geomagneticfield_z_muT = 0.;
    
    hadronic_model_low = 0;
    hadronic_model_high = 0;
    hadronic_model_transition_energy_GeV = 0.;
    
    
}

void VCORSIKARunheader::printHeader( std::ostream& a )
{
    a << "RUN " << runnumber << endl;
    a << "DATE " << production_date << endl;
    a << "CORSIKAVERSION " << corsika_version << endl;
    a << "OBSLEVEL " << observation_level_m << endl;
    a << "PARTICLEID " << particleID << endl;
    a << "STARTING ALTITUDE " << startingaltitude_gcm2 << endl;
    a << "TSTART " << tstart << endl;
    a << "E_SLOPE " << energy_slope << endl;
    a << "E_MIN " << energy_min_GeV << endl;
    a << "E_MAX " << energy_max_GeV << endl;
    a << "ZENITH_MIN " << zenith_min_deg << endl;
    a << "ZENITH_MAX " << zenith_max_deg << endl;
    a << "AZIMUTH_MIN " << azimuth_min_deg << endl;
    a << "AZIMUTH_MAX " << azimuth_max_deg << endl;
    a << "VIEWCONE_MIN " << viewcone_min_deg << endl;
    a << "VIEWCONE_MAX " << viewcone_max_deg << endl;
    a << "NSCATT " << nscatt << endl;
    a << "XSCATT " << xscatt_m << endl;
    a << "YSCATT " << yscatt_m << endl;
    a << "CBUNCH " << cherenkov_bunchsize << endl;
    a << "C_WMIN " << cherenkov_bandwidth_min_nm << endl;
    a << "C_WMAX " << cherenkov_bandwidth_max_nm << endl;
    a << "GEO_A " << geomagneticfield_arrang_deg << endl;
    a << "GEO_X " << geomagneticfield_x_muT << endl;
    a << "GEO_Z " << geomagneticfield_z_muT << endl;
    a << "HAD_LOW " << hadronic_model_low << endl;
    a << "HAD_HIGH " << hadronic_model_high << endl;
    a << "HAD_TRANS " << hadronic_model_transition_energy_GeV << endl;
    
    // print CHERENKOV FLAG
    a << "CFLAG " << cherenkov_flag << endl;
    bitset<32> EVTH76 = cherenkov_flag;
    for( unsigned int i = 0; i < 10; i++ )
    {
        EVTH76.set( i, 0 );
    }
    unsigned int iATM_tab = EVTH76.to_ulong() / 1024;
    EVTH76 = cherenkov_flag - 1024 * iATM_tab;
    a << "CERENKOV " << EVTH76.test( 0 ) << endl;
    a << "IACT " << EVTH76.test( 1 ) << endl;
    a << "CEFFIC " << EVTH76.test( 2 ) << endl;
    a << "ATMEXT " << EVTH76.test( 3 ) << endl;
    a << "ATMEXT with refraction " << EVTH76.test( 4 ) << endl;
    a << "VOLUMEDET " << EVTH76.test( 5 ) << endl;
    a << "CURVED " << EVTH76.test( 6 ) << endl;
    a << "SLANT " << EVTH76.test( 8 ) << endl;
    a << "ATM table " << iATM_tab << endl;
}
