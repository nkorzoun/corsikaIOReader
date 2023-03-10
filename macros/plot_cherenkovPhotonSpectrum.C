/*=============================================================================
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
/*
    simple ROOT macro to plot average simulated Cherenkov spectrum

    requires corsikaIOreader to run with -histo option

*/

#include <string>


/*
 * read spectrum histograms from file and average over all events
 *
*/
TH1D* getCherenkovSpectrum( string ifile, string iHistoName = "hSLambda" )
{
    TFile *f = new TFile( ifile.c_str() );
    if( f->IsZombie() )
    {
        cout << "error while reading file " << ifile << endl;
	return 0;
    }

    TTree *t = (TTree*)f->Get( "tcors" );
    if( !t ) return 0;

    TH1D *hI = 0;
    TH1D *h = 0;
    unsigned int iCounter = 0;

    t->SetBranchAddress( iHistoName.c_str(), &hI );

    cout << "reading histogram " << iHistoName << " from " << ifile << endl;
    cout << "total number of showers " << t->GetEntries() << endl;
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );

	if( hI )
	{
	    if( i == 0 )
	    {
	       h = (TH1D*)hI->Clone( "hOnMirror" );
            }
	    h->Add( hI );
            iCounter++;
        }
     }

     if( iCounter > 0 )
     {
          h->Scale( 1./iCounter );
     }
     return h;
}
	        
/*
 * plot average Cherenkov photon spectrum
 * before and after application of atmospheric 
 * absorption/scattering
 * 
*/
void plot_CherenkovPhotonSpectrum( string iHistoRootFile )
{
    TH1D *hSimulatedSpectrum = getCherenkovSpectrum( iHistoRootFile, "hGLambda" );
    if( !hSimulatedSpectrum) return;
    TH1D *hSimulatedAbsorbedSpectrum = getCherenkovSpectrum( iHistoRootFile, "hSLambda" );
    if( !hSimulatedAbsorbedSpectrum ) return;


    TCanvas *cCherenkovSpectrum = new TCanvas( "cCherenkovSpectrum", "Cherenkov Photon Spectrum", 100, 100, 600, 600 );
    cCherenkovSpectrum->SetGridx( 0 );
    cCherenkovSpectrum->SetGridy( 0 );

    hSimulatedSpectrum->Rebin( 2 );
    hSimulatedSpectrum->SetTitle( "" );
    hSimulatedSpectrum->SetStats( 0 );
    hSimulatedSpectrum->SetLineColor( 1 );
    hSimulatedSpectrum->SetLineWidth( 3 );
    hSimulatedSpectrum->SetAxisRange( 200., 800. );
    hSimulatedSpectrum->Draw();

    hSimulatedAbsorbedSpectrum->Rebin( 2 );
    hSimulatedAbsorbedSpectrum->SetLineColor(2 );
    hSimulatedAbsorbedSpectrum->Draw( "same" );
}
