/*
Vertex fitting code
*/
#include <TMath.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include "VertexFitFast.h"
//
// Constructors
//
//
// Empty construction (to be used when adding tracks later with AddTrk() )
VertexFitFast::VertexFitFast()
{
	fNtr = 0;
	fRstart = -1.0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fCovCstInv.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
}
//
// Build from list of parameters and covariances
VertexFitFast::VertexFitFast(Int_t Ntr, TVectorD** trkPar, TMatrixDSym** trkCov)
{
	fNtr = Ntr;
	fRstart = -1.0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fCovCstInv.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
	//
	Bool_t Charged = kTRUE;
	for (Int_t i = 0; i < fNtr; i++)
	{
		fPar.push_back(new TVectorD(*trkPar[i]));
		fParNew.push_back(new TVectorD(*trkPar[i]));
		fCov.push_back(new TMatrixDSym(*trkCov[i]));
		fCovNew.push_back(new TMatrixDSym(*trkCov[i]));
		fCharged.push_back(Charged);
	}
	fChi2List.ResizeTo(fNtr);
	//
}
//
// Build from list of parameters and covariances with charge tag
VertexFitFast::VertexFitFast(Int_t Ntr, TVectorD** trkPar, TMatrixDSym** trkCov, Bool_t* Charged)
{
	fNtr = Ntr;
	fRstart = -1.0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fCovCstInv.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
	
	for (Int_t i = 0; i < fNtr; i++)
	{
		fPar.push_back(new TVectorD(*trkPar[i]));
		fParNew.push_back(new TVectorD(*trkPar[i]));
		fCov.push_back(new TMatrixDSym(*trkCov[i]));
		fCovNew.push_back(new TMatrixDSym(*trkCov[i]));
		fCharged.push_back(Charged[i]);
	}
	fChi2List.ResizeTo(fNtr);
}

// Build from ObsTrk list of tracks
VertexFitFast::VertexFitFast(Int_t Ntr, ObsTrk** track)
{
	fNtr = Ntr;
	fRstart = -1.0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fCovCstInv.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
	//
	fChi2List.ResizeTo(fNtr);
	Bool_t Charged = kTRUE;
	for (Int_t i = 0; i < fNtr; i++)
	{
		fPar.push_back(new TVectorD(track[i]->GetObsPar()));
		fParNew.push_back(new TVectorD(track[i]->GetObsPar()));
		fCov.push_back(new TMatrixDSym(track[i]->GetCov()));
		fCovNew.push_back(new TMatrixDSym(track[i]->GetCov()));
		fCharged.push_back(Charged);
	}
}
//
// Destructor
//
void VertexFitFast::ResetWrkArrays()
{
	Int_t N = (Int_t)fdi.size();
	if(N > 0){
		for (Int_t i = 0; i < N; i++)
		{
			if (fx0i[i])  { fx0i[i]->Clear();	delete fx0i[i]; }
			if (fai[i])   { fai[i]->Clear();	delete fai[i]; }
			if (fdi[i])   { fdi[i]->Clear();	delete fdi[i]; }
			if (fAti[i])  { fAti[i]->Clear();	delete fAti[i]; }
			if (fDi[i])   { fDi[i]->Clear();	delete fDi[i]; }
			if (fWi[i])   { fWi[i]->Clear();	delete fWi[i]; }
			if (fWinvi[i]){ fWinvi[i]->Clear();     delete fWinvi[i]; }
		}
		fa2i.clear();
		fx0i.clear();
		fai.clear();
		fdi.clear();
		fAti.clear();
		fDi.clear();
		fWi.clear();
		fWinvi.clear();
	}
}
VertexFitFast::~VertexFitFast()
{	
	fxCst.Clear();		
	fCovCst.Clear();
	fCovCstInv.Clear();
	fXv.Clear();		
	fcovXv.Clear();		
	fChi2List.Clear();
	//
	for (Int_t i = 0; i < fNtr; i++)
	{
		fPar[i]->Clear();	delete fPar[i];
		fParNew[i]->Clear();	delete fParNew[i];
		fCov[i]->Clear();	delete fCov[i];
		fCovNew[i]->Clear();	delete fCovNew[i];
	}	
	fPar.clear();
	fParNew.clear();
	fCov.clear();
	fCovNew.clear();		
	//
	ResetWrkArrays();
	ffi.clear();	
	fCharged.clear();
	fNtr = 0;
}

TVectorD VertexFitFast::Fill_x(TVectorD par, Double_t phi, Bool_t Charged)
{
	//
	// Calculate track 3D position for a given phase, phi
	//
	TVectorD x(3);
	TVector3 xt;
	if(Charged) xt = Xtrack(par, phi);
	else        xt = Xtrack_N(par,phi);
	for(Int_t i = 0; i<3; i++)x(i) = xt(i);
	return x;
}

void VertexFitFast::VtxFitNoSteer()
{
	// Initialize
	std::vector<TVectorD*> x0i;				// Tracks at ma
	std::vector<TVectorD*> ni;				// Track derivative wrt phase
	std::vector<TMatrixDSym*> Ci;				// Position error matrix at fixed phase
	std::vector<TVectorD*> wi;				// Ci*ni
	std::vector<Double_t> s_in;				// Starting phase
	
    // Track loop
	for (Int_t i = 0; i < fNtr; i++)
	{
		TVectorD par = *fPar[i];
		TMatrixDSym Cov = *fCov[i];
		Double_t s = 0.;
		// Case when starting radius is provided
		if(fRstart > TMath::Abs(par(0))){
			if(fCharged[i])s = 2.*TMath::ASin(par(2)*TMath::Sqrt((fRstart*fRstart-par(0)*par(0))/(1.+2.*par(2)*par(0))));
			else s = TMath::Sqrt(fRstart*fRstart-par(0)*par(0));
		}
		x0i.push_back(new TVectorD(Fill_x(par, s, fCharged[i])));
		TMatrixD A(3,5);
		if(fCharged[i]){
			ni.push_back(new TVectorD(derXds(par, s)));
			A = derXdPar(par, s);}
		else{
			ni.push_back(new TVectorD(derXds_N(par, s)));
			A = derXdPar_N(par, s);}
		Ci.push_back(new TMatrixDSym(Cov.Similarity(A)));
		TMatrixDSym Cinv = RegInv(*Ci[i]);
		wi.push_back(new TVectorD(Cinv * (*ni[i])));
		s_in.push_back(s);
	}
	//
	// Get fit vertex
	//
	TMatrixDSym D(3); D.Zero();
	TVectorD Dx(3); Dx.Zero();
	for (Int_t i = 0; i < fNtr; i++)
	{
		TMatrixDSym Cinv = RegInv(*Ci[i]);
		TMatrixDSym W(3);
		W.Rank1Update(*wi[i], 1. / Ci[i]->Similarity(*wi[i]));
		TMatrixDSym Dd = Cinv - W;
		D += Dd;
		Dx += Dd * (*x0i[i]);
	}
	if(fVtxCst){
		D  += fCovCstInv;
		Dx += fCovCstInv*fxCst;
	}
	fXv = RegInv(D) * Dx;
    //
    // Chi2 computation (stage 1 only)
    //
    fChi2 = 0.0;
    for (Int_t i = 0; i < fNtr; i++){
        // residual r_i = x0_i - xv
        TVectorD r = (*x0i[i]) - fXv;

        // same Dd as in the vertex construction
        TMatrixDSym Cinv = RegInv(*Ci[i]);
        TMatrixDSym W(3);
        W.Rank1Update(*wi[i], 1. / Ci[i]->Similarity(*wi[i]));
        TMatrixDSym Dd = Cinv - W;

        // chi2 += r^T * Dd * r
        double chi2 = r * (Dd * r);
        fChi2 += chi2;
        fChi2List(i) = chi2;
    }
	//
	// Get fit phases
	//
	for (Int_t i = 0; i < fNtr; i++){
		Double_t si = Dot(*wi[i], fXv - (*x0i[i])) / Ci[i]->Similarity(*wi[i]);
		ffi.push_back(si+s_in[i]);
	}
	//
	// Cleanup
	for (Int_t i = 0; i < fNtr; i++)
	{
		x0i[i]->Clear();	delete x0i[i];
		ni[i]->Clear();		delete ni[i];
		Ci[i]->Clear();		delete Ci[i];
		wi[i]->Clear();		delete wi[i];
	}
	x0i.clear();;
	ni.clear();
	Ci.clear(); 
	wi.clear();
	s_in.clear();
}

void VertexFitFast::VertexFitter()
{
	if (fNtr < 2 && !fVtxCst){
		std::cout << "VertexFitFast::VertexFitter - Method called with less than 2 tracks - Aborting " << std::endl;
		std::exit(1);
	}
	//
	// Vertex fit
	//
	// Initial variable definitions
	TVectorD x(3);
	TMatrixDSym covX(3);
	Double_t Chi2 = 0;
	VtxFitNoSteer();	// Fast vertex finder on first pass (set ffi and fXv)
	TVectorD x0 = fXv;
	fVtxDone = kTRUE;		// Set fit completion flag
}
//
// Return fit vertex
TVectorD VertexFitFast::GetVtx()
{
	if (!fVtxDone) VertexFitter();
	return fXv;
}
//
// Return fit vertex covariance
TMatrixDSym VertexFitFast::GetVtxCov()
{
	if (!fVtxDone) VertexFitter();
	return fcovXv;
}
//
// Return fit vertex chi2
Double_t VertexFitFast::GetVtxChi2()
{
	if (!fVtxDone) VertexFitter();
	return fChi2;
}
//
// Return array of chi2 contributions from each track
TVectorD VertexFitFast::GetVtxChi2List()
{
	if (!fVtxDone) VertexFitter();
	return fChi2List;
}
