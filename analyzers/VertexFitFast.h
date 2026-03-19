// Faster version of standard Delphes VertexFit used in FCCAnalyses.

// The speedup is obtained by disabling the second stage where the track parameters are allowed to vary.
// This fast method might be suitable for quick potential vertex finding, but not for detailed analysis.

// Originally copied from here:
// https://github.com/delphes/delphes/blob/master/external/TrackCovariance/VertexFit.h

#ifndef G__VERTEXFITFAST_H
#define G__VERTEXFITFAST_H

#include <TMath.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <vector>
#include <iostream>

#include "TrackCovariance/TrkUtil.h" // from Delphes
#include "TrackCovariance/ObsTrk.h" // from Delphes

// Class for vertex fitting

class VertexFitFast: public TrkUtil
{

private:

	// Inputs
	Int_t fNtr;				// Number of tracks
	std::vector<TVectorD*> fPar;		// Input parameter array
	std::vector<TVectorD*> fParNew;		// Updated parameter array
	std::vector<TMatrixDSym*> fCov;		// Input parameter covariances
	std::vector<TMatrixDSym*> fCovNew;	// Updated parameter covariances
	std::vector<Bool_t>fCharged;		// Charge tag

	// Constraints
	Bool_t fVtxCst;				// Vertex constraint flag
	TVectorD fxCst;				// Constraint value
	TMatrixDSym fCovCst;			// Constraint 
	TMatrixDSym fCovCstInv;			// Inverse of constraint covariance
	
    // Results
	Bool_t fVtxDone;			// Flag vertex fit completed
	Double_t fRstart;			// Starting value of vertex radius (0 = none)
	TVectorD fXv;				// Found vertex
	TMatrixDSym fcovXv;			// Vertex covariance
	Double_t fChi2;				// Vertex fit Chi2
	TVectorD fChi2List;			// List of Chi2 contributions
	
    // Work arrays
	std::vector<Double_t> ffi;			// Fit phases
	std::vector<TVectorD*> fx0i;			// Track expansion points
	std::vector<TVectorD*> fai;			// dx/dphi
	std::vector<TVectorD*> fdi;			// x-shift
	std::vector<Double_t> fa2i;			// a'Wa
	std::vector<TMatrixD*> fAti;			// A transposed
	std::vector<TMatrixDSym*> fDi;			// W-WBW
	std::vector<TMatrixDSym*> fWi;			// (ACA')^-1
	std::vector<TMatrixDSym*> fWinvi;		// ACA'
	
    // Service routines
	void ResetWrkArrays();				// Clear work arrays
	TVectorD Fill_x(TVectorD par, Double_t phi, Bool_t Q);	// Track position at given phase
	void VtxFitNoSteer();				// Vertex fitter routine w/o parameter steering
	void VertexFitter();				// Vertex fitter routine wrapper

public:
	
    // Constructors
	VertexFitFast(); // Initialize waiting for tracks
	VertexFitFast(Int_t Ntr, ObsTrk** tracks); // Initialize with ObsTrk tracks
	VertexFitFast(Int_t Ntr, TVectorD** trkPar, TMatrixDSym** trkCov); // Initialize with parameters and covariances
	VertexFitFast(Int_t Ntr, TVectorD** trkPar, TMatrixDSym** trkCov, Bool_t* Charged);	// Initialize with parameters and covariances
	
	// Destructor
	~VertexFitFast();
	
    // Accessors also trigger calculations when needed
	Bool_t IsCharged(Int_t i) { return fCharged[i]; };
	Int_t GetNtrk() { return fNtr; };
	TMatrixDSym GetOldCov(Int_t i) { return *fCov[i]; }; // Input track covariance
	TVectorD GetVtx();
	TMatrixDSym GetVtxCov();
	Double_t GetVtxChi2();
	TVectorD GetVtxChi2List();
	
};

#endif
