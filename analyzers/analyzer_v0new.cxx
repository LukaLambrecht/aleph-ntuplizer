// Standalone two-tier V0 (Ks/Lambda) reconstruction module.
//
// Ported from the "two-tier V0 module" added in
// https://github.com/bistapf/Aleph/pull/1 (analyzer_v0new.h / analyzer_truth.h),
// a sibling ALEPH EDM4HEP ntuplizer project. That PR reports (on their own
// validation samples): Ks purity 76% -> 91% at equal efficiency, usable-Lambda
// yield x1.57 at purity 35% -> 87%, and mass resolution 5-10x narrower than the
// existing V0 finder (analyzer_svfinder.cxx::get_V0s).
//
// What changed vs a straight copy-paste:
// - analyzer_v0new.h's findV0s() calls a PATCHED fork of FCCAnalyses
//   (Apranikstar/FCCAnalyses, vendored as bistapf's own git submodule), whose
//   VertexFitterSimple::VertexFitter_Tk has 3 extra parameters (solenoidBz,
//   rescale_cm_mm, fast) not present in the vanilla FCCAnalyses this repo uses.
//   Tracing the actual source diff between the patched and vanilla versions
//   showed the patch's "rescale_cm_mm=false" branch (the one findV0s actually
//   uses) is EXACTLY what vanilla VertexFitterSimple::VertexFitter_Tk already
//   does: a single consistent fit, no second re-fit with rescaled units (that
//   defect -- which the patch's OTHER branch, rescale_cm_mm=true, exists to
//   replicate for comparison -- is specific to this repo's own hand-written
//   VertexFitterMod in analyzer_svfinder.cxx, not to vanilla FCCAnalyses).
//   The only genuine correction needed on top of the vanilla call is documented
//   at getV0Pair() below.
// - Only the 5 truth-FREE candidate accessors from analyzer_truth.h are ported
//   (candChi2, candDxyz, candP, candCosPointing, candVtxPos, here renamed with
//   a v0_ prefix for this repo's naming convention); the actual MC truth-matching
//   machinery in that file is not needed for production wiring and is not
//   ported (per the source PR's own README, --newV0 only ever uses those 5
//   accessors from analyzer_truth.h).
// - Everything else (cut constants, pointing/band helper functions, the
//   findV0s pair loop and tiered selection logic, the ML-pull accessors) is a
//   direct, unmodified port: it is pure math on FCCAnalysesVertex/FCCAnalysesV0
//   objects with no FCCAnalyses-version dependency.
//
// Track/vertex object compatibility with this repo (checked, not assumed):
// - Input tracks: this repo's SecondaryTracks (TrackTools::getSecondaryTracks)
//   are built from EFlowTrack_1 (TrackTools::getModifiedTrackStates), which
//   already applies the same D0/omega sign flip the source PR calls
//   "flipD0_copy'ed" -- same sign convention, no adaptation needed.
// - Primary vertex: this repo's PrimaryVertexObject (PrimaryVertexTools::
//   fitRecoPrimaryVertex) is already a FCCAnalyses::VertexingUtils::
//   FCCAnalysesVertex, built via the same vanilla VertexFitterSimple::
//   VertexFitter_Tk -- exactly the PV type/source findV0s expects.
// - Output type: FCCAnalyses::VertexingUtils::FCCAnalysesV0, the same type
//   this repo's own get_V0s (analyzer_svfinder.cxx) already returns, so all
//   existing V0Candidates_* derivations in analysis.py (dxy, dxyz, cosPointing,
//   correctedMass, prel, thetarel, phirel, chi2, ...) work unchanged regardless
//   of which finder produced the collection.
//
// Units (checked empirically against this repo, not just read off comments):
// - Vertex position and chi2 need NO correction: verified bit-for-bit identical
//   between vanilla VertexFitterSimple::VertexFitter_Tk and this repo's own
//   VertexFitterMod for the same track pairs, on real ALEPH MC. Both apply the
//   same (no-op) Units_mm=true passthrough to the track parameters AND their
//   covariance, so position/chi2 stay self-consistently numeric-cm throughout
//   (mislabeled "mm" by various FCCAnalyses/analyzer comments, but harmless --
//   same conclusion as the DStarCandidates_*vtx_dxyz unit note in
//   studies/mesons/plot_dstar.py).
// - Momentum DOES need a correction, applied once in getV0Pair() below: see
//   that function's docstring for the empirically-verified factor (7.5,
//   confirmed via ~1000 real secondary-track pairs, median ratio 7.4996-7.5001
//   against this repo's own VertexFitterMod).

#include <algorithm>
#include <cmath>
#include <numeric>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"


#ifndef V0NewFinder_H
#define V0NewFinder_H

namespace V0NewFinder{

const double m_pi_ = 0.13957039;
const double m_p_  = 0.93827208;
const double MKS   = 0.497611;
const double MLAM  = 1.115683;

// ---------------------------------------------------------------------------
// Cut package: single named source. TIGHT = the adopted physics package;
// findV0s's parameter defaults reference these, and candTight re-evaluates
// them offline. LOOSE = the ML-training superset tier: flat pointing, widened
// AP bands, relaxed Lambda qT veto. Mass windows, chi2 and displacement window
// are COMMON to both tiers. (Unchanged from the source PR.)
// ---------------------------------------------------------------------------
constexpr double KS_M_LO = 0.40, KS_M_HI = 0.60;
constexpr double LAM_M_LO = 1.08, LAM_M_HI = 1.20;
constexpr double DIS_LO = 0.1, DIS_HI = 150.;
constexpr double CHI2_CUT = 10.;
constexpr double TIGHT_COS_KS_LOWP = 0.999, TIGHT_COS_KS_MIDP = 0.9995,
                 TIGHT_COS_KS_HIGHP = 0.9999;
// Lambda tight pointing: p-tiered -- {x1, x2, x2} loosening in 1-cos relative
// to the earlier flat cut. Deliberately NOT a mirror of the Ks ladder (see
// source PR: a wider mid tier measured worse purity in a heavy-flavour
// closure study).
constexpr double TIGHT_COS_LAM_LOWP = 0.99995, TIGHT_COS_LAM_MIDP = 0.9999,
                 TIGHT_COS_LAM_HIGHP = 0.9999;
constexpr double TIGHT_QT_MIN_LAM = 0.04;
constexpr double AP_BAND_KS = 0.05, AP_LAM_LO = 0.10, AP_LAM_HI = 0.20;
constexpr double TIGHT_NSIG_KS_LOWP = 3., TIGHT_NSIG_KS_HIGHP = 4.;
constexpr double LOOSE_COS_POINT = 0.999;
constexpr double LOOSE_QT_MIN_LAM = 0.02;
constexpr double LOOSE_NSIG_KS = 6.;
constexpr double LOOSE_LAM_BAND_LO = 0.20, LOOSE_LAM_BAND_HI = 0.40;
// widened loose Lambda band for the tail-measurement variant (2x nominal)
constexpr double WIDE_LAM_BAND_LO = 2. * LOOSE_LAM_BAND_LO,
                 WIDE_LAM_BAND_HI = 2. * LOOSE_LAM_BAND_HI;
constexpr double LAM_P_LO = 8., LAM_P_HI = 20.;
// Lambda AP-band ellipse resolution: quadrature fit to the measured 68% width
// of the band variable vs p (truth-matched Lambda, p>2.5, light-flavour MC,
// per the source PR). Note: "nsig" is in units of that 68% width.
constexpr double SIG_ELL_LAM_A = 0.01622, SIG_ELL_LAM_B = 0.0033748,
                 SIG_ELL_LAM_C = 0.00015544;
// Ks AP-band ellipse resolution, linear model -- single source for the
// tight/loose band cuts AND the bandSig pull.
constexpr double SIG_ELL_KS_A = 0.007, SIG_ELL_KS_B = 0.0015;
constexpr double TIGHT_LAM_NSIG = 3.;
// Mass resolution sigma_m(p) [GeV]: quadrature fits to the measured mass
// resolution vs p over the full mass window (p>2.5; binned-fit estimator,
// descriptive, per the source PR) -- used by candMassSig.
constexpr double SIG_M_KS_A = 2.658e-3, SIG_M_KS_B = 0.5214e-3,
                 SIG_M_KS_C = 0.01418e-3;
constexpr double SIG_M_LAM_A = 1.045e-3, SIG_M_LAM_B = 0.2357e-3,
                 SIG_M_LAM_C = 0.005511e-3;

// Shared per-hypothesis acceptance helpers, used by findV0s (both tiers) and
// candTight (offline re-evaluation of the booked hypothesis).
inline double ksPointThr(double pmag, double lowp, double midp, double highp) {
  return (pmag < 2.) ? lowp : (pmag < 4.) ? midp : highp;
}
// Lambda tight pointing tiers -- same tier boundaries as Ks, different low-p
// value; single source for findV0s AND candTight.
inline double lamPointThr(double pmag, double lowp = TIGHT_COS_LAM_LOWP,
                          double midp = TIGHT_COS_LAM_MIDP,
                          double highp = TIGHT_COS_LAM_HIGHP) {
  return (pmag < 2.) ? lowp : (pmag < 4.) ? midp : highp;
}
inline double ksBandEll(double alpha, double qt, double pmag) {
  const double PSTAR_K = 0.20582, ESTAR_K = 0.248806;
  double beta = pmag / std::sqrt(pmag * pmag + MKS * MKS);
  double amax = PSTAR_K / (beta * ESTAR_K);
  return std::sqrt(std::pow(alpha / amax, 2) + std::pow(qt / PSTAR_K, 2));
}
inline double sigmaEllKs(double pmag) {
  return SIG_ELL_KS_A + SIG_ELL_KS_B * pmag;
}
inline double ksBandThr(double pmag, double floor_, double nsig_lo, double nsig_hi) {
  // resolution-scaled width: sigma_ell ~ 0.007+0.0015p;
  // floor_ acts as the low-p floor (bit-identical below ~9.5 GeV at nsig=3)
  double nsig = (pmag < 15.) ? nsig_lo : nsig_hi;
  return std::max(floor_, nsig * sigmaEllKs(pmag));
}
inline double lamBandEll(double alpha, double qt, double pmag) {
  const double PSTAR_L = 0.1005, ALPHA0_L = 0.69157;
  double beta = pmag / std::sqrt(pmag * pmag + MLAM * MLAM);
  double amp = 2. * PSTAR_L / (beta * MLAM);
  return std::sqrt(std::pow((std::abs(alpha) - ALPHA0_L) / amp, 2) +
                   std::pow(qt / PSTAR_L, 2));
}
inline double lamBandThr(double pmag, double lo, double hi) {
  return (pmag < LAM_P_LO) ? lo : (pmag < LAM_P_HI) ? lo + (hi - lo) * (pmag - LAM_P_LO) / (LAM_P_HI - LAM_P_LO) : hi;
}
inline double sigmaEllLam(double pmag) {
  return std::sqrt(SIG_ELL_LAM_A * SIG_ELL_LAM_A +
                   std::pow(SIG_ELL_LAM_B * pmag, 2) +
                   std::pow(SIG_ELL_LAM_C * pmag * pmag, 2));
}
// TIGHT Lambda AP band: resolution-scaled, floored at the earlier fixed low-p
// edge and CAPPED at the loose storage edge (protects the tight package from
// depending on the loose-tier config; see source PR for the closure-study
// numbers behind this choice).
inline double lamBandThrTight(double pmag, double floor_ = AP_LAM_LO,
                              double nsig = TIGHT_LAM_NSIG) {
  return std::min(std::max(floor_, nsig * sigmaEllLam(pmag)),
                  lamBandThr(pmag, LOOSE_LAM_BAND_LO, LOOSE_LAM_BAND_HI));
}

// momenta of the two tracks at the fitted vertex (already corrected to the
// true GeV scale by getV0Pair, so every downstream consumer sees consistent
// values)
inline void pairMomenta(const FCCAnalyses::VertexingUtils::FCCAnalysesVertex& v,
                        TVector3& p1, TVector3& p2) {
  p1 = v.updated_track_momentum_at_vertex[0];
  p2 = v.updated_track_momentum_at_vertex[1];
}

inline double invMass(const TVector3& p1, double m1, const TVector3& p2, double m2) {
  double e1 = std::sqrt(p1.Mag2() + m1 * m1);
  double e2 = std::sqrt(p2.Mag2() + m2 * m2);
  TVector3 p = p1 + p2;
  double e = e1 + e2;
  return std::sqrt(std::max(0., e * e - p.Mag2()));
}

// ---------------------------------------------------------------------------
// Single-consistent-fit wrapper around vanilla FCCAnalyses::VertexFitterSimple::
// VertexFitter_Tk (the "alltracks" overload, so reco_ind gets filled).
//
// Momentum correction (see file docstring for the full unit investigation):
// vertex position and chi2 need no correction (verified bit-for-bit identical
// to this repo's own VertexFitterMod on real data), but the fitted momentum
// magnitude needs a factor of exactly 10 * (solenoidBz / 2.0):
//  - x10: the same cm-native-tracks-passed-through-as-if-mm convention that
//    leaves position untouched (a pure relabeling) does NOT leave momentum
//    untouched, because the fitted curvature -> momentum conversion is
//    length-scale dependent. Treating cm-valued curvature as if unscaled
//    ("mm") makes the fitted momentum magnitude come out 10x too small.
//  - x(solenoidBz/2.0): VertexMore (inside VertexFitter_Tk) hardcodes B=2T;
//    ALEPH's real field is 1.5T, so solenoidBz=1.5 gives the standard 0.75
//    correction. This repo's own VertexFitterMod independently discovered and
//    fixes the identical issue (analyzer_svfinder.cxx: "par[2] *= (2/1.5)"),
//    applied to the curvature parameter pre-fit rather than to the output
//    momentum post-fit -- mathematically equivalent, same factor.
// Empirically verified: median momentum ratio of 7.4996-7.5001 (== 10*0.75)
// between this wrapper (uncorrected) and this repo's own VertexFitterMod,
// across ~1000 real secondary-track pairs from ALEPH MC.
// ---------------------------------------------------------------------------
inline FCCAnalyses::VertexingUtils::FCCAnalysesVertex getV0Pair(
    const ROOT::VecOps::RVec<edm4hep::TrackState>& tr_pair,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& np_tracks,
    double solenoidBz) {
  auto v = FCCAnalyses::VertexFitterSimple::VertexFitter_Tk(
      0, tr_pair, np_tracks, false, 0., 0., 0., 0., 0., 0.);
  double momCorrection = 10. * (solenoidBz / 2.0);
  for (auto& tp : v.updated_track_momentum_at_vertex) tp *= momCorrection;
  return v;
}

// ---------------------------------------------------------------------------
// The finder.
//   np_tracks : non-primary trackstates (this repo's SecondaryTracks --
//               already D0/omega sign-flipped via TrackTools::getModifiedTrackStates)
//   PV        : fitted primary vertex (this repo's PrimaryVertexObject)
// TWO-TIER selection: every pair is evaluated against the TIGHT (adopted)
// package first; only tight-failing pairs enter the LOOSE training tier.
// Tight candidates claim tracks first, so filtering the output to tight
// candidates (candTight==1) reproduces the tight-only module output EXACTLY,
// as if the loose tier had never run. Returns candidates in claim order (tight
// block first, chi2 ascending within each tier); pdgAbs holds the single best
// hypothesis (310 or 3122); invM the mass under that hypothesis.
// ---------------------------------------------------------------------------
inline FCCAnalyses::VertexingUtils::FCCAnalysesV0 findV0s(
    const ROOT::VecOps::RVec<edm4hep::TrackState>& np_tracks,
    const FCCAnalyses::VertexingUtils::FCCAnalysesVertex& PV,
    double solenoidBz = 1.5,
    double cos_point_ks = TIGHT_COS_KS_HIGHP,
    double ks_m_lo = KS_M_LO, double ks_m_hi = KS_M_HI,
    double lam_m_lo = LAM_M_LO, double lam_m_hi = LAM_M_HI,
    double dis_lo = DIS_LO, double dis_hi = DIS_HI,
    double cos_point_lam = TIGHT_COS_LAM_LOWP,
    double qt_min_lam = TIGHT_QT_MIN_LAM,
    double cos_ks_lowp = TIGHT_COS_KS_LOWP,
    double cos_ks_midp = TIGHT_COS_KS_MIDP,
    double ap_band_ks = AP_BAND_KS,
    double ap_lam_lo = AP_LAM_LO, double ap_lam_hi = AP_LAM_HI,
    double chi2_cut = CHI2_CUT,
    double trk_chi2_cut = -1.,
    bool lam_point_ks_tiers = false,
    double loose_lam_lo = LOOSE_LAM_BAND_LO,
    double loose_lam_hi = LOOSE_LAM_BAND_HI) {

  FCCAnalyses::VertexingUtils::FCCAnalysesV0 result;
  const int nTr = np_tracks.size();
  if (nTr < 2) return result;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);

  struct Cand {
    FCCAnalyses::VertexingUtils::FCCAnalysesVertex vtx;
    int i, j;
    int pdg;       // best hypothesis
    double m;      // mass under best hypothesis
    double chi2;
    bool tight;    // passed the tight (adopted) package
  };
  std::vector<Cand> cands;

  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair(2);
  for (int i = 0; i < nTr - 1; ++i) {
    tr_pair[0] = np_tracks[i];
    for (int j = i + 1; j < nTr; ++j) {
      if (np_tracks[i].omega * np_tracks[j].omega > 0) continue; // same charge
      tr_pair[1] = np_tracks[j];

      auto v = getV0Pair(tr_pair, np_tracks, solenoidBz);
      if (v.updated_track_momentum_at_vertex.size() != 2) continue;
      double chi2 = v.vertex.chi2; // normalised, ndf=1
      if (chi2 >= chi2_cut || !(chi2 == chi2)) continue;
      if (trk_chi2_cut > 0 && v.reco_chi2.size() == 2 &&
          (v.reco_chi2[0] > trk_chi2_cut || v.reco_chi2[1] > trk_chi2_cut))
        continue;

      // displacement window (cm) + pointing
      TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
      TVector3 d = x - pv;
      double dis = d.Mag();
      if (dis < dis_lo || dis > dis_hi) continue;

      TVector3 p1, p2;
      pairMomenta(v, p1, p2);
      TVector3 p = p1 + p2;
      double pmag = p.Mag();
      if (pmag <= 0) continue;
      double cp = d.Dot(p) / (dis * pmag);
      double qt = p1.Cross(p.Unit()).Mag(); // Armenteros qT (same for either daughter)
      // Armenteros alpha: physical charge = -sign(omega) (flipped collection)
      double la = p1.Dot(p) / pmag, lb = p2.Dot(p) / pmag;
      double q1 = (np_tracks[i].omega < 0) ? 1. : -1.;
      double lplus = (q1 > 0) ? la : lb, lminus = (q1 > 0) ? lb : la;
      double alpha = (lplus + lminus != 0.) ? (lplus - lminus) / (lplus + lminus) : 0.;

      // hypothesis masses: Ks(pipi), Lambda(p pi) with proton = higher-|p| track
      // (in a Lambda decay the baryon carries most of the momentum)
      double mks = invMass(p1, m_pi_, p2, m_pi_);
      double mlam = (p1.Mag() > p2.Mag()) ? invMass(p1, m_p_, p2, m_pi_)
                                          : invMass(p1, m_pi_, p2, m_p_);

      // TIGHT (adopted) package first. Arbitration among the tight-passing
      // hypotheses only, so the tight subset is EXACTLY what the module would
      // output with the loose tier switched off.
      bool inWinKs = (mks > ks_m_lo && mks < ks_m_hi);
      bool inWinLam = (mlam > lam_m_lo && mlam < lam_m_hi);
      bool okKs = inWinKs && cp > ksPointThr(pmag, cos_ks_lowp, cos_ks_midp, cos_point_ks);
      if (okKs && ap_band_ks > 0)
        okKs = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
               ksBandThr(pmag, ap_band_ks, TIGHT_NSIG_KS_LOWP, TIGHT_NSIG_KS_HIGHP);
      double cos_lam_thr = lam_point_ks_tiers
          ? ksPointThr(pmag, cos_ks_lowp, cos_ks_midp, cos_point_ks)
          : lamPointThr(pmag, cos_point_lam);
      bool okLam = inWinLam && cp > cos_lam_thr && qt > qt_min_lam;
      if (okLam && ap_lam_lo > 0)
        okLam = std::abs(lamBandEll(alpha, qt, pmag) - 1.) <
                lamBandThrTight(pmag, ap_lam_lo);
      bool tight = okKs || okLam;
      if (!tight) {
        // LOOSE training tier: flat pointing, widened AP bands, relaxed
        // Lambda qT veto; windows/chi2/displacement common.
        okKs = inWinKs && cp > LOOSE_COS_POINT;
        if (okKs && ap_band_ks > 0)
          okKs = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
                 ksBandThr(pmag, ap_band_ks, LOOSE_NSIG_KS, LOOSE_NSIG_KS);
        okLam = inWinLam && cp > LOOSE_COS_POINT && qt > LOOSE_QT_MIN_LAM;
        if (okLam && ap_lam_lo > 0)
          okLam = std::abs(lamBandEll(alpha, qt, pmag) - 1.) <
                  lamBandThr(pmag, loose_lam_lo, loose_lam_hi);
        if (!okKs && !okLam) continue;
      }
      double dks = std::abs(mks - MKS) / (0.5 * (ks_m_hi - ks_m_lo));
      double dlam = std::abs(mlam - MLAM) / (0.5 * (lam_m_hi - lam_m_lo));
      int pdg; double m;
      if (okKs && (!okLam || dks <= dlam)) { pdg = 310; m = mks; }
      else                                  { pdg = 3122; m = mlam; }

      cands.push_back({v, i, j, pdg, m, chi2, tight});
    }
  }

  // quality-ranked global claiming: tight candidates claim first (preserving
  // the tight-only output), then loose; best chi2 first within a tier
  std::vector<size_t> order(cands.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (cands[a].tight != cands[b].tight) return cands[a].tight;
    return cands[a].chi2 < cands[b].chi2;
  });

  std::vector<bool> used(nTr, false);
  for (size_t k : order) {
    const Cand& c = cands[k];
    if (used[c.i] || used[c.j]) continue;
    used[c.i] = true;
    used[c.j] = true;
    result.vtx.push_back(c.vtx);
    result.pdgAbs.push_back(c.pdg);
    result.invM.push_back(c.m);
  }
  return result;
}

// ---------------------------------------------------------------------------
// Truth-free per-candidate diagnostics (work on data; reco_ind is filled by
// this module, momenta are already at the true GeV scale via getV0Pair).
// ---------------------------------------------------------------------------

inline ROOT::VecOps::RVec<float> candAlpha(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& secondaries) {
  ROOT::VecOps::RVec<float> out;
  for (const auto& v : v0s.vtx) {
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(-99.);
      continue;
    }
    TVector3 pa = v.updated_track_momentum_at_vertex[0];
    TVector3 pb = v.updated_track_momentum_at_vertex[1];
    TVector3 p = pa + pb;
    double la = pa.Dot(p) / p.Mag(), lb = pb.Dot(p) / p.Mag();
    double q1 = (secondaries[v.reco_ind[0]].omega < 0) ? 1. : -1.;
    double lp = (q1 > 0) ? la : lb, lm = (q1 > 0) ? lb : la;
    out.push_back((lp + lm != 0.) ? (lp - lm) / (lp + lm) : -99.);
  }
  return out;
}

// Offline tight-package flag: 1 if the candidate's BOOKED hypothesis passes
// the adopted tight package (Ks pointing tiers + resolution-scaled band;
// Lambda pointing tiers + 3sigma-capped AP band), 0 if it entered via the
// loose training tier. Uses the same shared helpers/constants as findV0s
// (single source), so this encodes exactly that package and any tighter
// selection remains re-derivable offline from the stored loose tier. Assumes
// the adopted package (no variant override in the production). Works on data.
inline ROOT::VecOps::RVec<int> candTight(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s,
    const FCCAnalyses::VertexingUtils::FCCAnalysesVertex& PV,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& secondaries) {
  ROOT::VecOps::RVec<int> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& v = v0s.vtx[c];
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(0);
      continue;
    }
    TVector3 p1 = v.updated_track_momentum_at_vertex[0];
    TVector3 p2 = v.updated_track_momentum_at_vertex[1];
    TVector3 p = p1 + p2;
    double pmag = p.Mag();
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 d = x - pv;
    double dis = d.Mag();
    if (pmag <= 0 || dis <= 0) { out.push_back(0); continue; }
    double cp = d.Dot(p) / (dis * pmag);
    double qt = p1.Cross(p.Unit()).Mag();
    double la = p1.Dot(p) / pmag, lb = p2.Dot(p) / pmag;
    double q1 = (secondaries[v.reco_ind[0]].omega < 0) ? 1. : -1.;
    double lplus = (q1 > 0) ? la : lb, lminus = (q1 > 0) ? lb : la;
    double alpha = (lplus + lminus != 0.) ? (lplus - lminus) / (lplus + lminus) : 0.;
    double m = v0s.invM[c];
    bool ok;
    if (v0s.pdgAbs[c] == 310) {
      ok = (m > KS_M_LO && m < KS_M_HI) &&
           cp > ksPointThr(pmag, TIGHT_COS_KS_LOWP, TIGHT_COS_KS_MIDP, TIGHT_COS_KS_HIGHP);
      if (ok)
        ok = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
             ksBandThr(pmag, AP_BAND_KS, TIGHT_NSIG_KS_LOWP, TIGHT_NSIG_KS_HIGHP);
    } else {
      ok = (m > LAM_M_LO && m < LAM_M_HI) && cp > lamPointThr(pmag) && qt > TIGHT_QT_MIN_LAM;
      if (ok)
        ok = std::abs(lamBandEll(alpha, qt, pmag) - 1.) <
             lamBandThrTight(pmag);
    }
    out.push_back(ok ? 1 : 0);
  }
  return out;
}

// ML-input pulls: the selection variables the module cuts on, expressed in
// resolution units for training. Both SIGNED; -999 = undefined candidate.
inline ROOT::VecOps::RVec<float> candBandSig(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& secondaries) {
  ROOT::VecOps::RVec<float> out;
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& v = v0s.vtx[c];
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(-999.);
      continue;
    }
    TVector3 p1 = v.updated_track_momentum_at_vertex[0];
    TVector3 p2 = v.updated_track_momentum_at_vertex[1];
    TVector3 p = p1 + p2;
    double pmag = p.Mag();
    if (pmag <= 0) { out.push_back(-999.); continue; }
    double qt = p1.Cross(p.Unit()).Mag();
    double la = p1.Dot(p) / pmag, lb = p2.Dot(p) / pmag;
    double q1 = (secondaries[v.reco_ind[0]].omega < 0) ? 1. : -1.;
    double lplus = (q1 > 0) ? la : lb, lminus = (q1 > 0) ? lb : la;
    double alpha = (lplus + lminus != 0.) ? (lplus - lminus) / (lplus + lminus) : 0.;
    if (v0s.pdgAbs[c] == 310)
      out.push_back((ksBandEll(alpha, qt, pmag) - 1.) / sigmaEllKs(pmag));
    else
      out.push_back((lamBandEll(alpha, qt, pmag) - 1.) / sigmaEllLam(pmag));
  }
  return out;
}

inline ROOT::VecOps::RVec<float> candMassSig(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s) {
  ROOT::VecOps::RVec<float> out;
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& v = v0s.vtx[c];
    if (v.updated_track_momentum_at_vertex.size() < 2) {
      out.push_back(-999.);
      continue;
    }
    TVector3 p = v.updated_track_momentum_at_vertex[0] +
                 v.updated_track_momentum_at_vertex[1];
    double pmag = p.Mag(), p2 = pmag * pmag;
    bool isKs = (v0s.pdgAbs[c] == 310);
    double a = isKs ? SIG_M_KS_A : SIG_M_LAM_A;
    double b = isKs ? SIG_M_KS_B : SIG_M_LAM_B;
    double cc = isKs ? SIG_M_KS_C : SIG_M_LAM_C;
    double sig = std::sqrt(a * a + b * b * p2 + cc * cc * p2 * p2);
    out.push_back((v0s.invM[c] - (isKs ? MKS : MLAM)) / sig);
  }
  return out;
}

// Pointing significance: chi2-like significance of the displacement component
// PERPENDICULAR to the candidate momentum, using candidate + PV position
// covariance (both self-consistently numeric-cm, see file docstring). A
// well-pointing candidate has sig ~ O(1) regardless of how precisely it is
// measured. Returns -1 if the transverse covariance is singular/non-positive.
inline ROOT::VecOps::RVec<float> candPointSig(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s,
    const FCCAnalyses::VertexingUtils::FCCAnalysesVertex& PV) {
  ROOT::VecOps::RVec<float> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  const auto& cP = PV.vertex.covMatrix; // packed lower-tri: xx,xy,yy,xz,yz,zz
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
    TVector3 d = x - pv;
    if (p.Mag() <= 0 || d.Mag() <= 0) { out.push_back(-1.); continue; }
    TVector3 ph = p.Unit();
    TVector3 u1 = ph.Orthogonal().Unit();
    TVector3 u2 = ph.Cross(u1);
    const auto& cV = v.vertex.covMatrix;
    double C[3][3] = {
      {double(cV[0]) + cP[0], double(cV[1]) + cP[1], double(cV[3]) + cP[3]},
      {double(cV[1]) + cP[1], double(cV[2]) + cP[2], double(cV[4]) + cP[4]},
      {double(cV[3]) + cP[3], double(cV[4]) + cP[4], double(cV[5]) + cP[5]}};
    auto quad = [&](const TVector3& a, const TVector3& b) {
      double s = 0.;
      double av[3] = {a.X(), a.Y(), a.Z()}, bv[3] = {b.X(), b.Y(), b.Z()};
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) s += av[i] * C[i][j] * bv[j];
      return s;
    };
    double c11 = quad(u1, u1), c22 = quad(u2, u2), c12 = quad(u1, u2);
    double det = c11 * c22 - c12 * c12;
    if (det <= 0. || c11 <= 0. || c22 <= 0.) { out.push_back(-1.); continue; }
    double d1 = d.Dot(u1), d2 = d.Dot(u2);
    double sig2 = (d1 * (c22 * d1 - c12 * d2) + d2 * (c11 * d2 - c12 * d1)) / det;
    out.push_back(sig2 > 0. ? std::sqrt(sig2) : 0.);
  }
  return out;
}

inline ROOT::VecOps::RVec<float> candQt(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s) {
  ROOT::VecOps::RVec<float> out;
  for (const auto& v : v0s.vtx) {
    TVector3 pa = v.updated_track_momentum_at_vertex[0];
    TVector3 p = pa + v.updated_track_momentum_at_vertex[1];
    out.push_back(pa.Cross(p.Unit()).Mag());
  }
  return out;
}

// ---------------------------------------------------------------------------
// Truth-free candidate accessors ported from analyzer_truth.h (only the 5
// this module actually needs; see file docstring). Renamed with a v0_ prefix
// to avoid confusion with the cand* selection-variable accessors above.
// ---------------------------------------------------------------------------

inline ROOT::VecOps::RVec<float> v0_candChi2(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s) {
  ROOT::VecOps::RVec<float> out;
  for (const auto& v : v0s.vtx) out.push_back(v.vertex.chi2);
  return out;
}

inline ROOT::VecOps::RVec<float> v0_candDxyz(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s,
    const FCCAnalyses::VertexingUtils::FCCAnalysesVertex& PV) {
  ROOT::VecOps::RVec<float> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    out.push_back((x - pv).Mag());
  }
  return out;
}

inline ROOT::VecOps::RVec<float> v0_candP(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s) {
  ROOT::VecOps::RVec<float> out;
  for (const auto& v : v0s.vtx) {
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
    out.push_back(p.Mag());
  }
  return out;
}

inline ROOT::VecOps::RVec<float> v0_candCosPointing(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s,
    const FCCAnalyses::VertexingUtils::FCCAnalysesVertex& PV) {
  ROOT::VecOps::RVec<float> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
    TVector3 d = x - pv;
    out.push_back((d.Mag() > 0 && p.Mag() > 0) ? d.Dot(p) / (d.Mag() * p.Mag()) : -2.);
  }
  return out;
}

// Fitted-vertex position component (axis: 0=x, 1=y, 2=z).
inline ROOT::VecOps::RVec<float> v0_candVtxPos(
    const FCCAnalyses::VertexingUtils::FCCAnalysesV0& v0s, int axis) {
  ROOT::VecOps::RVec<float> out;
  for (const auto& v : v0s.vtx) out.push_back(v.vertex.position[axis]);
  return out;
}

}

#endif
