#ifndef EWKANA_NTUPLER_TEVENTINFO_HH
#define EWKANA_NTUPLER_TEVENTINFO_HH

#include <TObject.h>

namespace mithep 
{
  class TEventInfo : public TObject
  {
    public:
      TEventInfo(){}
      ~TEventInfo(){}

      UInt_t runNum; 			    // run number in data
      UInt_t evtNum; 			    // event number in data
      UInt_t lumiSec;			    // lumi section      
      UInt_t nTracks0;			    // number of reconstructed tracks in event
      UInt_t nLeptons0;			    // number of reconstructed muons/electrons in event
      UInt_t nCaloTowers0;  		    // number of reconstructed calorimeter towers in event
      UInt_t nPV0;                          // number of reconstructed primary vertices in event
      UInt_t nPUEvents;                     // number of pileup events.
      ULong_t triggerBits;		    // HLT trigger bits 
      UInt_t l1triggerBits;		    // L1 trigger bits 
      Float_t pvx, pvy, pvz;		    // primary vertex with the most associated tracks 
      Float_t bsx, bsy, bsz;		    // beamspot					  
      Float_t tcMEx, tcMEy, tcSumET;	    // track-corrected MET
      Float_t pfMEx, pfMEy, pfSumET;	    // particle flow MET
      Float_t pfTrackMEx, pfTrackMEy, pfTrackSumET; // particle flow track MET
      Float_t pfNeutralMEx, pfNeutralMEy, pfNeutralSumET; // particle flow neutral MET
      Float_t pfNeutralNoFwdMEx, pfNeutralNoFwdMEy, pfNeutralNoFwdSumET; // particle flow neutral MET, up to eta 2.5
      Float_t PileupEnergyDensity;
      Float_t PileupEnergyDensityHighEta;

      Bool_t dyVeto;                        // Drell-Yan veto for W's					 
      Bool_t VGammaEvent;                   // 
      Float_t eventweight;
      Bool_t hasGoodPV;                     // event has a good PV?
      Float_t ZeroJettinessPFCandidates;
      Float_t ZeroJettinessTracks;
      Float_t ZeroJettinessCalotowers;
      Float_t BeamThrustPFCandidates;
      Float_t BeamThrustTracks;
      Float_t BeamThrustCalotowers;
      UInt_t nPUMinusOne;
      UInt_t nPUPlusOne;

    ClassDef(TEventInfo,1)
  };
}
#endif
