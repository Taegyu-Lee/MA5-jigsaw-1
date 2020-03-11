#include "SampleAnalyzer/User/Analyzer/ATLAS_SUSY_2018_06.h"
using namespace MA5;
using namespace std;

template<typename T1,typename T2> std::vector<const T1*>
Removal(std::vector<const T1*> &v1, std::vector<const T2*> &v2, const MAdouble64 &drmin);

template<typename T1,typename T2> std::vector<const T1*>
RemovalLeptons(std::vector<const T1*> &v1, std::vector<const T2*> &v2);

vector<int> LeptonPick(vector<const RecLeptonFormat*> leptons);
MAdouble64 GetHBoost(vector<const RecLeptonFormat*> leptons,const RecParticleFormat* missing);



// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool ATLAS_SUSY_2018_06::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;

  //INFO << "ATLAS" << endmsg;

  //Add Signal Region
  Manager()->AddRegionSelection("SR-low");
  Manager()->AddRegionSelection("SR-ISR");

  Manager()->AddCut("3Leptons");                              // Exact three leptons
  Manager()->AddCut("SFOS");                                  // Same flavor, opposite-charge pair
  //Manager()->AddCut("Dilepton Trigger");
  Manager()->AddCut("B-veto");                                // b-jet veto
  Manager()->AddCut("Three Lepton Mass");                     // Invariant mass of three leptons system
  Manager()->AddCut("low-LeptonPT > 60, 40, 30","SR-low");
  Manager()->AddCut("ISR-LeptonPT > 25, 25, 20","SR-ISR");
  Manager()->AddCut("Dilepton InvMass");                      // Invariant mass of SFOS pair

  Manager()->AddCut("low-Jet-veto","SR-low");
  Manager()->AddCut("low-HBoost","SR-low");
  Manager()->AddCut("low-PTsoft/(PTsoft+Meff)","SR-low");
  Manager()->AddCut("low-Meff/HBoost","SR-low");

  Manager()->AddCut("ISR-Njet > 0","SR-ISR");
  Manager()->AddCut("ISR-Njet < 4","SR-ISR");
  Manager()->AddCut("ISR-DeltaPhi(MET,Jets) > 2.0","SR-ISR");
  Manager()->AddCut("ISR-R(MET,Jets)","SR-ISR");
  Manager()->AddCut("ISR-JetsPT > 100","SR-ISR");
  Manager()->AddCut("ISR-MET > 80","SR-ISR");
  Manager()->AddCut("Transverse Mass");
  Manager()->AddCut("ISR-PTsoft < 25","SR-ISR");


  Manager()->AddHisto("SR-low-MT",20,0,500,"SR-low");
  Manager()->AddHisto("SR-low-HBoost",20,200,700,"SR-low");
  Manager()->AddHisto("SR-low-R(Meff,HBoost)",15,0.3,1.05,"SR-low");
  Manager()->AddHisto("SR-low-R(PTsoft,(PTsoft+Meff))",15,0,0.15,"SR-low");
  
  Manager()->AddHisto("SR-ISR-MT",50,0,500,"SR-ISR");
  Manager()->AddHisto("SR-ISR-R(MET,Jets)",18,0.1,1,"SR-ISR");
  Manager()->AddHisto("SR-ISR-PTsoft",20,0,100,"SR-ISR");
  Manager()->AddHisto("SR-ISR-PTjets",18,0,500,"SR-ISR");

  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void ATLAS_SUSY_2018_06::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;

  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool ATLAS_SUSY_2018_06::Execute(SampleFormat& sample, const EventFormat& event)
{
  //Check Weight
  double myWeight;
  if(Configuration().IsNoEventWeight()) myWeight=1. ;
  else if(event.mc()->weight()!=0.) myWeight=event.mc()->weight();
  else { return false; }
  Manager()->InitializeForNewEvent(myWeight);

  if(event.rec()==0) return true;
  //Objects Isolation
  //Signal Electrons
  vector<const RecLeptonFormat*> SignalElectrons;
  for(unsigned int i=0; i<event.rec()->electrons().size();i++)
  {
      const RecLeptonFormat *Lep = &(event.rec()->electrons()[i]);

      //Kinamatics
      double eta = Lep->abseta();
      double pt = Lep->pt();

      double iso_dR = 0.2;
      double iso_tracks = PHYSICS->Isol->eflow->relIsolation(Lep,event.rec(),iso_dR,0,IsolationEFlow::TRACK_COMPONENT);
      double iso_all = PHYSICS->Isol->eflow->relIsolation(Lep,event.rec(),iso_dR,0,IsolationEFlow::ALL_COMPONENTS);
      bool iso = (iso_tracks<0.06 && iso_all<0.06);
      if(eta<2.47 && pt>10. && iso) { SignalElectrons.push_back(Lep);}
  }
  //Signal Muon
  vector<const RecLeptonFormat*> SignalMuons;
  for(unsigned int i=0; i<event.rec()->muons().size();i++)
  {
      const RecLeptonFormat *Lep = &(event.rec()->muons()[i]);

      //Kinamatics
      double eta = Lep->abseta();
      double pt = Lep->pt();

      double iso_dR = 0.2;
      if(pt<=33) {iso_dR=0.3;}
      else if(pt<50) {iso_dR=-0.0059*pt + 0.4941;}

      double iso_tracks = PHYSICS->Isol->eflow->relIsolation(Lep,event.rec(),iso_dR,0,IsolationEFlow::TRACK_COMPONENT);
      double iso_all = PHYSICS->Isol->eflow->relIsolation(Lep,event.rec(),0.2,0,IsolationEFlow::ALL_COMPONENTS);
      bool iso = (iso_tracks<0.04 && iso_all<0.15);

      if(eta<2.4 && pt>10. && iso){SignalMuons.push_back(Lep);}
  }
  //Signal Jet
  vector<const RecJetFormat*> Jets;
  for(unsigned int i=0; i<event.rec()->jets().size();i++)
  {
    const RecJetFormat *Jet = &(event.rec()->jets()[i]);
    double eta = Jet->abseta();
    double pt = Jet->pt();
    if(pt>20 && eta<2.4) {Jets.push_back(Jet);}
  }

  // Overlap removal 
  SignalElectrons = Removal(SignalElectrons, SignalMuons,0.2);
  Jets = PHYSICS->Isol->JetCleaning(Jets,SignalElectrons,0.2);
  Jets = PHYSICS->Isol->JetCleaning(Jets,SignalMuons,0.2);
  SignalElectrons = RemovalLeptons(SignalElectrons,Jets);
  SignalMuons = RemovalLeptons(SignalMuons,Jets);



  vector<const RecLeptonFormat*> SignalLeptons;
  for(unsigned int i=0; i<SignalMuons.size();i++)
  {SignalLeptons.push_back(SignalMuons[i]);}
  for(unsigned int i=0; i<SignalElectrons.size();i++)
  {SignalLeptons.push_back(SignalElectrons[i]);}
  SORTER->sort(SignalLeptons,PTordering);

  //Missing Vector
  const RecParticleFormat *MET = &(event.rec()->MET());

  //////////////////////////////////////////////////////////////////////
  //Cut
  //////////////////////////////////////////////////////////////////////

  // Nleptons = 3
  bool Nlepton = (SignalLeptons.size()==3);
  if(!Manager()->ApplyCut(Nlepton,"3Leptons")) return true;
  
  //SFOS
  vector<int> PairLeptons = LeptonPick(SignalLeptons);
  bool SFOS = (SignalLeptons[PairLeptons[0]]->isMuon()==SignalLeptons[PairLeptons[1]]->isMuon());
  SFOS = SFOS && (SignalLeptons[PairLeptons[0]]->charge()!=SignalLeptons[PairLeptons[1]]->charge());
  if(!Manager()->ApplyCut(SFOS,"SFOS")) return true;

  //if(!Manager()->ApplyCut("Dilepton Trigger")); // Not Applied
  
  //B-veto
  bool Bveto = true;
  for(unsigned int i=0;i<Jets.size();i++)
  if(Jets[i]->btag()) {Bveto = false;}
  if(!Manager()->ApplyCut(Bveto,"B-veto")) return true;

  // Mlll
  MALorentzVector ThreeLeptonsSystem;
  ThreeLeptonsSystem = SignalLeptons[0]->momentum()+SignalLeptons[1]->momentum()+SignalLeptons[2]->momentum();
  MAdouble64 mlll = ThreeLeptonsSystem.M();
  if(!Manager()->ApplyCut(mlll>105.,"Three Lepton Mass")) return true;
  
  //Leptons PT
  bool LowLeptonsPT = SignalLeptons[0]->pt()>60. && SignalLeptons[1]->pt()>40. && SignalLeptons[2]->pt()>30.;
  bool ISRLeptonsPT = SignalLeptons[0]->pt()>25. && SignalLeptons[1]->pt()>25. && SignalLeptons[2]->pt()>20.;
  if(!Manager()->ApplyCut(LowLeptonsPT,"low-LeptonPT > 60, 40, 30")) return true;
  if(!Manager()->ApplyCut(ISRLeptonsPT,"ISR-LeptonPT > 25, 25, 20")) return true;

  //Mll
  MAdouble64 mll = (SignalLeptons[PairLeptons[0]]->momentum()+SignalLeptons[PairLeptons[1]]->momentum()).M();
  if(!Manager()->ApplyCut(75. <= mll && mll <= 105.,"Dilepton InvMass")) return true;

  ///////Low region///////
  
  //Number of jets
  int Jet_size = Jets.size();
  if(!Manager()->ApplyCut(Jet_size==0,"low-Jet-veto")) return true;

  //low HBoost
  MAdouble64 HBoost = GetHBoost(SignalLeptons,MET);
  Manager()->FillHisto("SR-low-HBoost",HBoost);
  if(!Manager()->ApplyCut(HBoost>250.,"low-HBoost")) return true;

  //low - PT soft
  MAdouble64 low_PTsoft = (ThreeLeptonsSystem+MET->momentum()).Pt();
  MAdouble64 meff = SignalLeptons[0]->pt()+SignalLeptons[1]->pt()+SignalLeptons[2]->pt()+MET->pt();
  MAdouble64 R_PTsoft  = low_PTsoft/(low_PTsoft+meff);
  Manager()->FillHisto("SR-low-R(PTsoft,(PTsoft+Meff))",R_PTsoft);
  if(!Manager()->ApplyCut(R_PTsoft < 0.05,"low-PTsoft/(PTsoft+Meff)")) return true;
  Manager()->FillHisto("SR-low-R(Meff,HBoost)",(meff/HBoost));
  if(!Manager()->ApplyCut(meff/HBoost > 0.9,"low-Meff/HBoost")) return true;
  
  ///////ISR region///////
  
  //Number of jets
  if(!Manager()->ApplyCut(Jet_size>0,"ISR-Njet > 0")) return true;
  if(!Manager()->ApplyCut(Jet_size<4,"ISR-Njet < 4")) return true;
  
  //ISR Delta Phi between MET and Jets
  MALorentzVector Jets_momentum;
  for(int i=0;i<Jet_size;i++)
  {Jets_momentum += Jets[i]->momentum();}
  double DeltaPhi_MET_Jets = abs((MET->momentum()).DeltaPhi(Jets_momentum));
  if(!Manager()->ApplyCut(DeltaPhi_MET_Jets > 2.0,"ISR-DeltaPhi(MET,Jets) > 2.0")) return true;

  //ISR Ratio PT betwwen MET and Jets
  double R_MET_Jets;
  double R_MET_Jets_x = MET->momentum().Px()*Jets_momentum.Px();
  double R_MET_Jets_y = MET->momentum().Py()*Jets_momentum.Py();
  R_MET_Jets = abs(R_MET_Jets_x+R_MET_Jets_y)/pow(Jets_momentum.Pt(),2);
  Manager()->FillHisto("SR-ISR-R(MET,Jets)",R_MET_Jets);
  if(!Manager()->ApplyCut(0.55<=R_MET_Jets && R_MET_Jets <= 1.,"ISR-R(MET,Jets)")) return true;

  //PT jet
  Manager()->FillHisto("SR-ISR-PTjets",Jets_momentum.Pt());
  bool ISR_JetPT = (Jets_momentum.Pt() > 100.);
  if(!Manager()->ApplyCut(ISR_JetPT,"ISR-JetsPT > 100")) return true;

  //MET PT
  if(!Manager()->ApplyCut(MET->pt()>80.,"ISR-MET > 80")) return true;


  // MT
  const RecLeptonFormat* thirdlepton = SignalLeptons[PairLeptons[2]];
  double MT = sqrt(2*thirdlepton->pt()*MET->pt()*(1-cos(thirdlepton->dphi_0_pi(MET->momentum()))));
  Manager()->FillHisto("SR-low-MT",MT);
  Manager()->FillHisto("SR-ISR-MT",MT);
  if(!Manager()->ApplyCut(MT>100.,"Transverse Mass")) return true;

  //ISR PT soft
  MAdouble64 ISR_PTsoft = (ThreeLeptonsSystem+Jets_momentum+MET->momentum()).Pt();
  Manager()->FillHisto("SR-ISR-PTsoft",ISR_PTsoft);
  if(!Manager()->ApplyCut(ISR_PTsoft<25.,"ISR-PTsoft < 25")) return true;


  return true; // end of anlaysis
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
template<typename T1,typename T2> std::vector<const T1*>
Removal(std::vector<const T1*> &v1, std::vector<const T2*> &v2, const MAdouble64 &drmin)
{
    std::vector<bool> mask(v1.size(),false);
    for(MAuint32 j=0;j<v1.size();j++)
      {for(MAuint32 i=0;i<v2.size();i++)
        if(v2[i]->dr(v1[j]) < drmin)
        {
          mask[j]=true;
          break;
        }
      }
    std::vector<const T1*> cleaned_v1;
    for(MAuint32 i=0;i<v1.size();i++)
      {if(!mask[i]) cleaned_v1.push_back(v1[i]);}

    return cleaned_v1;
}


template<typename T1,typename T2> std::vector<const T1*>
RemovalLeptons(std::vector<const T1*> &v1, std::vector<const T2*> &v2)
{
  std::vector<bool> mask(v1.size(),false);
  for(MAuint32 j=0;j<v1.size();j++)
  { 
    double pt = v1[j]->pt();
    double drmin = 0.2;
    if(pt <= 25) drmin = 0.4;
    else if(pt <= 50) drmin = -0.008*pt+0.6 ;
    else drmin = 0.2;

    for(MAuint32 i=0;i<v2.size();i++)
      if(v2[i]->dr(v1[j]) < drmin)
      {
        mask[j]=true;
        break;
      }
  }

  std::vector<const T1*> cleaned_v1;
  for(MAuint32 i=0;i<v1.size();i++)
    if(!mask[i]) cleaned_v1.push_back(v1[i]);

  return cleaned_v1;
}

vector<int> LeptonPick(vector<const RecLeptonFormat*> leptons)
{
  int size = leptons.size();
  vector<int> pair(2);
  double Zmass=91.1876;
  double Invmass;
  double deltaMass = 1000000;
  for(int i=0;i<size-1;i++)
    for(int j=i+1;j<size;j++)
    {
      if(leptons[i]->isMuon()==leptons[j]->isMuon())
        if(leptons[i]->charge()!=leptons[j]->charge())
        {
          double temp_invmass = (leptons[i]->momentum()+leptons[j]->momentum()).M();
          if(abs(Zmass-temp_invmass)<deltaMass)
          {
            Invmass = temp_invmass;
            deltaMass = abs(Zmass-Invmass);
            pair[0]=i;
            pair[1]=j;
          }

        }
    }

  for(int i=0;i<size;i++)
    if(i!=pair[0] && i!=pair[1]) pair.push_back(i);

  return pair;
}

MAdouble64 GetHBoost(vector<const RecLeptonFormat*> leptons,const RecParticleFormat* missing)
{
  MAdouble64 HBoost=0;

  MALorentzVector LeptonsMomentum = leptons[0]->momentum()+leptons[1]->momentum()+leptons[2]->momentum();
  MAdouble64 Missing_Pz = (LeptonsMomentum.Pz()*missing->pt())/sqrt(pow(LeptonsMomentum.Pt(),2)+pow(LeptonsMomentum.M(),2));
  
  MAdouble64 Missing_E = sqrt(pow(missing->pt(),2)+pow(Missing_Pz,2));
  MALorentzVector MissingMomentum(missing->px(),missing->py(),Missing_Pz,Missing_E);
  //cout<<MissingMomentum.M()<<endl;
  MAdouble64 ETotal = LeptonsMomentum.E()+MissingMomentum.P();

  MAdouble64 Boost_x = (LeptonsMomentum.Px()+MissingMomentum.Px())/ETotal;
  MAdouble64 Boost_y = (LeptonsMomentum.Py()+MissingMomentum.Py())/ETotal;
  MAdouble64 Boost_z = (LeptonsMomentum.Pz()+MissingMomentum.Pz())/ETotal;

  MABoost Boost(-Boost_x,-Boost_y,-Boost_z);
  
  for(unsigned int i=0;i<leptons.size();i++)
  {
    MALorentzVector temp_mom = leptons[i]->momentum();
    Boost.boost(temp_mom);
    //HBoost +=(Boost*leptons[i]->momentum()).P();
    HBoost += temp_mom.P();
  }
  HBoost += (Boost*MissingMomentum).P();

  return HBoost;
}

