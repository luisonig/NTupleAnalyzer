#include <fstream>
#include <sstream>
#include <string>
#ifndef NDEBUG
  #include <iostream>
#endif

#include "TSelectorReader.h"
#include "TSelectorWrite.h"


extern "C" void heft_OLP_Start(const char *fname, int *ierr);
extern "C" void heft_OLP_EvalSubProcess(int, double*, double, double*, double*);
extern "C" void heft_OLP_Option(const char * assignment, int* success);
extern "C" void heft_OLP_PrintParameter(const char * filename);

extern "C" void full_OLP_Start(const char *fname, int *ierr);
extern "C" void full_OLP_EvalSubProcess(int, double*, double, double*, double*);
extern "C" void full_OLP_Option(const char * assignment, int* success);
extern "C" void full_OLP_PrintParameter(const char * filename);

// --------------------------------------------------------------------------- //
// Selector
// --------------------------------------------------------------------------- //

TSelectorWrite::TSelectorWrite()
  : debug(0)
{
  std::cout<<"*************************************"<<std::endl;
  std::cout<<"*** REWEIGHTING SELECTOR SELECTED ***"<<std::endl;
  std::cout<<"*************************************"<<std::endl;
}

TSelectorWrite::~TSelectorWrite()
{
//
}

int TSelectorWrite::Type()
{
  return 2;
}

void TSelectorWrite::SetFileName(string folder, string fname, string suffix)
{
  std::cout<<folder<<std::endl;
  std::cout<<fname<<std::endl;
  outfilename=fname.substr(0, fname.size()-5)+suffix+".root";
  outfilename=folder+"/"+outfilename;
  std::cout<<"Reweighted events save in file: ";
  std::cout<<outfilename<<std::endl;
  return;
}

void TSelectorWrite::SetParameters(double nf, double mt, double mb, double mbms)
{
  // Input parameters are VEV, mZ and alpha:
  alpha = 1./128.802223294837; // As computed in Sherpa
  vev   = 246.;
  Nf    = nf; 
  mZ    = 91.1876;
  mT    = mt;
  mB    = mb;
  mBMS  = mbms;

  // Derived parameters:
  e2 = alpha*4.*pi;
  mW = sqrt(mZ*mZ/2.+sqrt(mZ*mZ*mZ*mZ/4.-pi*alpha*mZ*mZ*vev*vev));
  sw = sqrt(1-mW*mW/(mZ*mZ));

  return;
}

void TSelectorWrite::Init(const TSelectorReader* reader)
{
  input_ed_ntuples = &reader->ed_ntuples;
  input_id = &reader->id;
  input_nparticle = &reader->nparticle;
  input_ncount = &reader->ncount;
  input_px = &reader->px[0];
  input_py = &reader->py[0];
  input_pz = &reader->pz[0];
  input_E  = &reader->E[0];
  input_px_f = &reader->px_f[0];
  input_py_f = &reader->py_f[0];
  input_pz_f = &reader->pz_f[0];
  input_E_f  = &reader->E_f[0];
  input_alphas = &reader->alphas;
  input_kf = &reader->kf[0];
  input_ps_wgt = &reader->ps_wgt;
  input_weight = &reader->weight;
  input_weight2 = &reader->weight2;
  input_me_wgt = &reader->me_wgt;
  input_me_wgt2 = &reader->me_wgt2;
  input_x1 = &reader->x1;
  input_x2 = &reader->x2;
  input_x1p = &reader->x1p;
  input_x2p = &reader->x2p;
  input_id1 = &reader->id1;
  input_id2 = &reader->id2;
  input_fac_scale = &reader->fac_scale;
  input_ren_scale = &reader->ren_scale;
  input_nuwgt = &reader->nuwgt;
  input_usr_wgts = &reader->usr_wgts[0];
  input_alphaspower = &reader->alphaspower;
  input_part = &reader->part[0];

  outputfile  = new TFile(outfilename.c_str(),"recreate");

  n_fChain = new TTree("t3","Reweighted Ntuple");
  n_fChain->Branch("id", &n_id, "id/I");
  n_fChain->Branch("nparticle", &n_nparticle, "nparticle/I");
  n_fChain->Branch("ncount", &n_ncount, "ncount/I");
  n_fChain->Branch("px", n_px, "px[nparticle]/D");
  n_fChain->Branch("py", n_py, "py[nparticle]/D");
  n_fChain->Branch("pz", n_pz, "pz[nparticle]/D");
  n_fChain->Branch("E", n_E, "E[nparticle]/D");
  n_fChain->Branch("alphas", &n_alphas, "alphas/D");
  n_fChain->Branch("kf", n_kf, "kf[nparticle]/I");
  n_fChain->Branch("ps_wgt", &n_ps_wgt, "ps_wgt/D");
  n_fChain->Branch("weight", &n_weight, "weight/D");
  n_fChain->Branch("weight2", &n_weight2, "weight2/D");
  n_fChain->Branch("me_wgt", &n_me_wgt, "me_wtg/D");
  n_fChain->Branch("me_wgt2", &n_me_wgt2, "me_wtg2/D");
  n_fChain->Branch("x1", &n_x1, "x1/D");
  n_fChain->Branch("x2", &n_x2, "x2/D");
  n_fChain->Branch("x1p", &n_x1p, "x1p/D");
  n_fChain->Branch("x2p", &n_x2p, "x2p/D");
  n_fChain->Branch("id1", &n_id1, "id1/I");
  n_fChain->Branch("id2", &n_id2, "id2/I");
  n_fChain->Branch("id1p", &n_id1p, "id1p/I");
  n_fChain->Branch("id2p", &n_id2p, "id2p/I");
  n_fChain->Branch("fac_scale", &n_fac_scale, "fac_scale/D");
  n_fChain->Branch("ren_scale", &n_ren_scale, "ren_scale/D");
  n_fChain->Branch("nuwgt", &n_nuwgt, "nuwgt/I");
  n_fChain->Branch("usr_wgts", &n_usr_wgts, "usr_wgts[nuwgt]/D");
  n_fChain->Branch("alphasPower", &n_alphaspower, "alphasPower/B");
  n_fChain->Branch("part", n_part, "part[2]/C");

  InitializeOLP();

}


void TSelectorWrite::Notify()
{
  return;
}

bool TSelectorWrite::Process()
{

  PrepareEvent();
  Reweight();

  return true;
}

void TSelectorWrite::SlaveBegin()
{
  // pass
}

void TSelectorWrite::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed.

  outputfile->Write();
  outputfile->Close();
  
  FinalizeStat();
  //analysis->analysis_finalize(this);
}


void TSelectorWrite::PrepareEvent()
{
  if (event_prev_id != get_event_id()) {
    event_prev_id = get_event_id();
    new_event = 1;
    event_groups++;
  }
  else {
    new_event = 0;
  }
  
  return;
}


int TSelectorWrite::InitializeOLP()
{
  // Setup varibles:
  const string fname("OLE_order.olc");
  int ierr(0);

  // Transmit to GoSam and initialize it:
  stringstream ss;
  ss<<fixed<<setprecision(12)<<mW;
  const string mW_str = "mW=" + ss.str();

  ss.str(std::string());
  ss<<fixed<<setprecision(12)<<mZ;
  const string mZ_str = "mZ=" + ss.str();

  ss.str(std::string());
  ss<<fixed<<setprecision(12)<<sw;
  const string sw_str = "sw=" + ss.str();
  
  ss.str(std::string());
  ss<<fixed<<setprecision(12)<<mT;
  const string mT_str = "mT=" + ss.str();  
  
  ss.str(std::string());
  ss<<fixed<<setprecision(12)<<mB;
  const string mB_str = "mB=" + ss.str();   
  
  ss.str(std::string());
  ss<<fixed<<setprecision(12)<<mBMS;
  const string mBMS_str = "mBMS=" + ss.str();

  ss.str(std::string());
  ss<<fixed<<setprecision(12)<<Nf;
  const string Nf_str = "Nf=" + ss.str();


#ifndef DISABLE_OLP
  heft_OLP_Option(mZ_str.c_str(), &ierr);
  heft_OLP_Option(mW_str.c_str(), &ierr);
  heft_OLP_Option(mT_str.c_str(), &ierr);    
  heft_OLP_Option(mBMS_str.c_str(), &ierr);
  
  full_OLP_Option(Nf_str.c_str(), &ierr);  
  full_OLP_Option(mZ_str.c_str(), &ierr);
  full_OLP_Option(mW_str.c_str(), &ierr);
  full_OLP_Option(mT_str.c_str(), &ierr);    
  full_OLP_Option(mB_str.c_str(), &ierr); 
  full_OLP_Option(mBMS_str.c_str(), &ierr);
  
  heft_OLP_Start(fname.c_str(),&ierr);
  full_OLP_Start(fname.c_str(),&ierr);

  if (debug) {
    heft_OLP_PrintParameter("heft_parameters.dat");
    full_OLP_PrintParameter("full_parameters.dat");
  }

#else  // DISABLE_OLP
  std::cout << "Fake calls to OLP_Option and OLP_Start"<<std::endl;
#endif // DISABLE_OLP
  
  FillSubProcessMap_h1j();

  return ierr;
}

void TSelectorWrite::Reweight()
{
  PseudoJetVector particles;
  double Etot = 0.0;
  double momenta[(get_nparticle()+2)*5];

  for (Int_t j=0; j<get_nparticle(); j++) {
    Etot+=get_E(j);
  }

  fastjet::PseudoJet vec1 = fastjet::PseudoJet(0., 0., get_x1()*Etot/(get_x1()+get_x2()), get_x1()*Etot/(get_x1()+get_x2()));
  vec1.set_user_index(get_id1());
  fastjet::PseudoJet vec2 = fastjet::PseudoJet(0., 0.,-get_x2()*Etot/(get_x1()+get_x2()), get_x2()*Etot/(get_x1()+get_x2()));
  vec2.set_user_index(get_id2());
  particles.push_back(vec1);
  particles.push_back(vec2);

  // Create and fill particle kinematic arrays:
  for (Int_t i=0; i<get_nparticle(); i++){

    fastjet::PseudoJet vec = fastjet::PseudoJet(get_px(i), get_py(i), get_pz(i), get_E(i));
    vec.set_user_index(get_kf(i));
    particles.push_back(vec);
  }

  SubProcessto2 flavlist;
  flavlist.fill(0);

  for (unsigned i=0; i<particles.size(); i++){
    momenta[0+i*5] = particles[i].E();
    momenta[1+i*5] = particles[i].px();
    momenta[2+i*5] = particles[i].py();
    momenta[3+i*5] = particles[i].pz();
    momenta[4+i*5] = particles[i].m();

    flavlist[i] = particles[i].user_index();
  } 

  std::map<SubProcessto2, int>::iterator it;
  it = h1j_SubProcesses.find(flavlist);
  if ( it != h1j_SubProcesses.end()){
    if (debug) {
      PrintEvent(particles);
      std::cout<<"subprocess = "<<h1j_SubProcesses[flavlist]<<std::endl;
    }
  }
  else {
    std::cerr<<"ERROR SUBPROCESS NOT FOUND!\n---> "
	     <<flavlist[0]<<" "<<flavlist[1]<<" -> "
	     <<flavlist[2]<<" "<<flavlist[3]<<std::endl;
    return;
  }

  double heft_wgt = 0.;
  double full_wgt = 0.;
  double pdfs = 0.;
  double params[10] = {};
  double res_heft[4] = {};
  double res_full[4] = {};

  // alphas for OLP:
  params[0]=1.0;

#ifndef DISABLE_OLP  
  heft_OLP_EvalSubProcess(h1j_SubProcesses[flavlist], momenta, orig_ren_scale(), params, res_heft);
  full_OLP_EvalSubProcess(h1j_SubProcesses[flavlist], momenta, orig_ren_scale(), params, res_full);

  double gs = sqrt(orig_alphas()*4.0*pi);
  
  // res_heft[0] = double pole;
  // res_heft[1] = single pole;
  // res_heft[2] = finite part;
  // res_heft[3] = born;

  // res_full[0] = double pole;
  // res_full[1] = single pole;
  // res_full[2] = finite part;

  // for heft take Born
  heft_wgt = orig_ps_wgt()*res_heft[3]*pow(gs,6)*e2;
 
  // for full take finite part of loop induced
  full_wgt = orig_ps_wgt()*res_full[2]*(2.*pi)*pow(gs,6)*e2/64.0/pow(pi,4);
  
  if(debug){
    std::cout<<"------ CHECK -----"<<std::endl;
    std::cout<<"ratio heft/orig = "<<heft_wgt/orig_me_wgt()<<std::endl;
    std::cout<<"ratio heft/full = "<<heft_wgt/full_wgt<<std::endl;
    std::cout<<"ratio full/orig = "<<full_wgt/orig_me_wgt()<<std::endl;
    std::cout<<" "<<std::endl;
  }
  

#else  // DISABLE_OLP
  std::cout << "Fake calls to OLP_EvalSubProcess"<<std::endl;
  heft_wgt = orig_me_wgt();
  full_wgt = orig_me_wgt();
#endif // DISABLE_OLP

  int check = CheckPoint(heft_wgt, full_wgt);
  if (check != 0){
    std::cout<<"Unstable point, set to zero."<<std::endl;
  }

  CopyEvent(0);
  
  pdfs      = orig_weight()/orig_me_wgt();
  n_weight  = full_wgt*pdfs;
  n_weight2 = n_weight;
  n_me_wgt  = full_wgt;
  n_me_wgt2 = n_me_wgt;

  n_fChain->Fill();

  return;
}

int TSelectorWrite::CheckPoint(double heft_wgt, double full_wgt)
{
  if (heft_wgt/orig_me_wgt()-1.0 > 1e-5){
    full_wgt = 0.0;
    std::cout<<"Born ratio unstable"<<std::endl;
    return 1;
  }

  if (full_wgt != full_wgt) {
    full_wgt = 0.0;
    std::cout<<"NaN detected"<<std::endl;
    return 1;
  }
  
  if (abs(full_wgt) == std::numeric_limits<double>::infinity()) {
    full_wgt = 0.0;
    std::cout<<"INF detected"<<std::endl;       
    return 1;
  }
  
  if (abs(full_wgt/full_wgt) > 50.0) {
    full_wgt = 0.0;
    std::cout<<"K-factor check failed"<<std::endl;
    return 1;
  }  
  return 0;
}


void TSelectorWrite::CopyEvent(bool fillchain)
{

  n_id        = get_event_id();
  n_nparticle = get_nparticle();
  n_ncount    = get_ncount();
  
  for (int i=0; i< n_nparticle; i++){
  n_px[i] = get_px(i);
  n_py[i] = get_py(i);
  n_pz[i] = get_pz(i);
  n_E[i]  = get_E(i);
  n_kf[i] = get_kf(i);
}

  n_alphas  = orig_alphas();
  n_ps_wgt  = orig_ps_wgt();
  n_weight  = orig_weight();
  n_weight2 = orig_weight2();
  n_me_wgt  = orig_me_wgt();
  n_me_wgt2 = orig_me_wgt2();

  n_x1  = get_x1();
  n_x2  = get_x2();
  n_x1p = get_x1p();
  n_x2p = get_x2p();
  n_id1 = get_id1();
  n_id2 = get_id2();
  
  n_fac_scale = orig_fac_scale();
  n_ren_scale = orig_ren_scale();
  n_nuwgt = get_nuwgt();

  for (int i=0; i<n_nuwgt; i++){
  n_usr_wgts[i] = orig_usr_wgts(i);  //[nuwgt]
  }
  n_alphaspower = get_alphaspower();

  for (int i=0; i<2; i++){
  n_part[i] = get_part(i);
  }

  if (fillchain) n_fChain->Fill();

  return;
}


void TSelectorWrite::FinalizeStat()
{
  std::cout << std::endl;
  std::cout << "--------------------" 
	    << std::endl;
  std::cout << "Finalize: "
	    << event_groups << " event groups"
	    << std::endl;
  return;
}


void TSelectorWrite::PrintEvent(PseudoJetVector particles)
{
  cout.precision(15);
  cout.setf(ios::scientific, ios::floatfield);

  std::cout<<"--------------------\n";
  std::cout<<"proc = "
	   <<particles[0].user_index()<<" "<<particles[1].user_index()<<" -> ";
  for(unsigned i=2; i<particles.size(); i++){
    std::cout<<particles[i].user_index()<<" ";
  }
  std::cout<<std::endl;
  for(unsigned i=0; i<particles.size(); i++){
    std::cout<<particles[i].E() <<"\t"
	     <<particles[i].px()<<"\t"
	     <<particles[i].py()<<"\t"
	     <<particles[i].pz()<<";\t m="
	     <<particles[i].m()<<std::endl;
  }
}

void TSelectorWrite::TestReweighting()
{

  std::cout<<"==================== REWEIGHTING ====================="<<std::endl;
  PseudoJetVector particles;

  Double_t Etot = 0.0;

  for (Int_t j=0; j<get_nparticle(); j++) {
    Etot+=get_E(j);
  }

  fastjet::PseudoJet vec1 = fastjet::PseudoJet(0., 0., get_x1()*Etot/(get_x1()+get_x2()), get_x1()*Etot/(get_x1()+get_x2()));
  vec1.set_user_index(get_id1());
  fastjet::PseudoJet vec2 = fastjet::PseudoJet(0., 0.,-get_x2()*Etot/(get_x1()+get_x2()), get_x2()*Etot/(get_x1()+get_x2()));
  vec2.set_user_index(get_id2());
  particles.push_back(vec1);
  particles.push_back(vec2);

  // Create and fill particle kinematic arrays:
  for (Int_t i=0; i<get_nparticle(); i++){

    fastjet::PseudoJet vec = fastjet::PseudoJet(get_px(i), get_py(i), get_pz(i), get_E(i));
    vec.set_user_index(get_kf(i));
    particles.push_back(vec);
  }

  PrintEvent(particles);
}

void TSelectorWrite::FillSubProcessMap_h1j()
{
  // h1j_SubProcesses[{ 1, -1, 25, 21}] = 0;
  // h1j_SubProcesses[{ 1, 21, 25,  1}] = 1;
  // h1j_SubProcesses[{-1,  1, 25, 21}] = 2;
  // h1j_SubProcesses[{-1, 21, 25, -1}] = 3;
  // h1j_SubProcesses[{21,  1, 25,  1}] = 4;
  // h1j_SubProcesses[{21, -1, 25, -1}] = 5; 
  // h1j_SubProcesses[{ 4, -4, 25, 21}] = 6;
  // h1j_SubProcesses[{ 4, 21, 25,  4}] = 7;
  // h1j_SubProcesses[{-4,  4, 25, 21}] = 8;
  // h1j_SubProcesses[{-4, 21, 25, -4}] = 9;
  // h1j_SubProcesses[{21,  4, 25,  4}] = 10;
  // h1j_SubProcesses[{21, -4, 25, -4}] = 11;
  // h1j_SubProcesses[{ 2, -2, 25, 21}] = 12;
  // h1j_SubProcesses[{ 2, 21, 25,  2}] = 13;
  // h1j_SubProcesses[{-2,  2, 25, 21}] = 14;
  // h1j_SubProcesses[{-2, 21, 25, -2}] = 15;
  // h1j_SubProcesses[{21,  2, 25,  2}] = 16;
  // h1j_SubProcesses[{21, -2, 25, -2}] = 17;
  // h1j_SubProcesses[{21, 21, 25, 21}] = 18;
  // h1j_SubProcesses[{ 3, -3, 25, 21}] = 19;
  // h1j_SubProcesses[{ 3, 21, 25,  3}] = 20;
  // h1j_SubProcesses[{-3,  3, 25, 21}] = 21;
  // h1j_SubProcesses[{-3, 21, 25, -3}] = 22;
  // h1j_SubProcesses[{21,  3, 25,  3}] = 23;
  // h1j_SubProcesses[{21, -3, 25, -3}] = 24;
  // h1j_SubProcesses[{ 5, -5, 25, 21}] = 0;
  // h1j_SubProcesses[{ 5, 21, 25,  5}] = 1;
  // h1j_SubProcesses[{-5,  5, 25, 21}] = 2;
  // h1j_SubProcesses[{-5, 21, 25, -5}] = 3;
  // h1j_SubProcesses[{21,  5, 25,  5}] = 4;
  // h1j_SubProcesses[{21, -5, 25, -5}] = 5;

  h1j_SubProcesses[{ 1, -1, 25, 21}] = 0;
  h1j_SubProcesses[{ 1, 21, 25,  1}] = 1;
  h1j_SubProcesses[{-1,  1, 25, 21}] = 2;
  h1j_SubProcesses[{-1, 21, 25, -1}] = 3;
  h1j_SubProcesses[{21,  1, 25,  1}] = 4;
  h1j_SubProcesses[{21, -1, 25, -1}] = 5; 
  h1j_SubProcesses[{ 4, -4, 25, 21}] = 0;
  h1j_SubProcesses[{ 4, 21, 25,  4}] = 1;
  h1j_SubProcesses[{-4,  4, 25, 21}] = 2;
  h1j_SubProcesses[{-4, 21, 25, -4}] = 3;
  h1j_SubProcesses[{21,  4, 25,  4}] = 4;
  h1j_SubProcesses[{21, -4, 25, -4}] = 5;
  h1j_SubProcesses[{ 2, -2, 25, 21}] = 0;
  h1j_SubProcesses[{ 2, 21, 25,  2}] = 1;
  h1j_SubProcesses[{-2,  2, 25, 21}] = 2;
  h1j_SubProcesses[{-2, 21, 25, -2}] = 3;
  h1j_SubProcesses[{21,  2, 25,  2}] = 4;
  h1j_SubProcesses[{21, -2, 25, -2}] = 5;
  h1j_SubProcesses[{21, 21, 25, 21}] = 6;
  h1j_SubProcesses[{ 3, -3, 25, 21}] = 0;
  h1j_SubProcesses[{ 3, 21, 25,  3}] = 1;
  h1j_SubProcesses[{-3,  3, 25, 21}] = 2;
  h1j_SubProcesses[{-3, 21, 25, -3}] = 3;
  h1j_SubProcesses[{21,  3, 25,  3}] = 4;
  h1j_SubProcesses[{21, -3, 25, -3}] = 5;
  h1j_SubProcesses[{ 5, -5, 25, 21}] = 0;
  h1j_SubProcesses[{ 5, 21, 25,  5}] = 1;
  h1j_SubProcesses[{-5,  5, 25, 21}] = 2;
  h1j_SubProcesses[{-5, 21, 25, -5}] = 3;
  h1j_SubProcesses[{21,  5, 25,  5}] = 4;
  h1j_SubProcesses[{21, -5, 25, -5}] = 5;

  return; 
}
