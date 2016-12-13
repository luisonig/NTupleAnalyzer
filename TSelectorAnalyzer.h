#ifndef TSELECTOR_ANALYZER_H
#define TSELECTOR_ANALYZER_H

#include <fastjet/ClusterSequence.hh>
#include "TSelectorMain.h"

class TSelectorReader;

class TSelectorAnalyzer : public TSelectorMain
{
 public :

  //--[ Overridable methods:

  int get_alphaspower() const { return Int_t(*input_alphaspower) - opt_extra_alphas; }

  int  Type();
  void Notify();
  void Init(const TSelectorReader* reader);
  bool Process();
  void SlaveBegin();
  void SlaveTerminate();

  TSelectorAnalyzer();
  ~TSelectorAnalyzer();

  //--] Overridable methods

  //--[ Analysis stuff:

  typedef std::vector<fastjet::PseudoJet> PseudoJetVector;
  fastjet::PseudoJet get_vec(int i) const;

  void TestAnalysis();
  void PrintEvent(const PseudoJetVector particles);

  //--] Analysis stuff:


  //--[ Member variables:
  int opt_extra_alphas;     // number of extra alphas powers

  //--] Member variables


  //--[ Reweighting variables:


  //--] Reweighting variables


  //--[ Counting events:
    Int_t event_prev_id;
    long  event_groups;
    double event_trials;
    bool new_event;
  //--] Counting events


};

#endif
