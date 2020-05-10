#ifndef POINTSOURCE_INCLUDE_REBASE_COPYTREESTRUCT_H_
#define POINTSOURCE_INCLUDE_REBASE_COPYTREESTRUCT_H_

#include "TLeaf.h"
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>

class CopyTreeStruct {
public:
  struct LeafInfo_ {
    TString leafname_;
    TString datatype_;
    int datanumber_;
    TString leafword_;
    Int_t data_int_[1000];
    Long64_t data_long64_[1000];
    Short_t data_short_[1000];
    Float_t data_float_[1000];
    Double_t data_double_[1000];
    Bool_t data_bool_[1000];
    // could add more!
  };
  std::vector<LeafInfo_> vLeafInfo_;
  CopyTreeStruct();
  CopyTreeStruct(const char *root_file_name);
  CopyTreeStruct(const char *root_file_name, const TString tree_name);
  ~CopyTreeStruct();
  int Init();
  int SetReadTree(TChain *&chain);
  int SetWriteTree(TTree *&tree);
  template <typename T> int SetBond(const char *leaf, T &var);
  // int SetBond(const char *leaf, T *&var);
  TString GetTreeName();
  //  public :
private:
  const char *template_name_;
  TString treename_;
  int init_status_ = -1;
  int setreadtree_status_ = -1;
  int setwritetree_status_ = -1;
};

template <typename T> int CopyTreeStruct::SetBond(const char *leaf, T &var) {
  // int CopyTreeStruct::SetBond(const char *leaf, T *&var) {
  int setbond_status_ = -1;
  while (init_status_ != 0) {
    if (init_status_ > 0) {
      setbond_status_ = 1;
      std::cerr << "CopyTreeStruct::SetBond: Init() failed! error code: "
                << init_status_ << std::endl;
      return setbond_status_;
    } else if (init_status_ < 0) {
      Init();
    }
  }
  for (auto &&v : vLeafInfo_) {
    if (leaf == v.leafname_) {
      if (v.datatype_ == "Int_t") {
        var.I = v.data_int_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Long64_t") {
        var.L = v.data_long64_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Short_t") {
        var.S = v.data_short_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Float_t") {
        var.F = v.data_float_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Double_t") {
        var.D = v.data_double_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Bool_t") {
        var.O = v.data_bool_;
        setbond_status_ = 0;
      } else
        setbond_status_ = 2;
      break;
    }
  }
  if (setbond_status_ == -1) {
    var.I = nullptr;
    var.L = nullptr;
    var.S = nullptr;
    var.F = nullptr;
    var.D = nullptr;
    var.O = nullptr;
  }
  return setbond_status_;
} // setbond_status_: 0 == success; 1 == init error; 2 == data type error; -1 ==
  // not set bond.

#endif // COPYTREESTRUCT_H_
