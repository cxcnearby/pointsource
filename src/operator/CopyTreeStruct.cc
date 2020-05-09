#include "operator/CopyTreeStruct.h"
#include <iostream>

CopyTreeStruct::CopyTreeStruct() {}
CopyTreeStruct::CopyTreeStruct(const char *root_file_name)
    : template_name_(root_file_name) {}
CopyTreeStruct::CopyTreeStruct(const char *root_file_name,
                               const TString tree_name)
    : template_name_(root_file_name), treename_(tree_name) {}
CopyTreeStruct::~CopyTreeStruct() {}

int CopyTreeStruct::Init() {
  TFile *file = TFile::Open(template_name_, "READ");
  if (!file) {
    init_status_ = 1;
    return init_status_;
  };
  // TObjArray *tobja_keys = (TObjArray *)file->GetListOfKeys();
  // if (!tobja_keys) {
  //   init_status_ = 2;
  //   return init_status_;
  // }
  // treename_ = tobja_keys->At(0)->GetName(); // An alternative way to get tree
  // name.

  // treename_ = file->GetListOfKeys()->At(0)->GetName(); // if you want to
  // check
  // more tree or hists,
  // use an TIter to get
  // all of them.
  //   treename_ = "trec"; // TODO: Here I use a deliberate name!!!
  if (treename_.IsNull()) {
    treename_ = file->GetListOfKeys()->At(0)->GetName();
    if (treename_.IsNull()) {
      init_status_ = 2;
      return init_status_;
    }
  }
  TTree *tree = (TTree *)file->Get(treename_);
  TIter next(tree->GetListOfLeaves());
  while (TLeaf *leaf = (TLeaf *)next()) {
    LeafInfo thisleaf;
    thisleaf.leafname = leaf->GetName();
    thisleaf.datatype = leaf->GetTypeName();
    thisleaf.datanumber = leaf->GetNdata();
    if (thisleaf.datatype == "Int_t") {
      if (thisleaf.datanumber == 1) {
        thisleaf.leafword = thisleaf.leafname + "/I";
      } else {
        thisleaf.leafword = thisleaf.leafname + "[" +
                            std::to_string(thisleaf.datanumber) + "]/I";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype == "Long64_t") {
      if (thisleaf.datanumber == 1) {
        thisleaf.leafword = thisleaf.leafname + "/L";
      } else {
        thisleaf.leafword = thisleaf.leafname + "[" +
                            std::to_string(thisleaf.datanumber) + "]/L";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype == "Short_t") {
      if (thisleaf.datanumber == 1) {
        thisleaf.leafword = thisleaf.leafname + "/S";
      } else {
        thisleaf.leafword = thisleaf.leafname + "[" +
                            std::to_string(thisleaf.datanumber) + "]/S";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype == "Float_t") {
      if (thisleaf.datanumber == 1) {
        thisleaf.leafword = thisleaf.leafname + "/F";
      } else {
        thisleaf.leafword = thisleaf.leafname + "[" +
                            std::to_string(thisleaf.datanumber) + "]/F";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype == "Double_t") {
      if (thisleaf.datanumber == 1) {
        thisleaf.leafword = thisleaf.leafname + "/D";
      } else {
        thisleaf.leafword = thisleaf.leafname + "[" +
                            std::to_string(thisleaf.datanumber) + "]/D";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype == "Bool_t") {
      if (thisleaf.datanumber == 1) {
        thisleaf.leafword = thisleaf.leafname + "/O";
      } else {
        thisleaf.leafword = thisleaf.leafname + "[" +
                            std::to_string(thisleaf.datanumber) + "]/O";
      }
      init_status_ = 0;
    } else {
      init_status_ = 3;
    }
    mLeafInfo_[thisleaf.leafname] = thisleaf;
  }
  if (file)
    file->Close();
  return init_status_;
} // init_status: 0 == success; 1 == open file failed; 2 == getlistofkeys
  // failed; 3 == datatype not listed.

int CopyTreeStruct::SetReadTree(TChain *&chain) {
  while (init_status_ != 0) {
    if (init_status_ > 0) {
      setreadtree_status_ = 1;
      std::cerr << "CopyTreeStruct::SetReadTree: Init() failed! error code: "
                << init_status_ << std::endl;
      return setreadtree_status_;
    } else if (init_status_ < 0) {
      Init();
    }
  }
  chain = new TChain(treename_);
  for (auto &&v : mLeafInfo_) {
    if (v.second.datatype == "Int_t") {
      chain->SetBranchAddress(v.second.leafname, &v.second.data_int);
      setreadtree_status_ = 0;
    } else if (v.second.datatype == "Long64_t") {
      chain->SetBranchAddress(v.second.leafname, &v.second.data_long64);
      setreadtree_status_ = 0;
    } else if (v.second.datatype == "Short_t") {
      chain->SetBranchAddress(v.second.leafname, &v.second.data_short);
      setreadtree_status_ = 0;
    } else if (v.second.datatype == "Float_t") {
      chain->SetBranchAddress(v.second.leafname, &v.second.data_float);
      setreadtree_status_ = 0;
    } else if (v.second.datatype == "Double_t") {
      chain->SetBranchAddress(v.second.leafname, &v.second.data_double);
      setreadtree_status_ = 0;
    } else if (v.second.datatype == "Bool_t") {
      chain->SetBranchAddress(v.second.leafname, &v.second.data_bool);
      setreadtree_status_ = 0;
    } else {
      setreadtree_status_ = 2;
    }
  }
  return setreadtree_status_;
} // setreadtree_status_: 0 == success; 1 == init error; 2 == data type error.

int CopyTreeStruct::SetWriteTree(TTree *&tree) {
  while (init_status_ != 0) {
    if (init_status_ > 0) {
      setwritetree_status_ = 1;
      std::cerr << "CopyTreeStruct::SetWriteTree: Init() failed! error code: "
                << init_status_ << std::endl;
      return setwritetree_status_;
    } else if (init_status_ < 0) {
      Init();
    }
  }
  tree = new TTree(treename_, "an auto-generated tree based on template");
  for (auto &&v : mLeafInfo_) {
    if (v.second.datatype == "Int_t") {
      tree->Branch(v.second.leafname, v.second.data_int, v.second.leafword);
      setwritetree_status_ = 0;
    } else if (v.second.datatype == "Long64_t") {
      tree->Branch(v.second.leafname, v.second.data_long64, v.second.leafword);
      setwritetree_status_ = 0;
    } else if (v.second.datatype == "Short_t") {
      tree->Branch(v.second.leafname, v.second.data_short, v.second.leafword);
      setwritetree_status_ = 0;
    } else if (v.second.datatype == "Float_t") {
      tree->Branch(v.second.leafname, v.second.data_float, v.second.leafword);
      setwritetree_status_ = 0;
    } else if (v.second.datatype == "Double_t") {
      tree->Branch(v.second.leafname, v.second.data_double, v.second.leafword);
      setwritetree_status_ = 0;
    } else if (v.second.datatype == "Bool_t") {
      tree->Branch(v.second.leafname, v.second.data_bool, v.second.leafword);
      setwritetree_status_ = 0;
    } else {
      setwritetree_status_ = 2;
    }
  }
  return setwritetree_status_;
} // setwritetree_status_: 0 == success; 1 == init error; 2 == data type error.

TString CopyTreeStruct::GetTreeName() {
  while (init_status_ != 0) {
    if (init_status_ > 0) {
      std::cerr << "CopyTreeStruct::GetTreeName: Init() failed! error code: "
                << init_status_ << std::endl;
      exit(1);
    } else if (init_status_ < 0) {
      Init();
    }
  }
  return treename_;
}

bool CopyTreeStruct::IsExists(const TString leaf_name) {
  auto iter = mLeafInfo_.find(leaf_name);
  if (iter != mLeafInfo_.end()) {
    return true;
  } else {
    return false;
  }
}

const double CopyTreeStruct::operator[](const TString leaf_name) const {
  double value;
  if (mLeafInfo_.at(leaf_name).datatype == "Int_t") {
    value = mLeafInfo_.at(leaf_name).data_int[0];
  } else if (mLeafInfo_.at(leaf_name).datatype == "Long64_t") {
    value = mLeafInfo_.at(leaf_name).data_long64[0];
  } else if (mLeafInfo_.at(leaf_name).datatype == "Short_t") {
    value = mLeafInfo_.at(leaf_name).data_short[0];
  } else if (mLeafInfo_.at(leaf_name).datatype == "Float_t") {
    value = mLeafInfo_.at(leaf_name).data_float[0];
  } else if (mLeafInfo_.at(leaf_name).datatype == "Double_t") {
    value = mLeafInfo_.at(leaf_name).data_double[0];
  } else if (mLeafInfo_.at(leaf_name).datatype == "Bool_t") {
    value = mLeafInfo_.at(leaf_name).data_bool[0];
  }
  return value;
}
