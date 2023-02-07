/*
    A spartial constraint tumour growth model
    Copyright (C) 2018 Timon Heide (timon.heide@icr.ac.uk)
                     & Kate Chkhaidze (Kate.Chkhaidze@icr.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef _DEBUG_
#define D(x) x
#else
#define D(x)
#endif

#include "extern_global_variables.h"
#include "Phylogeny.h"
#include "Cell.h"
#include "CellType.h"
#include <boost/random/poisson_distribution.hpp>
#include <string>  
#include <iostream> 
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
// [[Rcpp::depends(BH)]]
#include "Rcpp.h"
#pragma clang diagnostic pop


// Statics:
unsigned int PhylogenyNode::msNextId = 0;

// Constructors
PhylogenyNode::PhylogenyNode(Cell* pCell) //
  : mId(msNextId++),
    mpUp(0),
    mpLeft(0),
    mpRight(0),
    mNumMutsGeneration(0),
    mNumStagedCells(0)
  {
    mTypeId = pCell->Type()->Id();
    pCell->AssociatedNode(this);
  }

PhylogenyNode::PhylogenyNode(Cell* pCell, PhylogenyNode* pUp)
  : mId(msNextId++),
    mpUp(pUp),
    mpLeft(0),
    mpRight(0),
    mNumMutsGeneration(0),
    mNumStagedCells(0)
  {
    mTypeId = pCell->Type()->Id();
    pCell->AssociatedNode(this);
  }


//Destructor
PhylogenyNode::~PhylogenyNode() {
  D(Rcpp::Rcout << "PhylogenyNode::~PhylogenyNode()" << std::endl;)

  // Delete left node:
  if (mpLeft != 0) {
    delete mpLeft;
    mpLeft = 0;
  }

  // Delete right node:
  if (mpRight != 0) {
    delete mpRight;
    mpRight = 0;
  }
  
  // Clean pointers in the up node:
  if (mpUp != 0) {
    if (mpUp->LeftNode() == this) {
      mpUp->LeftNode(0);
    } else if (mpUp->RightNode() == this) {
      mpUp->RightNode(0);
    }
  }
}


// Getters:
unsigned int PhylogenyNode::Id() const {return mId;};
int PhylogenyNode::TypeId() const {return mTypeId;};
unsigned int PhylogenyNode::NumMutations() const {return mNumMutsGeneration;};
unsigned long PhylogenyNode::NumStagedCells() const {return mNumStagedCells;};
PhylogenyNode* PhylogenyNode::UpNode() const {return mpUp;};
PhylogenyNode* PhylogenyNode::LeftNode() const {return mpLeft;};
PhylogenyNode* PhylogenyNode::RightNode() const {return mpRight;};

std::vector <PhylogenyNode*> PhylogenyNode::NodeAncestry() {
  std::vector <PhylogenyNode*> result;
  PhylogenyNode *pCurrent = this;
  do {
    result.push_back(pCurrent);
  } while ((pCurrent = pCurrent->UpNode()) != 0);
  return result;
};

void PhylogenyNode::MutationsInAncestry(std::vector<std::string> &results) const {

  std::stringstream hex_node_id;  // node id to a hex
  std::stringstream hex_mutation_id; // mutation id to a hex
  hex_node_id << std::hex << mId;
  
  for (unsigned int i = 0; i < mNumMutsGeneration; i++) {
    hex_mutation_id.str("");
    hex_mutation_id << std::hex << i;
    std::string m_id = hex_node_id.str() + "X" + hex_mutation_id.str(); // concat above ids
    results.push_back(m_id);  // append mutation to output vectors:
  }
  
  if (mpUp) mpUp->MutationsInAncestry(results); // recursion along the nodes.
  
};

unsigned int PhylogenyNode::MutationNumberInAncestry() const{
  unsigned int n_total = mNumMutsGeneration;
  if (mpUp) n_total += mpUp->MutationNumberInAncestry();
  return n_total;
};

void PhylogenyNode::CollectStagedNodes(CellSample &sample) const {
  // Skip if there are no staged nodes at this level anymore:
  if (mNumStagedCells == 0)
    return;

  // Convert node id to a hex
  std::stringstream hex_node_id;
  hex_node_id << std::hex << mId;

  // Sample each mutation in current node once:
  for (unsigned i = 0; i < mNumMutsGeneration; i++) {
    // Convert mutation id to a hex string:
    std::stringstream hex_mutation_id;
    hex_mutation_id.str("");
    hex_mutation_id << std::hex << i;
    std::string full_mut_id = hex_node_id.str() + "X" + hex_mutation_id.str();

    // Append mutation information to output vectors:
    sample.AddMutation(full_mut_id, mTypeId, mNumStagedCells);
  }

  // Now parse the nodes at the next level if these nodes are not empty:
  if (mpLeft != 0)
    mpLeft->CollectStagedNodes(sample);

  if (mpRight != 0)
    mpRight->CollectStagedNodes(sample);

  return;
};

void PhylogenyNode::StagedNodesToStream(std::ostream& os) const {
  
  // Skip if there are no staged nodes at this level anymore:
  if (mNumStagedCells == 0)
    return;
  
  if (mpLeft == 0 && mpRight == 0) { // as leaf
    if (mLabel != "") { os << mLabel; } else { os << mId; }
    os << ":" << mNumMutsGeneration;
  } else { // as node
    bool putLeft = mpLeft != 0 && mpLeft->NumStagedCells() > 0;
    bool putRight = mpRight != 0 && mpRight->NumStagedCells() > 0;
    os << "(";
    if (putLeft) mpLeft->StagedNodesToStream(os);
    if (putLeft && putRight) os << ",";
    if (putRight) mpRight->StagedNodesToStream(os);
    os << ")" << mId << ":" << mNumMutsGeneration;
  }
}

PhylogenyNode* PhylogenyNode::GetTrimmingPosition() {

  if (mpLeft || mpRight || !mpUp) { // no ancestors and not root
    Rcpp::Rcout << 0 << std::endl;
    return 0;
  }
  
  PhylogenyNode* cNode = this;
  PhylogenyNode* nNode = mpUp;
  
  while (!(nNode->LeftNode() && nNode->RightNode())) {
    if (!nNode->UpNode()) break; // do not delete the root node. 
    cNode = nNode;
    nNode = nNode->UpNode();
  }

  return cNode;
}
  

// Setters:
void PhylogenyNode::AddNewMutations(unsigned int num_new_mutations) {
  mNumMutsGeneration += num_new_mutations;
}

void PhylogenyNode::SetMutations(unsigned int num_new_mutations) {
  mNumMutsGeneration = num_new_mutations;
}

void PhylogenyNode::StageNodeForSequencing(std::string label){
  mLabel = label;
  StageNodeForSequencing();
};

void PhylogenyNode::StageNodeForSequencing(){
  mNumStagedCells++;
  if (mpUp != 0) // Not root node
    mpUp->StageNodeForSequencing();
};

void PhylogenyNode::UnstageNodes() {
  mLabel = "";
  if (mNumStagedCells > 0) {
    if (mpLeft != 0) mpLeft->UnstageNodes();
    if (mpRight != 0) mpRight->UnstageNodes();
    mNumStagedCells = 0; 
  }
}

void PhylogenyNode::ScaleGenerationTree(double scale_factor, Rcpp::List fixed_vals, double dispersion) {
  
  std::string mId_str = std::to_string(mId);
  if (fixed_vals.containsElementNamed(mId_str.c_str())) {
    mNumMutsGeneration = fixed_vals[mId_str];
  } else {
    if (dispersion == 0.0) {
      // Sample number of mutations and update current node:
      double mean_rate = scale_factor * mNumMutsGeneration;
      boost::random::poisson_distribution<int> dist_mutations(mean_rate);
      mNumMutsGeneration = dist_mutations(rng);
    } else {
      // Sample number as sum of n negative_binomials:
      double mu = scale_factor ;
      double theta = 1 / dispersion;
      unsigned int n_muts = 0;
    
      boost::random::gamma_distribution<double> gamma_component(theta); // gamma component
      mNumMutsGeneration = mNumMutsGeneration * scale_factor * gamma_component(rng) / theta;
    }
  }
  
  // Call method on left and right of the tree, too:
  if (mpLeft != 0) mpLeft->ScaleGenerationTree(scale_factor, fixed_vals, dispersion);
  if (mpRight != 0) mpRight->ScaleGenerationTree(scale_factor, fixed_vals, dispersion);
}

void PhylogenyNode::UpNode(PhylogenyNode* pNode) {mpUp = pNode;};

void PhylogenyNode::LeftNode(PhylogenyNode* pNode) {mpLeft = pNode;};

void PhylogenyNode::RightNode(PhylogenyNode* pNode) {mpRight = pNode;};

void PhylogenyNode::TypeId(int NewType) {
  mTypeId = NewType;
  if (mpUp) mpUp->TypeId(0); // update the TypeId of all ancestral nodes to 0 (undertermined/mix type):
};

void PhylogenyNode::Branch(Cell* old_cell, Cell* new_cell) {
  
  //if (mpRight || mpLeft || old_cell->AssociatedNode() != this) 
  //  Rcpp::stop("Invalid branching of PhylogenyNode.");
  
  mpLeft = new PhylogenyNode(old_cell, this);
  
  if (new_cell) { 
    mpRight = new PhylogenyNode(new_cell, this); // new node to the right
  } 
  
};

