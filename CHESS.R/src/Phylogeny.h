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

#ifndef PHYLOGENY_H
#define PHYLOGENY_H

#include <vector>
#include "CellSample.h"
#include "Rcpp.h"
class Cell; // forward decleare




class PhylogenyNode {
  const unsigned int mId;
  static unsigned int msNextId; // Increments from one.

  PhylogenyNode* mpUp;
  PhylogenyNode* mpLeft;
  PhylogenyNode* mpRight;

  unsigned int mNumMutsGeneration;  
  int mTypeId;
  unsigned long mNumStagedCells;
  std::string mLabel;
  
  // Getters (private):
  unsigned int Id() const;
  int TypeId() const;
  unsigned int NumMutations() const;
  unsigned long NumStagedCells() const;
  std::vector <PhylogenyNode*> NodeAncestry();
  PhylogenyNode* UpNode() const;
  PhylogenyNode* LeftNode() const;
  PhylogenyNode* RightNode() const;
  
  // Setters (private):
  void UpNode(PhylogenyNode*);
  void LeftNode(PhylogenyNode*);
  void RightNode(PhylogenyNode*);

  
  public:
    
    // Constructors:
    PhylogenyNode(Cell*);
    PhylogenyNode(Cell*, PhylogenyNode*);

    //Destructor:
    ~PhylogenyNode();

    // Getters:
    void MutationsInAncestry(std::vector<std::string>&) const;
    unsigned int MutationNumberInAncestry() const;
    void CollectStagedNodes(CellSample&) const;
    void StagedNodesToStream(std::ostream& os) const;
    PhylogenyNode* GetTrimmingPosition();
      
    // Setters:
    void AddNewMutations(unsigned int);
    void SetMutations(unsigned int);
    void StageNodeForSequencing();
    void StageNodeForSequencing(std::string);
    void UnstageNodes();
    void ScaleGenerationTree(double, Rcpp::List, double);
    void TypeId(int);
    void Branch(Cell*, Cell*);
};

#endif // PHYLOGENY_H
