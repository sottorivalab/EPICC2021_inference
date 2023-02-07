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

#ifndef UNIVERSE_H
#define UNIVERSE_H

// Forward declerations: ///////////////////////////////////////////////////////
//class Phylogeny_Node;
class Cell;
class CellType;
class PhylogenyNode;
class Shape;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include "boost/multi_array.hpp"
#include "CellSample.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
// [[Rcpp::depends(BH)]]
#include "Rcpp.h"
#pragma clang diagnostic pop

#define EVENT_SUBCLONE 'c'
#define EVENT_KILLREGROW 'k'
#define EVENT_SNAPSHOT 's'
#define EVENT_NON '0'


// event struct ////////////////////////////////////////////////////////////////

struct event { 
  
  long double time;
  char event_type;
  int event_idx;
  
  event(long double t, char type, int idx) 
    : time(t), 
      event_type(type), 
      event_idx(idx) {}
  
  event() 
    : time(-1*std::numeric_limits<long double>::infinity()),
      event_type(EVENT_NON),
      event_idx(-1) {}
  
  void print() const {
    Rcpp::Rcout << "time: " << time 
                << ", event_type: " << event_type 
                << ", event_idx: " << event_idx 
                << std::endl;
  }
  
  friend bool operator<(const event a, const event b){
    return a.time < b.time;
  }
  
  friend bool operator<(const event a, long double t){
    return a.time < t;
  }
  
  friend bool operator<=(const event a, long double t){
    return a.time <= t;
  }
};


class event_list {
  
  std::vector<event> events;
  bool sorted;
  int idx;
  
public:
  
  event_list() 
    : sorted(false),
      idx(-1) {}
  
  event next() {
    
    if (!contains_events()) {
      return event(); // undefined element
    }
    
    if (!sorted) {
      std::sort(events.begin() + idx + 1, events.end());
      sorted=true;
    }
    
    idx++;
    return events[idx];
  }
  
  event last() {
    if (idx < 0) {
      return event(); // non, undefined element
    }
    return events[idx];
  }
  
  bool contains_events() const {
    return events.size() > idx + 1;
  }
  
  void add(long double t, char type, int idx) {
    add(event(t, type, idx));
  }
  
  void add(event e) {
    
    if (this->last().time > e.time) {
      Rcpp::stop("Can't add element prior to time of last element.");
    }
    
    events.push_back(e);
    
    std::sort(events.begin() + idx + 1, events.end());
    sorted=true;
  }
  
  void print() const {
    Rcpp::Rcout << "################################" << std::endl;
    Rcpp::Rcout << "      Event list:" << std::endl;
    Rcpp::Rcout << "################################" << std::endl;
    
    for (int i = 0; i < events.size(); i++) {
      Rcpp::Rcout << "Event #" << i;
      if (i == idx) {  Rcpp::Rcout << " (last)"; }
      if (i == idx + 1) { Rcpp::Rcout << " (next)"; }
      Rcpp::Rcout << ": ";
      events[i].print();
    }
    
    Rcpp::Rcout << "################################" << std::endl << std::endl;
  }
  
  
};




// Universe ////////////////////////////////////////////////////////////////////

class Universe {
  
    // Parameters
    std::vector <int> mvSizeX;
    std::vector <int> mvSizeY;
    std::vector <int> mvSizeZ;
    std::vector <std::array<bool, 3>> mIsLimit;
    unsigned int mNumberOfClones;
    unsigned int mNumberOfSpaces;
    
    unsigned int mNumberClonalMutations;

    std::vector <double> mvMutationrates;
    std::vector <double> mvBirthrates;
    std::vector <double> mvDeathrates;
    std::vector <double> mvAggressions;
    std::vector <unsigned int> mvPushPower;
    std::vector <double> mvCloneStartTimes;
    std::vector <double> mvKillRegrowTime;
    std::vector <int> mvFathers;
    std::vector <unsigned int> mvUniverses;
    
    int mSeed;
    bool mExploreNeighborhood;
    bool mFlagRecordHistory;
    bool mFlagReturnGenerationTree;
    bool mAltPush;
    
    // Other variables:
    long double mTime;
    std::vector< boost::multi_array<Cell*, 3> > mvSpace; // the space
    std::vector<CellType*> mpTypes;
    std::vector<PhylogenyNode*> mpPhylogenies;
    bool mLimitReached;
    unsigned int mStopPopulationSize;
    event_list events;
    
    // Universe history:
    std::vector<long double> mvHistoryTime;
    std::vector< std::vector <unsigned long> > mvCellCounts;
    
    
    // Universe snapshots:
    int mSnapshotIndex;
    std::vector <double> mvSnapshotTimes;
    Rcpp::List mUniverseSnapshots;
    
    
    // Transformation events:
    std::vector<int> mvTransformX;
    std::vector<int> mvTransformY;
    std::vector<int> mvTransformZ;
    std::vector<long double> mvTransformTime;
    std::vector<int> mvTransformTypeFrom;
    std::vector<int> mvTransformTypeTo;
    std::vector<int> mvDistToEdge;
    
    
    
    // Functions /////////////////////
    
    // Getter functions:
    bool LimitIsReached() const;
    void Size(int, int&, int&, int&) const;
    class Cell* GetCell (int i, int x, int y, int z);
    class CellType* NextReactionType(long double*, int*);
    class CellSample TakeSample(int, Shape*);
    std::vector <CellSample> TakeMultipleSamples(int);
    std::vector <std::string> SingleCellTrees(int, Shape*);
    std::vector <std::string> SingleCellTrees(int, Shape*, bool);
    
    // Setter functions:
    void IncrementTimeBy(long double);
    void RecordToHistory();
    void RecordSnapshot();
    
  public:
    // Constructor:
    Universe(int, int, int, Rcpp::NumericMatrix, int, int);
    
    // Destructor:
    ~Universe();
    
    // Setter functions:
    void SetSnapshotTimes(Rcpp::NumericVector);
    void SetExploreLocally(bool);
    void SetRecordHistory(bool);
    void SetReturnGenerationTree(bool);
    void ScaleMutationTree(double, Rcpp::List, double);
    void SetClonalMutations(unsigned int);
    void SetStopPopulationSize(unsigned int);
    void SetAltPush(bool);
    
    // Sampling functions (getters):
    Rcpp::List TakeSampleBox(int,int,int, int,int,int, double,int,int,double, int, double);
    Rcpp::List SingleCellTree(std::vector<int>, std::vector<int>, std::vector<int>);
    Rcpp::List SingleCellTreeWithNodeIds(std::vector<int>, std::vector<int>, std::vector<int>);
    Rcpp::List SingleCellTreeInternal(std::vector<int>, std::vector<int>, std::vector<int>, bool);
    Rcpp::List SingleCellTreeN(int, int, int, int);
    
    
    // Other getters:
    long double Time() const;
    Rcpp::NumericVector GetTime() const;
    unsigned int GetCellType(int, int, int);
    unsigned int GetMutationBurden(int, int, int);
    Rcpp::StringVector GetMutations(int, int, int);
    Rcpp::DataFrame GetHistory() const;
    Rcpp::NumericVector GetSnapshotTimes() const;
    Rcpp::DataFrame GetTransformations() const;
    Rcpp::List GetSnapshots() const;
    std::vector <unsigned long> CellCountsPerType() const;
    unsigned long CellCounts() const;
    void FindClosestFreeGrid(int, int&, int&, int&, double) const;
    bool GetExploreLocally() const;
    bool GetRecordHistory() const;
    bool GetReturnGenerationTree() const;
    double FindClosestEdge(int, int&, int&, int&, double&, double) const;
    double GetDistanceToEdge(int, int &, int &, int &, double) const;
    double GetDistanceToEdgeRcpp(int, int &, int &, int &, double) const;
    double GetDistanceToEdgeRcppV2(int, int &, int &, int &, double) const;
    double GetDistanceToEdge(int, int &, int &, int &, double, double) const;
    double GetDistanceToEdge(int, int &, int &, int &, double, double, bool) const;
    unsigned int GetClonalMutations() const;
    unsigned int GetStopPopulationSize() const;
    bool GetAltPush() const;
    
    Rcpp::IntegerVector FindClosestEdgeRcpp(int, int, int, int) const;
      
    
    // Functions for simulations
    class Cell* RemoveCell (Cell*);
    class Cell* RemoveCell (int, int, int, int);
    std::vector< std::array<int, 3> > FreeNeighbours(int, int, int, int) const;
    void MarkLimitIsReached();
    bool InsertCell(int, int, int, int, Cell*);
    bool ContainsCell(int, int, int, int) const;
    bool ContainsCoordinate(int, int, int, int) const;
    bool TouchesHardLimit(int, int, int, int) const;
    bool PushCell(int, double, double, double, double, double, double, unsigned int);
    bool PushCell(int, int, int, int, int, int, int);

    
    // Output Functions:
    void PrintSimulationParameters() const;
      
    // Other functions:
    bool RunSimulation(bool);

};

#endif // UNIVERSE_H
