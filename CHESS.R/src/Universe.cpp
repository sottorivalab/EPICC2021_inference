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


#define TICK_SOURCE_STARTTIME CellCounts() // (Time(): gillespie time, CellCounts(): cell numbers)
#define TICK_MIN_INCREMENT 1.0


#define CENTER(X) ((X + 1) / 2) - 1
#define KILL_REGROW_PERCENT 0.99
#define GROW_TO_CLOSEST_EDGE false

#define SPACE_OUPUT_FILE_PREFIX "space"
#define GENOT_OUPUT_FILE_PREFIX "space_generation"
#define TYPES_OUPUT_FILE_PREFIX "types"
#define PHYLO_OUPUT_FILE_PREFIX "phylo"
#define HISTORY_OUPUT_FILE_PREFIX "timestep_cellcount"
#define DIPLAY_FILE "space_image"
#define PARAMS_OUTPUT_FILE_PREFIX "sim_params"

#include "extern_global_variables.h"
#include "CellType.h"
#include "Shapes.h"
#include "Cell.h"
#include "Universe.h"
#include "Phylogeny.h"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <cmath>
#include <limits.h>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_real.hpp>
// [[Rcpp::depends(BH)]]
#include "Rcpp.h"


// Universe ////////////////////////////////////////////////////////////////////


//Constructor:
Universe::Universe(
  int size_x, 
  int size_y, 
  int size_z,
  Rcpp::NumericMatrix clone_params,
  int seed, 
  int clonal_mutations
):
  mvSizeX(1, size_x),
  mvSizeY(1, size_y),
  mvSizeZ(1, size_z),
  mNumberOfSpaces(1),
  mNumberClonalMutations(clonal_mutations),
  mvCloneStartTimes(0.0),
  mSeed(seed),
  mExploreNeighborhood(false),
  mFlagRecordHistory(false),
  mFlagReturnGenerationTree(false),
  mTime(0.0L),
  mLimitReached(false),
  mvCellCounts(clone_params.ncol()),
  mSnapshotIndex(0), 
  mStopPopulationSize(0)
{
  
  mNumberOfClones = clone_params.ncol();
  
  for (int i = 0; i < clone_params.ncol(); i++) {
    mvBirthrates.push_back(clone_params(0,i));
    mvDeathrates.push_back(clone_params(1,i));
    mvAggressions.push_back(clone_params(2,i));
    mvPushPower.push_back(clone_params(3,i));
    mvMutationrates.push_back(clone_params(4,i));
    mvCloneStartTimes.push_back(clone_params(5,i));
    mvKillRegrowTime.push_back(clone_params(6,i));
    mvFathers.push_back(clone_params(7,i));
    mvUniverses.push_back(0);
  }

  
  // -----------------------------------------
  // Input validation
  // -----------------------------------------
  
  if (!mNumberOfClones) { // n_clones > 0
    Rcpp::stop("At least one clone has to be defined.");
  }
    
  // check that start times are ordered
  double last_start_time = mvCloneStartTimes[0]; // at least one exists, tested above
  for (unsigned int i = 1; i < mNumberOfClones; i++) { 
    if (mvCloneStartTimes[i] < last_start_time) {
      Rcpp::stop("Clone start times have to be sorted.");
    }
    last_start_time = mvCloneStartTimes[i];
  }
  
  // first father is -1 (i.e. no other type)
  if (mvFathers[0] != -1){
    Rcpp::stop("Fathers of first clone has to be -1.");
  }
  
  if (mvCloneStartTimes[0] != 0.0){
    Rcpp::stop("Start time of first clone has to be 0.");
  }
    
  // check that fathers are defined when sub clones are introduced
  for (int i = 1; i < mNumberOfClones; i++) { 
    if (mvFathers[i] >= i || mvFathers[i] < -1){
      Rcpp::stop("Fathers of clones have to be introduced first!");
    }
  }
  
  // check rates and values are within bounds
  for (unsigned int i = 0; i < mNumberOfClones; i++) { 
    
    if (mvBirthrates[i] <= 0.0) {
      Rcpp::stop("Birthrate must be >0.");
    }
    
    if (mvDeathrates[i] < 0.0) {
      Rcpp::stop("Deathrate must be >=0.");
    }
    
    if (mvAggressions[i] < 0.0) {
      Rcpp::stop("Aggression must be >=0.");
    }
    
    if (mvMutationrates[i] < 0.0) {
      Rcpp::stop("Mutation rate must be >=0.");
    }
    
    if (mvCloneStartTimes[i] < 0.0) {
      Rcpp::stop("Start time rate must be >=0.");
    }
    
    if (mvKillRegrowTime[i] < 0.0) { //?
      Rcpp::stop("Kill regrowth time must be >=0.");
    }
  }
  
  
  // -----------------------------------------
  // Initialisation of remaining variables
  // -----------------------------------------
  
  // Construct the space vector:
  for (unsigned int n = 0; n < mNumberOfSpaces; n++) {
    boost::multi_array<Cell*, 3> space(boost::extents[mvSizeX[n]][mvSizeY[n]][mvSizeZ[n]]);
    
    for(size_t i = 0; i < space.shape()[0]; i++) {
      for(size_t j = 0; j < space.shape()[1]; j++) {
        for(size_t k = 0; k < space.shape()[2]; k++) {
          space[i][j][k] = 0; // Make the space empty:
        }
      }
    }
    mvSpace.push_back(space);
  }

  // default limiting edges:
  for (unsigned int n = 0; n < mNumberOfSpaces; n++) {
    
    unsigned int max_size = 0;

    for (int i = 0; i< 3; i++) {
      if (mvSpace[n].shape()[i] > max_size) 
        max_size = mvSpace[n].shape()[i];
    }

    std::array<bool, 3> hard_limit;
    for (int i = 0; i < 3; i++) {
      hard_limit[i] = mvSpace[n].shape()[i] == max_size ? 1 : 0;
    }
    mIsLimit.push_back(hard_limit);
  }


  // Create all cell types
  CellType* pType = 0;
    
  for (unsigned int i = 0; i < mNumberOfClones; i++) {
    
    pType = new CellType(i+1, // New cell type:
                         mvBirthrates[i],
                         mvDeathrates[i],
                         mvAggressions[i],
                         mvPushPower[i],
                         mvMutationrates[i]);
    
    mpTypes.push_back(pType);
  }
  
  // Insert cell of the first type into the center
  Cell* pCell = new Cell(mpTypes[0]);
  
  InsertCell(mvUniverses[0],
             CENTER(mvSizeX[mvUniverses[0]]),
             CENTER(mvSizeY[mvUniverses[0]]),
             CENTER(mvSizeZ[mvUniverses[0]]),
             pCell);
  
  pCell->MutateCellFixedNumber(mNumberClonalMutations); // clonal variants
  
  
  // record inital transformation event:
  mvTransformX.push_back(CENTER(mvSizeX[mvUniverses[0]]));
  mvTransformY.push_back(CENTER(mvSizeY[mvUniverses[0]]));
  mvTransformZ.push_back(CENTER(mvSizeZ[mvUniverses[0]]));
  mvTransformTime.push_back(TICK_SOURCE_STARTTIME);
  mvTransformTypeFrom.push_back(-1);
  mvTransformTypeTo.push_back(0);
  
  
  // fill event list
  for (unsigned int i = 1; i < mNumberOfClones; i++) {
    events.add(mvCloneStartTimes[i], EVENT_SUBCLONE, i);
  }
  
  for (unsigned int i = 0; i < mNumberOfClones; i++) {
    if (mvKillRegrowTime[i] > 0) {
      events.add(mvKillRegrowTime[i], EVENT_KILLREGROW, i);
    }
  }
  
}

// Destructor:
Universe::~Universe(){
  D(Rcpp::Rcout << "Universe::~Universe() " << this << std::endl;)

  // Remove all Cells:
  for (size_t n = 0; n < mvSpace.size(); n++) {
    for (size_t x = 0; x < mvSpace[n].shape()[0]; x++) {
      for (size_t y = 0; y < mvSpace[n].shape()[1]; y++) {
        for (size_t z = 0; z < mvSpace[n].shape()[2]; z++) {
          if (ContainsCell(n, x, y, z)) {
            delete mvSpace[n][x][y][z];
          }
        }
      }
    }
  }

  // Remove all CellTypes:
  for(std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++) {
    delete mpTypes[i];
  }

  // Remove the complete phylogeny:
  for(std::vector<PhylogenyNode*>::size_type i = 0;
     i < mpPhylogenies.size(); i++) {
    delete mpPhylogenies[i];
  }
}


// Getter functions:
long double Universe::Time() const {return mTime;};

Rcpp::NumericVector Universe::GetTime() const { 
  return Rcpp::NumericVector::create(TICK_SOURCE_STARTTIME);
}

inline bool Universe::LimitIsReached() const {return mLimitReached;};

inline bool  Universe::ContainsCell(int i, int x, int y, int z) const {
  // Returns true if the position i, x, y, z contains a cell.
  if (!ContainsCoordinate(i,x,y,z)) return 0; 
  return mvSpace[i][x][y][z] != 0;
}

bool Universe::ContainsCoordinate(int i, int x, int y, int z) const {
  // Returns true if the position x, y, z is within the limits of the universe.

  return x < mvSpace[i].shape()[0] &&
         y < mvSpace[i].shape()[1] &&
         z < mvSpace[i].shape()[2] &&
         x >= 0 &&
         y >= 0 &&
         z >= 0;
}

bool Universe::TouchesHardLimit(int n, int x, int y, int z) const {
  // Returns true if the position is touching a edge that leads to termination.
  
  int coord[3] = {x,y,z};

  for (int d = 0; d < 3; d++) { // for x, y, z
    if (mIsLimit[n][d]) { // limiting dimension
      if (coord[d] <= 0 || coord[d] >= (mvSpace[n].shape()[d] - 1)) {
        return true;
      }
    }
  }
  
  return false;
}

/*
void Universe::SetHardLimits(bool lx, bool, ly, bool lz) const {
  // Returns true if the position x, y, z is within the limits of the universe.
  mIsLimit = {lx, ly, lz};
}
*/


void  Universe::Size(int n, int& x, int& y, int& z) const {
  x = mvSpace[n].shape()[0];
  y = mvSpace[n].shape()[1];
  z = mvSpace[n].shape()[2];
}

Cell *Universe::GetCell(int i, int x, int y, int z) {
  return mvSpace[i][x][y][z];
}

Cell* Universe::RemoveCell (Cell* pCell) {
  return RemoveCell(pCell->N(), pCell->X(), pCell->Y(), pCell->Z());
}

Cell* Universe::RemoveCell (int i, int x, int y, int z) {
  // Returns a cell after removal from the space of the universe.

  Cell* p_removed_cell = mvSpace[i][x][y][z];
  mvSpace[i][x][y][z] = 0;
  p_removed_cell->Location(0, 0, 0, 0);
  return p_removed_cell;
}

CellType* Universe::NextReactionType(long double *r_delta_time, int *action){
  // Samples and returns the next reaction type that occures.

  //// Return null if the universe contains no types:
  //if (mpTypes.size() == 0) {
  //  Rcpp::Rcerr << "Error in Universe::NextReactionType(double *, int *):\n";
  //  Rcpp::Rcerr << "   > No types to choose from!" << std::endl;
  //  Rcpp::stop("No type to pick for next reaction.");
  //}
   
  // Debug messages:
  D(Rcpp::Rcout << std::endl;)
  D(Rcpp::Rcout << "########## Sampling ###########" << std::endl;)
  D(Rcpp::Rcout << "  Next reaction type:" << std::endl;)
  D(Rcpp::Rcout << std::endl;)

  // Store reaction specs in these vectors:
  std::vector<long double> deltas;
  std::vector<int> reaction_types;
  std::vector<std::vector<CellType*>::size_type> cell_types;
  D(int n_react = 0;)

  boost::uniform_real<double> runi_dist(0.0, 1.0);

  // Calculate dt for all types in the universe:
  for(std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++) {
    
    unsigned long s_number = mpTypes[i]->NumMembers();
    if (s_number == 0) continue; // skip if number of members is zero
    
    for(int reaction_type = 1; reaction_type <= 2; reaction_type++) {

      double runi = runi_dist(rng);

      double reaction_rate;
      switch(reaction_type){
        case 1:
          reaction_rate = mpTypes[i]->Deathrate();
          break;
        case 2:
          reaction_rate = mpTypes[i]->Birthrate();
          break;
        default: //undefined action requested
          Rcpp::Rcerr << "Error in Universe::NextReactionType:\n";
          Rcpp::Rcerr << "   > Action " <<  reaction_type << " undefined.\n";
          Rcpp::Rcerr << std::endl;
          Rcpp::stop("Bad reaction.");
      }

      long double delta_t = (log(1) - log(runi)) / (s_number * reaction_rate);

      // Debug messages:
      D(Rcpp::Rcout << "    Candidate reaction: " << ++n_react << std::endl;)
      D(Rcpp::Rcout << "      Species: " << i << std::endl;)
      D(Rcpp::Rcout << "      Species ID: " << mpTypes[i]->Id() << std::endl;)
      D(Rcpp::Rcout << "      Number of members : " << s_number << std::endl;)
      D(Rcpp::Rcout << "      Reaction type: " << reaction_type << std::endl;)
      D(Rcpp::Rcout << "      Reaction rate: " << reaction_rate << std::endl;)
      D(Rcpp::Rcout << "      Random uint: " << runi << std::endl;)
      D(Rcpp::Rcout << "      Delta t: " << delta_t << std::endl;)
      D(Rcpp::Rcout << std::endl;)

      // Collect candidate reactions:
      deltas.push_back(delta_t);
      reaction_types.push_back(reaction_type);
      cell_types.push_back(i);
    }
  }

  // Find minimum delta t of all reactions:
  int i_min;
  long double delta_time_minimum;
  for(std::vector<double>::size_type i = 0; i < deltas.size(); i++){
    if(i == 0 || delta_time_minimum > deltas[i]){
      i_min = i;
      delta_time_minimum = deltas[i];
    }
  }

  // Debug messages:
  D(Rcpp::Rcout << "  Next reaction:" << std::endl;)
  D(Rcpp::Rcout << "    Reaction id: " << reaction_types[i_min] << std::endl;)
  D(Rcpp::Rcout << "    Species id: " << cell_types[i_min] << std::endl;)
  D(Rcpp::Rcout << "    dt: " << deltas[i_min] << std::endl;)
  D(Rcpp::Rcout << "###############################" << std::endl;)

  // Return results:
  *r_delta_time = deltas[i_min];
  *action = reaction_types[i_min];
  
  // set reaction type to 3 if nearest neighboor should be found.
  if (reaction_types[i_min] == 2 && (mExploreNeighborhood || GROW_TO_CLOSEST_EDGE)) {
    *action = 3;
  }
    
  return(mpTypes[cell_types[i_min]]);
}

std::vector< std::array<int, 3> > Universe::FreeNeighbours(int i, int x, int y, int z) const {
  std::vector< std::array<int, 3> > result;
  std::array<int, 3> cords;

  int xstart, xend;
  if (mvSpace[i].shape()[0] > 1) {
    xstart = x - 1;
    xend = x + 1;
  } else {
    xstart = x;
    xend = x;
  }

  int ystart, yend;
  if (mvSpace[i].shape()[1] > 1) {
    ystart = y - 1;
    yend = y + 1;
  } else {
    ystart = y;
    yend = y;
  }

  int zstart, zend;
  if (mvSpace[i].shape()[2] > 1) {
    zstart = z - 1;
    zend = z + 1;
  } else {
    zstart = z;
    zend = z;
  }

  // Get all free neighbours:
  for (int xn = xstart; xn <= xend; xn++) {
    for (int yn = ystart; yn <= yend; yn++) {
      for (int zn = zstart; zn <= zend; zn++) {
        if (!ContainsCell(i, xn, yn, zn) && 
            ContainsCoordinate(i, xn, yn, zn) && 
            (xn != x || yn != y || zn != z))
        {
          cords[0] = xn; cords[1] = yn; cords[2] = zn;
          result.push_back(cords);
        }
      }
    }
  }
  return result;
}


void Universe::FindClosestFreeGrid(int i, int &x, int &y, int &z, double max_dist) const {

  class SpiralOut{
  protected:
    unsigned layer;
    unsigned leg;
  public:
    int x, y; //read these as output from next, do not modify.
    SpiralOut():layer(1),leg(0),x(0),y(0){}
    void goNext(){
      switch(leg){
      case 0: ++x; if(x  == layer)  ++leg;                break;
      case 1: ++y; if(y  == layer)  ++leg;                break;
      case 2: --x; if(-x == layer)  ++leg;                break;
      case 3: --y; if(-y == layer){ leg = 0; ++layer; }   break;
      }
    }
  };
  
  SpiralOut spiral;
  int x_t = 0;
  int y_t = 0;
  int d = 0;
  
  do {
    x_t = x + spiral.x;
    y_t = y + spiral.y;
    d = abs(spiral.x) + abs(spiral.y); // manhattan distance
    spiral.goNext();
  } while ((!ContainsCoordinate(i,x_t,y_t,z) || ContainsCell(i,x_t,y_t,z)) && d < max_dist);

  if (!ContainsCell(i,x_t,y_t,z)) {
    x = x_t;
    y = y_t;
  }
}



// get total cell counts per cell type
std::vector <unsigned long> Universe::CellCountsPerType() const {
  std::vector <unsigned long> result;
  
  for (size_t i = 0; i < mpTypes.size(); i++) {
    result.push_back(mpTypes[i]->NumMembers());
  }
  
  return result;
}

unsigned long Universe::CellCounts() const { 
  unsigned long cell_counts = 0;
  
  for (size_t i = 0; i < mpTypes.size(); i++) {
    cell_counts += mpTypes[i]->NumMembers();
  }

  return cell_counts;
}

Rcpp::List Universe::GetSnapshots() const {
  return mUniverseSnapshots;
}

bool Universe::GetExploreLocally() const {return mExploreNeighborhood; };

bool Universe::GetRecordHistory() const {return mFlagRecordHistory; };

bool Universe::GetReturnGenerationTree() const {return mFlagReturnGenerationTree; };

unsigned int Universe::GetClonalMutations() const { return mNumberClonalMutations; }

unsigned int Universe::GetStopPopulationSize() const { return mStopPopulationSize; }
bool Universe::GetAltPush() const{ return mAltPush; }


// time
Rcpp::NumericVector Universe::GetSnapshotTimes() const { 
 
  Rcpp::NumericVector stime = Rcpp::NumericVector(mvSnapshotTimes.size());
  
  for (std::size_t i = 0; i < mvSnapshotTimes.size(); i++) {
    stime[i] = mvSnapshotTimes[i];
  }
  
  return stime; 
}


// Sampling functions (getters):
CellSample Universe::TakeSample(int n, Shape* pShape){

  // Counter for the number of sampled nodes:
  unsigned long int sampled_cells_cnt = 0;

  // Ints to store the next cell to pick:
  int x, y, z;
  
  // Stage all samples for sequencing:
  while (pShape->next_coordinate(x, y, z)) { // Get next location from shape.
    if (ContainsCell(n,x,y,z)) { //selected location is not empty
      
      std::stringstream label;
      label << std::hex << mvSpace[n][x][y][z]->Id()
            << std::dec<<" ("<<x<<","<<y<<","<<z<<")";
      
      sampled_cells_cnt++; // Increment counters of cells:
      mvSpace[n][x][y][z]->AssociatedNode()->StageNodeForSequencing(label.str());
    }
  }

  // The resulting sample will be stored here:
  CellSample sample;
  sample.SetSampledCellNumber(sampled_cells_cnt);

  // Now sequence the staged cells:
  for (size_t i = 0; i < mpPhylogenies.size(); i++)
  {
    mpPhylogenies[i]->CollectStagedNodes(sample);
    mpPhylogenies[i]->UnstageNodes();
  }

  return sample;
}

std::vector <CellSample> Universe::TakeMultipleSamples(int sample_strategy) {

  switch(sample_strategy){
    case 1: // Take squares low depth (20X)
      // ...
      break;
    case 2: // Take squares high depth (100X)
      // ...
      //TakeSquares(radius, sim_id, size_x, size_y, size_z, depth, min_reads,
      //  pois_depth, sampgrid_output_file, bulk_output_file);
      break;
    case 20:
      // to do: TakeStripes
    break;
  }

  std::vector <CellSample> pEmptyVector;
  return pEmptyVector;
}


std::vector <std::string> Universe::SingleCellTrees(int n, Shape *pShape) {
  return SingleCellTrees(n, pShape, true);
}


std::vector <std::string> Universe::SingleCellTrees(int n, Shape *pShape, bool add_label) {
  
  std::vector <std::string> result;
  std::stringstream outstream;
  
  // Ints to store the next cell to pick:
  int x, y, z;
  
  // Stage all samples for sequencing:
  while (pShape->next_coordinate(x, y, z)) { // Get next location from shape.
    if (ContainsCell(n,x,y,z)) { //selected location is not empty
      
      Cell *mpCell = mvSpace[n][x][y][z];
      
      if (add_label) {
        
        std::stringstream label;
        label << "cell" << mpCell->Id() << "_"
              << "type" << mpCell->Type()->Id() << "_"
              << "x" << x<< "_" << "y" << y << "_" << "z" << z;
        mvSpace[n][x][y][z]->AssociatedNode()->StageNodeForSequencing(label.str());
        
      } else {
        mvSpace[n][x][y][z]->AssociatedNode()->StageNodeForSequencing();
      }

    }
  }
  
  // Now sequence the staged cells:
  for (size_t i = 0; i < mpPhylogenies.size(); i++)
  {
    outstream.str("");
    outstream << "tree" << i;
    mpPhylogenies[i]->StagedNodesToStream(outstream);
    outstream << ";";
    mpPhylogenies[i]->UnstageNodes();
    result.push_back(outstream.str());
  }
  
  return result;
}


// Setter functions:
void Universe::IncrementTimeBy(long double Delta){ mTime += Delta; }

void Universe::MarkLimitIsReached() { mLimitReached = true; }

bool Universe::InsertCell(int n, int x, int y, int z, Cell* pCell){

  // Debug messages:
  D(Rcpp::Rcout << std::endl;)
  D(Rcpp::Rcout << "########## Insertion ##########" << std::endl;)
  D(Rcpp::Rcout << "   ID: " << pCell->Id() << std::endl;)
  D(Rcpp::Rcout << "   Location:" << std::endl;)
  D(Rcpp::Rcout << "       Universe: " << this << std::endl;)
  D(Rcpp::Rcout << "       x: " << x << std::endl;)
  D(Rcpp::Rcout << "       y: " << y << std::endl;)
  D(Rcpp::Rcout << "       z: " << z << std::endl;)
  D(Rcpp::Rcout << "###############################" << std::endl;)

  
  // Dont allow inserts outside of the array extends:
  if (ContainsCoordinate(n,x,y,z) == false) {
    Rcpp::Rcout << "Error: Insertion of cell failed:" << std::endl;
    Rcpp::Rcout << "Error: Coordinates (i+1) were: " << n << x + 1 << ",";
    Rcpp::Rcout << y + 1 << "," << z + 1 << std::endl;
    Rcpp::Rcout << "Error:  Universe size was: "
              << mvSpace[n].shape()[0] << ","
              << mvSpace[n].shape()[1] << ","
              << mvSpace[n].shape()[2] << std::endl;
    Rcpp::stop("Cell outside of valid positions.");

  }

  // Dont allow inserts on top of cell:
  if (ContainsCell(n,x,y,z)) {
    Rcpp::Rcout << "Error: Insertion of cell failed:" << std::endl;
    Rcpp::Rcout << "Error: Coordinates (i) were: " << x << ",";
    Rcpp::Rcout << y << "," << z  << std::endl;
    Rcpp::Rcout << "Error:  Position was taken by" << mvSpace[n][x][y][z] << std::endl;
    Rcpp::stop("Cell inserted in location of other cell.");
  }

  // Update position and universe of the cell:
  pCell->Location(this, n, x, y, z);
  mvSpace[n][x][y][z] = pCell;

  // Register new phylogeny:
  if (pCell->AssociatedNode() == 0) {
    PhylogenyNode* pNewPhylo = new PhylogenyNode(pCell);
    mpPhylogenies.push_back(pNewPhylo);
  }

  return true;
}


void Universe::RecordToHistory() {
  for (size_t i = 0; i < mpTypes.size(); i++) {
    mvCellCounts[i].push_back(mpTypes[i]->NumMembers());
  }
  mvHistoryTime.push_back(mTime);
}


void Universe::SetSnapshotTimes(Rcpp::NumericVector new_times) {
  
  if (mvSnapshotTimes.size()) {
    Rcpp::stop("Already set snapshot times can't be changed.");
  }
  
  // sort and remove duplicates:
  Rcpp::NumericVector tmp_copy = clone(new_times);
  std::sort(tmp_copy.begin(), tmp_copy.end());
  for (std::size_t i = new_times.size() - 1; i > 0; i--) {
    if (tmp_copy[i] == tmp_copy[i-1]) tmp_copy.erase(i);
  }
  
  // store snapshot times in event list:
  for (std::size_t i = 0; i < tmp_copy.size(); i++) {
    events.add(tmp_copy[i], EVENT_SNAPSHOT, i);
  }
  
};

void Universe::SetExploreLocally(bool new_value) {
  mExploreNeighborhood = new_value;  
};

void Universe::SetRecordHistory(bool new_value) {
  mFlagRecordHistory = new_value;  
};

void Universe::SetReturnGenerationTree(bool new_value) {
  mFlagReturnGenerationTree = new_value;  
};

void Universe::SetClonalMutations(unsigned int n) {
  mpPhylogenies[0]->SetMutations(n);
  mNumberClonalMutations = n;
}

void Universe::SetStopPopulationSize(unsigned int n) {
  mStopPopulationSize = n;
}

void Universe::SetAltPush(bool value) { 
  mAltPush = value; 
}




void Universe::ScaleMutationTree(double scale_factor, Rcpp::List fixed_vals, double dispersion) {
  
  // Scale all trees:
  for (size_t i = 0; i < mpPhylogenies.size(); i++){
    mpPhylogenies[i]->ScaleGenerationTree(scale_factor, fixed_vals, dispersion);
  }
}


// Output Functions:

void Universe::PrintSimulationParameters() const {
  
  // Print arguments:
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "########## Options #############" << std::endl;
  
  for (unsigned int i = 0; i < mNumberOfSpaces; i++) {
    Rcpp::Rcout << "  Size " << i << ": "
                << mvSizeX[i] << "x" << mvSizeY[i] << "x" << mvSizeZ[i] << "\n";
  }
  
  // clone level data
  for (unsigned int i = 0; i < mNumberOfClones; i++) { // for each clone seperatly:
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "  Mutation rate: " << mvMutationrates[i] << std::endl;
    Rcpp::Rcout << "  Birth rate: " << mvBirthrates[i] << std::endl;
    Rcpp::Rcout << "  Death rate: " << mvDeathrates[i] << std::endl;
    Rcpp::Rcout << "  Aggression: " << mvAggressions[i] << std::endl;
    Rcpp::Rcout << "  Push power: " << mvPushPower[i] << std::endl;
    Rcpp::Rcout << "  Clone start time: " << mvCloneStartTimes[i] << std::endl;
    Rcpp::Rcout << "  Kill & regrow time: " << mvKillRegrowTime[i] << std::endl;
    Rcpp::Rcout << "  Father: " << mvFathers[i] << std::endl;
    Rcpp::Rcout << "  Universe: " << mvUniverses[i] << std::endl;
  } // end printing of clone level data
  
  Rcpp::Rcout << "  Random seed: " << mSeed << std::endl;
  Rcpp::Rcout << "################################" << std::endl;
  Rcpp::Rcout << std::endl;
}

void  Universe::RecordSnapshot() {
  unsigned int n = 0;
  
  //Rcpp::Rcout << "creating snapeshot at time " << TICK_SOURCE_STARTTIME << std::endl;
  
  // collect types of all non empty locations
  std::vector<int> vx,vy,vz,vid, vtype, vmburden;
  unsigned int type = 0;
  
  for(size_t i = 0; i < mvSpace[n].shape()[0]; i++) {
    for(size_t j = 0; j < mvSpace[n].shape()[1]; j++) {
      for(size_t k = 0; k < mvSpace[n].shape()[2]; k++) {
        if ((type = GetCellType(i,j,k))) {
          
          // cell type data
          vx.push_back(i);
          vy.push_back(j);
          vz.push_back(k);
          vid.push_back(mvSpace[n][i][j][k]->Id());
          vtype.push_back(type);
          vmburden.push_back(mvSpace[n][i][j][k]->TotalMutNum());
        }
      }
    }
  }
  
  // construct data frame and store in snapshot list
  Rcpp::DataFrame df = 
    Rcpp::DataFrame::create(
      Rcpp::Named("x") = vx,
      Rcpp::Named("y") = vy,
      Rcpp::Named("z") = vz,
      Rcpp::Named("id") = vid,
      Rcpp::Named("type") = vtype,
      Rcpp::Named("burden") = vmburden);
  
  mUniverseSnapshots.push_back(df);
  mvSnapshotTimes.push_back(mTime);
    
}
 

// Other functions:
bool Universe::PushCell(int n, double newX, double newY, double newZ,
                        double dx, double dy, double dz,
                        unsigned int push) {


  // Determine integer coordinates of current cell:
  int intCurX = round(newX);
  int intCurY = round(newY);
  int intCurZ = round(newZ);
  int intNewX;
  int intNewY;
  int intNewZ;

  // Increase till new integer coords jumps:
  do {
    newX += dx;
    newY += dy;
    newZ += dz;
    intNewX = round(newX);
    intNewY = round(newY);
    intNewZ = round(newZ);
  } while (intNewX == intCurX && intNewY == intCurY && intNewZ == intCurZ);


  // If you hit a border stop pushes.
  if (ContainsCoordinate(n, intNewX, intNewY, intNewZ) == false) {
    D(Rcpp::Rcout << "  Hit border during push." << std::endl;)
    //this->MarkLimitIsReached();
    return false;
  }

  // Test if next space is free.
  if (ContainsCell(n, intNewX, intNewY, intNewZ) == false ||
      (push > 0 && PushCell(n, newX, newY, newZ, dx, dy, dz, push - 1)))
  {
    // Debug messages:
    D(Rcpp::Rcout << "  Push sucessfull:" << std::endl;)
    D(Rcpp::Rcout << "    id: " << mvSpace[n][intCurX][intCurY][intCurZ] << std::endl;)
    D(Rcpp::Rcout << "    x: " << intCurX << "->" << intNewX << std::endl;)
    D(Rcpp::Rcout << "    y: " << intCurY << "->" << intNewY << std::endl;)
    D(Rcpp::Rcout << "    z: " << intCurZ << "->" << intNewZ << std::endl;)
    D(Rcpp::Rcout << "    push: " << push << "->" << intNewZ << std::endl;)
    D(Rcpp::Rcout << std::endl;)

    // Move cell:
    mvSpace[n][intNewX][intNewY][intNewZ] = mvSpace[n][intCurX][intCurY][intCurZ];
    mvSpace[n][intCurX][intCurY][intCurZ]->Location(intNewX, intNewY, intNewZ);
    mvSpace[n][intCurX][intCurY][intCurZ] = 0;
    return true;
  }

  return false;
}


bool Universe::PushCell(int N, int x1, int y1, int z1, int x2, int y2, int z2) {
  
  D(std::cout << "Moving cells from (" << x1 << "," << y1 << "," << z1 << ") to (" << x2 << "," << y2 << "," << z2 << ")" << std::endl;)
  
  /*
  // Check that target coordinates are valid. 
  if (!ContainsCoordinate(N, x2, y2, z2)) {
    Rcpp::stop("Target coordinates are outside of simulation space!");
  }
  
  // Check that target coordinates are empty.
  if (ContainsCell(N, x2, y2, z2)) {
    Rcpp::stop("Target coordinates are not empty!");
  }
  */
  
  ////////////////////////////////////////////////////////////////////////
  // Code from Bresenham3D (https://gist.github.com/yamamushi/5823518)  //
  ////////////////////////////////////////////////////////////////////////
  
  // shared variables
  int i, err_1, err_2;
  Cell *last_cell = 0, *this_cell;
  
  // current point
  int point[3] = {x1, y1, z1};;

  // differences
  int dx = x2 - x1;
  int dy = y2 - y1;
  int dz = z2 - z1;

  // increments
  int x_inc = (dx < 0) ? -1 : 1;
  int y_inc = (dy < 0) ? -1 : 1;
  int z_inc = (dz < 0) ? -1 : 1;
    
  // total length
  int l = abs(dx);
  int m = abs(dy);
  int n = abs(dz);
  
  // dn2
  int dx2 = l << 1;
  int dy2 = m << 1;
  int dz2 = n << 1;

  
  if ((l >= m) && (l >= n)) {
    
    err_1 = dy2 - l;
    err_2 = dz2 - l;
    
    for (i = 0; i < l; i++) {
    
    // Move cell:
    if (last_cell) {
      last_cell->Location(point[0], point[1], point[2]);
    } else {
      if (i > 0) {
        return true;
      }
    }
    this_cell = mvSpace[N][point[0]][point[1]][point[2]];
    mvSpace[N][point[0]][point[1]][point[2]] = last_cell;
    last_cell = this_cell;
    D(std::cout << " - Moving cell to (" << point[0] << "," << point[1] << "," << point[2] << ")" << std::endl;)
      
      
      
      if (err_1 > 0) {
        point[1] += y_inc;
        err_1 -= dx2;
      }
      
      if (err_2 > 0) {
        point[2] += z_inc;
        err_2 -= dx2;
      }
      
      err_1 += dy2;
      err_2 += dz2;
      point[0] += x_inc;
    }
    
  } else if ((m >= l) && (m >= n)) {
    
    err_1 = dx2 - m;
    err_2 = dz2 - m;
    
    for (i = 0; i < m; i++) {
      
      // Move cell:
      if (last_cell) {
        last_cell->Location(point[0], point[1], point[2]);
      } else {
        if (i > 0) {
          return true;
        }
      }
      this_cell = mvSpace[N][point[0]][point[1]][point[2]];
      mvSpace[N][point[0]][point[1]][point[2]] = last_cell;
      last_cell = this_cell;
      D(std::cout << " - Moving cell to (" << point[0] << "," << point[1] << "," << point[2] << ")" << std::endl;)
      
      
      if (err_1 > 0) {
        point[0] += x_inc;
        err_1 -= dy2;
      }
      
      if (err_2 > 0) {
        point[2] += z_inc;
        err_2 -= dy2;
      }
      
      err_1 += dx2;
      err_2 += dz2;
      point[1] += y_inc;
    }
    
  } else {
    
    err_1 = dy2 - n;
    err_2 = dx2 - n;
    
    for (i = 0; i < n; i++) {
      
      // Move cell:
      if (last_cell) {
        last_cell->Location(point[0], point[1], point[2]);
      } else {
        if (i > 0) {
          return true;
        }
      }
      this_cell = mvSpace[N][point[0]][point[1]][point[2]];
      mvSpace[N][point[0]][point[1]][point[2]] = last_cell;
      last_cell = this_cell;
      D(std::cout << " - Moving cell to (" << point[0] << "," << point[1] << "," << point[2] << ")" << std::endl;)
        
      if (err_1 > 0) {
        point[1] += y_inc;
        err_1 -= dz2;
      }
        
      if (err_2 > 0) {
        point[0] += x_inc;
        err_2 -= dz2;
      }
        
      err_1 += dy2;
      err_2 += dx2;
      point[2] += z_inc;
    }
  }
  
  // put last cell into target location  
  D(std::cout << " - Moving cell to (" << point[0] << "," << point[1] << "," << point[2] << ")" << std::endl;)
  if (last_cell) last_cell->Location(point[0], point[1], point[2]);
  mvSpace[N][point[0]][point[1]][point[2]] = last_cell;
  
  return true;
}




// Simulate tumours:
bool Universe::RunSimulation(bool verbose) {
  
  if (verbose) { PrintSimulationParameters(); } // Print arguments:
  //if (verbose) { events.print(); };
  
  rng.seed(mSeed); // Set random seed
  RecordToHistory();
  
  // -----------------------------------------
  // Variables
  // -----------------------------------------
  long double dt = 0.0L;
  int action;
  Cell* pCell = 0;
  unsigned int max_cells = mStopPopulationSize;
  // -----------------------------------------

  
  // -----------------------------------------
  // Simulation code
  // -----------------------------------------

  unsigned int n_steps = 0; // used to check for interrupts
  
  try { // block that can be interrupted by user
    
    event next_event = events.next(); 
    while (next_event.event_type != EVENT_NON && !LimitIsReached()) {
      
      if (n_steps++ % 100000 == 0) Rcpp::checkUserInterrupt(); // check for user interrupt
      
      while (next_event.event_type != EVENT_NON && next_event.time <= TICK_SOURCE_STARTTIME) {
        
        switch(next_event.event_type) {
        
        // Event subclone ******************************************************
        case EVENT_SUBCLONE:
          
        {
          
          // delay events for subclones with only one member:
          if (mpTypes[mvFathers[next_event.event_idx]]->NumMembers() < 2) { 
          
          D(Rcpp::Rcout << "Only one member for cell type. Delaying event." << std::endl;)
          
          events.add(
            next_event.time + TICK_MIN_INCREMENT, 
            next_event.event_type, 
            next_event.event_idx
          );
          
          D(events.print();)
            
        } else {
          
          // Debug messages:
          if (verbose) {
            Rcpp::Rcout << std::endl;
            Rcpp::Rcout << "########## Comment #############" << std::endl;
            Rcpp::Rcout << " New clone: " << next_event.event_idx << std::endl;
            Rcpp::Rcout << " Time: " << TICK_SOURCE_STARTTIME << std::endl;
            Rcpp::Rcout << "################################" << std::endl;
            Rcpp::Rcout << std::endl;
          }
          
          // introduce a new cell or change a existing one
          if (mvFathers[next_event.event_idx] == -1) { // new cell of a given type.
            pCell = new Cell(mpTypes[next_event.event_idx]);
            Rcpp::stop("Feature does not exist. Insert new cell.");
            // Need cells locations for this.
          } else { // existing cell that changes it's type.
            pCell = mpTypes[mvFathers[next_event.event_idx]]->RandomMember(); // random cell to transform
            pCell->Type(mpTypes[next_event.event_idx]);
          }
          
          // potentially move the cell into a new space.
          int new_space = mvUniverses[next_event.event_idx];
          int current_space = pCell->N();
          
          if (current_space != new_space) {
            
            pCell = RemoveCell(pCell);
            
            InsertCell(new_space,
                       CENTER(mvSizeX[new_space]),
                       CENTER(mvSizeY[new_space]),
                       CENTER(mvSizeZ[new_space]),
                       pCell);
          }
          
          
          // record the transformation event:
          mvTransformX.push_back(pCell->X());
          mvTransformY.push_back(pCell->Y());
          mvTransformZ.push_back(pCell->Z());
          mvTransformTime.push_back(TICK_SOURCE_STARTTIME);
          mvTransformTypeFrom.push_back(mvFathers[next_event.event_idx]);
          mvTransformTypeTo.push_back(next_event.event_idx);
        }
        }
          
          break;
          
          // Event kill regrowth *************************************************
        case EVENT_KILLREGROW:
          
        {
          int curr_popsize = 0;
          for (size_t i = 0; i < mpTypes.size(); i++) {
            curr_popsize += mpTypes[i]->NumMembers();
          }
          
          int popsize_tokill = KILL_REGROW_PERCENT * curr_popsize;
          
          int popsize_killed = 0;
          while (popsize_killed < popsize_tokill) {
            
            // sample random type to kill
            boost::random::uniform_int_distribution<int> rand_celltype(0, mpTypes.size());
            int rand_celltype_i = rand_celltype(rng);
            
            if (mpTypes[rand_celltype_i]->NumMembers() > 1) {
              
              popsize_killed += 1;
              mpTypes[rand_celltype_i]->RandomMember()->Kill();
            }
          }
        }
          
          break;
          
        case EVENT_SNAPSHOT:
          
        {
          RecordSnapshot();
          mSnapshotIndex++;
        }
          
          break;
          
          
        case EVENT_NON:
          
          break;
          
        default:
          Rcpp::stop("Unkown event type %c", next_event.event_type);
        }
        
        next_event = events.next();
        
      }
      
      if (mFlagRecordHistory && n_steps % 10 == 0) RecordToHistory();
      NextReactionType(&dt, &action)->RandomMember()->DoAction(&action);
      IncrementTimeBy(dt);
      if (max_cells > 0 && CellCounts() >= max_cells) MarkLimitIsReached();
    }
     
  
    while(!LimitIsReached()) { // run till limit reached:
      if (n_steps++ % 100000 == 0) Rcpp::checkUserInterrupt(); // check for user interrupt
      if (mFlagRecordHistory && n_steps % 10 == 0) RecordToHistory();
      NextReactionType(&dt, &action)->RandomMember()->DoAction(&action);
      IncrementTimeBy(dt);
      if (max_cells > 0 && CellCounts() >= max_cells) MarkLimitIsReached();
    }
    
  } catch (Rcpp::internal::InterruptedException& e) { // catch interruptions
    return false; // maybe we could return the unfinished data too? 
  }
    
  return true;
}



Rcpp::List Universe::TakeSampleBox(
    int minx, int miny, int minz, int maxx, int maxy, int maxz,
    double dp, int dp_model, int min_reads, double min_vaf,
    int seed, double purity)
{
  
  if (seed != 0) {
    rng.seed(seed); // Set random seed, if its non 0.
  }
  
  // Sample
  Shape* sample_shape = new Box(minx, miny, minz, maxx, maxy, maxz);
  CellSample sample = TakeSample(0, sample_shape);
  SequencingResult result = sample.Sequence(dp, dp_model, min_reads, min_vaf, purity);
  delete sample_shape;
  
  // Get std::vectors containing result elements:
  std::vector <unsigned int> alt_cpp_type = result.AltVector();
  std::vector <unsigned int> clone_id_cpp_type = result.CloneIdVector();
  std::vector <unsigned int> depth_cpp_type = result.DepthsVector();
  std::vector <std::string> mutation_id_cpp_type = result.MutationIdVector();
  std::vector <double> true_ccf_cpp_type = result.CCFVector();
  
  // Allocate R vector types:
  Rcpp::NumericVector alt(alt_cpp_type.size());
  Rcpp::NumericVector clone_id(clone_id_cpp_type.size());
  Rcpp::NumericVector depth_res(depth_cpp_type.size());
  Rcpp::CharacterVector mutation_id(mutation_id_cpp_type.size());
  Rcpp::NumericVector true_ccf(true_ccf_cpp_type.size());
  
  // Fill R vectors with the data:
  for (size_t i = 0; i < alt_cpp_type.size(); i++) { // vectors should all have the same length ...
    alt[i] = alt_cpp_type[i];
    depth_res[i] = depth_cpp_type[i];
    clone_id[i] = clone_id_cpp_type[i];
    mutation_id[i] = mutation_id_cpp_type[i];
    true_ccf[i] = true_ccf_cpp_type[i];
  }
  
  std::vector <unsigned long> cellcounts;
  cellcounts = CellCountsPerType();
  
  Rcpp::NumericVector cell_counts(cellcounts.size());
  for (size_t i = 0; i < cellcounts.size(); i++){
    cell_counts[i] = cellcounts[i];
  }
  
  // Return results as list containing a data.frame:
  return Rcpp::List::create(
    Rcpp::Named("number_mutations") =
      alt_cpp_type.size(),
      Rcpp::Named("mutation_data") =
        Rcpp::DataFrame::create(
          Rcpp::Named("clone") = clone_id,
          Rcpp::Named("alt") = alt,
          Rcpp::Named("depth") = depth_res,
          Rcpp::Named("id") = mutation_id,
          Rcpp::Named("ccf") = true_ccf,
          Rcpp::_["stringsAsFactors"] = false),
          Rcpp::Named("cell_counts") =
            Rcpp::DataFrame::create(
              Rcpp::Named("cellcounts") = cell_counts
            )
  );
}


Rcpp::List Universe::SingleCellTree(
    std::vector <int> x, std::vector <int> y, std::vector <int> z)
{
  return SingleCellTreeInternal(x, y, z, false);
}

Rcpp::List Universe::SingleCellTreeWithNodeIds(
    std::vector <int> x, std::vector <int> y, std::vector <int> z)
{
  return SingleCellTreeInternal(x, y, z, true);
}


Rcpp::List Universe::SingleCellTreeInternal(
    std::vector <int> x, std::vector <int> y, std::vector <int> z, bool include_labels)
{
  
  // Sample
  Shape* sample_shape = new PositionVector(x, y, z);
  std::vector <std::string> result = SingleCellTrees(0, sample_shape, !include_labels);
  delete sample_shape;
  
  // Allocate R vector types:
  Rcpp::CharacterVector tree_strings(result.size());
  
  // Fill R vectors with the data:
  for (size_t i = 0; i < result.size(); i++) {
    tree_strings[i] = result[i];
  }
  
  std::vector <unsigned long> cellcounts;
  cellcounts = CellCountsPerType();
  
  Rcpp::NumericVector cell_counts(cellcounts.size());
  for (size_t i = 0; i < cellcounts.size(); i++){
    cell_counts[i] = cellcounts[i];
  }
  
  // Return results as list containing a data.frame:
  //return tree_strings;
  return Rcpp::List::create(
    Rcpp::Named("tree_string") = tree_strings,
    Rcpp::Named("cell_counts") =
      Rcpp::DataFrame::create(
        Rcpp::Named("cellcounts") = cell_counts
      )
  );
}


Rcpp::List Universe::SingleCellTreeN(
    int nsamples, int grid_x, int grid_y, int grid_z)
{
  std::vector <int> x;
  std::vector <int> y;
  std::vector <int> z;
  
  int cnt = 0;
  while (cnt < nsamples) {
    
    boost::random::uniform_int_distribution<int> unif_coord_x(0, grid_x);
    int coord_x = unif_coord_x(rng);
    boost::random::uniform_int_distribution<int> unif_coord_y(0, grid_y);
    int coord_y = unif_coord_y(rng);
    boost::random::uniform_int_distribution<int> unif_coord_z(0, grid_z);
    int coord_z = unif_coord_z(rng);
    
    if (ContainsCell(0, coord_x, coord_y, coord_z)) {
      
      x.push_back(coord_x);
      y.push_back(coord_y);
      z.push_back(coord_z);
      cnt++;
    }
  }
  
  // Sample
  Shape* sample_shape = new PositionVector(x, y, z);
  std::vector <std::string> result = SingleCellTrees(0,sample_shape);
  delete sample_shape;
  
  // Allocate R vector types:
  Rcpp::CharacterVector tree_strings(result.size());
  
  // Fill R vectors with the data:
  for (size_t i = 0; i < result.size(); i++) {
    tree_strings[i] = result[i];
  }
  
  std::vector <unsigned long> cellcounts;
  cellcounts = CellCountsPerType();
  
  Rcpp::NumericVector cell_counts(cellcounts.size());
  for (size_t i = 0; i < cellcounts.size(); i++){
    cell_counts[i] = cellcounts[i];
  }
  
  // Return results as list containing a data.frame:
  //return tree_strings;
  return Rcpp::List::create(
    Rcpp::Named("tree_string") = tree_strings,
    Rcpp::Named("cell_counts") =
      Rcpp::DataFrame::create(
        Rcpp::Named("cellcounts") = cell_counts
      )
  );
}

unsigned int Universe::GetCellType(int x, int y, int z){
  int n = 0; // space number
  if (!ContainsCell(n, x, y, z)) return 0;
  return GetCell(n, x, y, z)->Type()->Id();
  
}

unsigned int Universe::GetMutationBurden(int x, int y, int z){
  int n = 0; // space number
  if (!ContainsCell(n, x, y, z))  return 0;
  return GetCell(n, x, y, z)->AssociatedNode()->MutationNumberInAncestry();
}



Rcpp::StringVector Universe::GetMutations(int x, int y, int z){
  
  int n = 0; // space number
  Rcpp::StringVector result_vector(0); // empty result vector
  if (!ContainsCell(n, x, y, z)) return result_vector;

  // get mutations of cell
  std::vector<std::string> mids;  // mutation ids
  GetCell(n, x, y, z)->AssociatedNode()->MutationsInAncestry(mids);
  
  // move to rcpp struc
  for (size_t i = 0; i < mids.size(); i++ ){
    result_vector.push_back(mids[i]);
  }
  
  return(result_vector);
}




Rcpp::DataFrame Universe::GetHistory() const {
  
  // temporary list and name vector
  Rcpp::List result_list(1 + mvCellCounts.size());
  Rcpp::CharacterVector result_names(1 + mvCellCounts.size());
  
  // time
  Rcpp::NumericVector time = Rcpp::NumericVector(mvHistoryTime.size());
  for (std::size_t i = 0; i < mvHistoryTime.size(); i++) {
    time[i] = mvHistoryTime[i];
  }
  result_list[0] = clone(time);
  result_names[0] = "time";
  
  // types
  for (std::size_t i = 0; i < mvCellCounts.size(); i++) { // vector over types
    
    Rcpp::IntegerVector type = Rcpp::IntegerVector(mvHistoryTime.size());
    for (std::size_t j = 0; j < mvCellCounts[i].size(); j++) {
      type[j] = mvCellCounts[i][j];
    }
    
    result_list[i+1] = clone(type);
    result_names[i+1] = "type"+std::to_string(i+1);
  }
  
  
  // create data frame
  Rcpp::DataFrame result(result_list);
  result.attr("names") = result_names;
  return result;
}


Rcpp::DataFrame Universe::GetTransformations() const {
  
  Rcpp::DataFrame df = 
    Rcpp::DataFrame::create(
      Rcpp::Named("x") = mvTransformX,
      Rcpp::Named("y") = mvTransformY,
      Rcpp::Named("z") = mvTransformZ,
      Rcpp::Named("time") = mvTransformTime,
      Rcpp::Named("from") = mvTransformTypeFrom,
      Rcpp::Named("to") = mvTransformTypeTo);
  
  return df;
}



Rcpp::IntegerVector Universe::FindClosestEdgeRcpp(int n, int x, int y, int z) const {
  int n_x = x, n_y = y, n_z = z;
  double phi = 0.0;
  FindClosestEdge(n, n_x, n_y, n_z, phi, 1e6);
  Rcpp::IntegerVector new_pos = {n_x, n_y, n_z};//(3);
  //new_pos[0] = n_x; new_pos[1] = n_y; new_pos[2] = n_z;
  return(new_pos);
}


double Universe::GetDistanceToEdge(int n, int &x, int &y, int &z, double phi) const {
  return GetDistanceToEdge(n, x, y, z, phi, sqrt(mvSizeX[n] * 1.0 * mvSizeY[n]), false);
}

double Universe::GetDistanceToEdgeRcpp(int n, int &x, int &y, int &z, double phi) const {
  return GetDistanceToEdge(n, x, y, z, phi, sqrt(mvSizeX[n] * 1.0 * mvSizeY[n]), false);
}

double Universe::GetDistanceToEdgeRcppV2(int n, int &x, int &y, int &z, double phi) const {
  return GetDistanceToEdge(n, x, y, z, phi, sqrt(mvSizeX[n] * 1.0 * mvSizeY[n]), true);
}

double Universe::GetDistanceToEdge(int n, int &x, int &y, int &z, double phi, double max_dist) const {
  if (max_dist < 0.0) {
    return GetDistanceToEdge(n, x, y, z, phi, sqrt(mvSizeX[n] * 1.0 * mvSizeY[n]), false);
  } else{
    return GetDistanceToEdge(n, x, y, z, phi, max_dist, false);
  }
}
  
double Universe::GetDistanceToEdge(int n, int &x, int &y, int &z, double phi, double max_dist, bool check_not_empty_loc) const {

  D(std::cout << "GetDistanceToEdge " << " phi = " << phi << std::endl;)
  
  // Check max distance first;
  double best_dist = max_dist; // maximum distance
  int c_x = round(x + cos(phi) * max_dist);
  int c_y = round(y + sin(phi) * max_dist);
  int c_z = z; 
  
  if (ContainsCoordinate(n, c_x, c_y, c_z) && ContainsCell(n,c_x,c_y,c_z)) {
    return DBL_MAX; // large distance ... 
  }
  
  // find distance to edge for phi:
  double dist = best_dist / 2.0;
  double delta_dist = dist / 2.0;
  int b_x = c_x, b_y = c_y, b_z = c_z;
  
  while (delta_dist >= 0.5) {
    
    // find current position to check:
    D(std::cout << "Current distance from center x = " << x << ", y = " << y << ", z = " << z << ", d = " << dist << " phi = " << phi << std::endl;)
    c_x = round(x + cos(phi) * dist);
    c_y = round(y + sin(phi) * dist);
    D(std::cout << "Current position: x = " << c_x << ", y = " << c_y << ", z = " << c_z << "d = " << dist << ", delta d = " << delta_dist << std::endl;)
    D(std::cout << "Current best: x = " << b_x << ", y = " << b_y << ", z = " << b_z << " d = " << best_dist << std::endl;)
    
    // check if position is valid:
    if (ContainsCoordinate(n, c_x, c_y, c_z)) { 
      D(std::cout << " => Valid  position." << std::endl;)
      
      // a valid position
      bool contains_cell = ContainsCell(n, c_x, c_y, c_z);
      
      if (check_not_empty_loc && !contains_cell) { // check this is not a random empty position by moving
                                                // further away and checking that grids are empty
        double dist_add = 0;
        while (dist_add < 7.5 && !contains_cell) { // ~5 grid points along diag?
          dist_add += 0.5;
          int c_x_check = round(x + cos(phi) * (dist+dist_add));
          int c_y_check = c_y = round(y + sin(phi) * (dist+dist_add));
          if (ContainsCoordinate(n, c_x_check, c_y_check, c_z)) {
            contains_cell = ContainsCell(n, c_x_check, c_y_check, c_z);
          }
        }
      }
      
      if (contains_cell) { 
        D(std::cout << " => Contains cell." << std::endl;)
        // still containing a cell (move further away):
        dist += delta_dist;
      } else {                  
        D(std::cout << " => Empty." << std::endl;)

        if (dist <= best_dist) {      
          // but maybe a better distance?
          best_dist = dist;
          b_x = c_x;
          b_y = c_y;
          b_z = c_z;
        }
        
        // does not contain a cell (move closer):
        dist -= delta_dist;
        
      }
    } else { 
      D(std::cout << " => Invalid  position." << std::endl;)
      // not a valid position (move closer)
      dist -= delta_dist;
    }
    
    // next step half distance difference
    delta_dist = delta_dist / 2.0;
  }
  
  D(std::cout << "Distance to edge along phi = " << (phi) << " is d = " << dist << "(x = " << c_x - CENTER(mvSizeX[n]) << ", y = " << c_y - CENTER(mvSizeY[n]) << ")" << std::endl;)
    
  x = b_x;
  y = b_y;
  z = b_z;
  
  return best_dist;
}

  
// find position of closest edge:
double Universe::FindClosestEdge(int n, int &x, int &y, int &z, double &phi_ret, double max_dist) const {

  D(std::cout << "Finding closest edge for x=" << x << ", y=" << y << ", z=" << z << "." << std::endl;)
  
  // code does not work for 3d at the moment
  if (mvSizeZ[n] > 1)
    Rcpp::stop("Universe::FindClosestEdge only works for 2d simulations ...");
  
  // code does not work for 1d at the moment
  if (mvSizeX[n] != mvSizeY[n])
    Rcpp::stop("Universe::FindClosestEdge only works for 2d simulations ...");
  
  // First guess for a good direction: opposite of center
  boost::random::uniform_real_distribution<double> rand_offset_phi(-PI/0.025, PI/0.025);
  double phi = atan2(y - CENTER(mvSizeY[n]), x - CENTER(mvSizeX[n])) + rand_offset_phi(rng);
  D(std::cout << "The arc tangent for (x=" << x << ", y=" << y << ") is " << (phi) << " degrees" << std::endl;)
  
  // Do a quick sweep (10 lines) and give up if all edges are taken up
  bool found_empty_position = false;
  double max_dist_sweep = max_dist * 1.5 + 5.0;
  double phi_os[10] = {0.0,PI*0.2,PI*-0.2,PI*0.4,PI*-0.4,PI*0.6,PI*-0.6,PI*0.8,PI*-0.8,PI};
  double phi_o = phi;
  
  for (int i = 0; i < 10; i++) {

      int c_x = round(x + cos(phi_o+phi_os[i]) * max_dist_sweep);
      int c_y = round(y + sin(phi_o+phi_os[i]) * max_dist_sweep);
      int c_z = z;
      
    if (!ContainsCoordinate(n,c_x,c_y,c_z) || !ContainsCell(n,c_x,c_y,c_z)) {
      found_empty_position = true;
      phi = phi_o+phi_os[i]; // use the found edge as start point ... 
      break;
    }
  }

  if (!found_empty_position) { // give up  
    D(std::cout << "Giving up on it." << std::endl;)
    return max_dist + 1.0;
  }

  // find distance to edge for inital value of phi:
  int n_x = x, n_y = y, n_z = z; // position of closest edge 
  double dist = GetDistanceToEdge(n, n_x, n_y, n_z, phi, max_dist + 5.0);

  
  // now test different angles:
  double delta_phi = PI;

  do {

    double last_phi = phi;
    
    for (double i = -1.0; i <= 1.0; i = i + 1.0) {

      D(std::cout << " Multiplicator:  i = " << i << std::endl;)

      // search again for the distance to edge along phi
      double c_phi = last_phi + i * delta_phi;;
      int c_n_x = x, c_n_y = y, c_n_z = z;
      double c_dist = GetDistanceToEdge(n, c_n_x, c_n_y, c_n_z, c_phi, dist + 1.0);
      D(std::cout << " => Tested new phi = " << c_phi << " is d = " << c_dist << "(x = " << c_n_x << ", y = " << c_n_y << ")" << std::endl;)

      // update distance and phi if current one is better
      if (c_dist < dist) {
        D(std::cout << " ** New better distance **" << std::endl;)
        dist = c_dist;
        phi = c_phi;
        n_x = c_n_x;
        n_y = c_n_y;
        n_z = c_n_z;
      }
    }

    delta_phi = delta_phi / 2.0;

  D(std::cout << "New distance to edge along phi = " << (phi) << " is d = " << dist << "(x = " << n_x << ", y = " << n_y << ")" << std::endl;)

  } while (delta_phi >= 2 * PI * 5.0 / 360.0);
  
  // set to best values found:
  phi_ret = phi;
  x = n_x; y = n_y; z = n_z;
  return dist;
}
