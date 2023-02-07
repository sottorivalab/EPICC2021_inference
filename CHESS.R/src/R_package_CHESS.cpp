#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "boost/multi_array.hpp"
#include <iostream>
#include <cassert>
#include <time.h>
#include "Universe.h"
#include "extern_global_variables.h"
// [[Rcpp::depends(BH)]]
#include "Rcpp.h"

boost::random::mt19937_64 rng;


RCPP_MODULE(CHESS_Universe_object) {
  using namespace Rcpp;

  class_<Universe>("CHESS_Universe_object")
  .constructor<int, int, int, Rcpp::NumericMatrix, int, int>()
  .method("TakeSample", &Universe::TakeSampleBox)
  .method("SingleCellTree", &Universe::SingleCellTree)
  .method("SingleCellTreeWithNodeIds", &Universe::SingleCellTreeWithNodeIds)
  .method("SingleCellTreeN", &Universe::SingleCellTreeN)
  .method("ScaleMutationTree", &Universe::ScaleMutationTree)
  
  .method("CellType", &Universe::GetCellType)
  .method("CellMutations", &Universe::GetMutations)
  .method("CellMutationBurden", &Universe::GetMutationBurden)
  .method("ClosestEdge", &Universe::FindClosestEdgeRcpp)
  .method("GetDistanceToEdge", &Universe::GetDistanceToEdgeRcpp)
  .method("GetDistanceToEdgeV2", &Universe::GetDistanceToEdgeRcppV2)
  
  .method("RunSimulation", &Universe::RunSimulation)
  .property("History", &Universe::GetHistory, "Record of cell numbers over time")
  .property("Snapshots", &Universe::GetSnapshots, "Record of space at past time points")
  .property("SnapshotTimes", &Universe::GetSnapshotTimes, &Universe::SetSnapshotTimes, "Time points at which do record the space" )
  .property("Time", &Universe::GetTime, "Current universe time")
  .property("GillespieTime", &Universe::Time, "Current universe time")
  .property("Transformations", &Universe::GetTransformations, "Transformations" )
  .property("CellCountsPerType", &Universe::CellCountsPerType, "Number of cells per type" )
  .property("CellCount", &Universe::CellCounts, "Total cell count" )
  .property("ExploreLocally", &Universe::GetExploreLocally, &Universe::SetExploreLocally, "T" )
  .property("AltPush", &Universe::GetAltPush, &Universe::SetAltPush, "T" )
  .property("RecordHistory", &Universe::GetRecordHistory, &Universe::SetRecordHistory, "T" )
  .property("GenerationTree", &Universe::GetReturnGenerationTree, &Universe::SetReturnGenerationTree, "T" )
  .property("ClonalMutations", &Universe::GetClonalMutations, &Universe::SetClonalMutations, "T" )
  .property("StopPopulationSize", &Universe::GetStopPopulationSize, &Universe::SetStopPopulationSize, "Population size at which simulations are stopped.")
  ;
}
