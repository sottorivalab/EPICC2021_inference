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

#include "CellSample.h"
#include "extern_global_variables.h"
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/beta_distribution.hpp>
#include <math.h>       /* pow */
#include <stdlib.h>     /* abs */
#include <iomanip>      // std::setw
#include <fstream>
#include <iostream>
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
// [[Rcpp::depends(BH)]]
#include "Rcpp.h"
#pragma clang diagnostic pop

#define DEPTH_OVERDISPERSION 0.08


namespace ngs_simulation {
  int simulate_depth(double depth, int depth_model) {
    // Sample depth depending upon selected depth model:
    int seq_depth;
    switch(depth_model) {

      case 1: {// Poisson distributed depth model.
        boost::random::poisson_distribution<int> dist_depth1(depth);
        seq_depth = dist_depth1(rng);
        break;
      }

      case 2: {// overdispersed beta binomial model
        double mu = 0.6;
        double rho = DEPTH_OVERDISPERSION;
        double sh1 = mu * (1.0 / rho - 1.0); // scale param 1 (beta)
        double sh2 = (mu - 1.0) * (rho - 1.0) / rho; // scale param 2 (beta)

        // beta component
        boost::random::beta_distribution <double> beta_component_depth2(sh1, sh2);
        double cbd = beta_component_depth2(rng);

        // binomial component
        boost::random::binomial_distribution<int> dist_depth2(depth / mu, cbd);
        seq_depth = dist_depth2(rng);
        break;
      }

      case 3: {// fixed depth model
        seq_depth = depth;
        break;
      }

      default: {
        Rcpp::Rcerr << "Error in function AddNoiseToVaf:" << std::endl;
        Rcpp::Rcerr << "    Unknow depth_noise_model." << std::endl;
        Rcpp::stop("Unknown depth model.\n");
        break;
      }
    }
    return seq_depth;
  }


  int simulate_sequencing(int seq_depth, double seq_vaf) {
    // Do bionomial sampling at vaf:
    boost::random::binomial_distribution<int> binom(seq_depth, seq_vaf);
    int seq_reads = binom(rng);
    return seq_reads;
  }
}

// Base class TumourSample /////////////////////////////////////////////////////

CellSample::CellSample()
  : mTotalNumberCells(0)
  {}

void CellSample::AddMutation(std::string mutation_id,
                             unsigned int clone_id,
                             unsigned int number_mutated)
{
  mvMutationId.push_back(mutation_id);
  mvCloneId.push_back(clone_id);
  mvNumberMutatedCells.push_back(number_mutated);
}

void CellSample::SetSampledCellNumber(unsigned int n_cells) {
  mTotalNumberCells = n_cells;
}


SequencingResult CellSample::Sequence(double depth, int depth_model,
                                      int min_reads, double min_vaf, 
                                      double purity) const
{
  SequencingResult results;

  for (std::vector<unsigned int>::size_type i = 0;
       i < mvNumberMutatedCells.size();
       i++)
  {

    // Simulate ngs:
    double ccf = mvNumberMutatedCells[i] * 1.0 / mTotalNumberCells;
    double exp_vaf = ccf / 2.0 * purity;
    int sim_depth = ngs_simulation::simulate_depth(depth, depth_model);
    int sim_reads = ngs_simulation::simulate_sequencing(sim_depth, exp_vaf);

    // Add to result object if more equal min reads:
    if (sim_reads >= min_reads && ((sim_reads * 1.0 / sim_depth) >= min_vaf || (sim_depth == 0 && min_reads == 0))){
      results.AddMutation(mvMutationId[i], mvCloneId[i], sim_reads, sim_depth, ccf);
    }
  }

  return results;
}



// Derived class SequencingResult //////////////////////////////////////////////


SequencingResult::SequencingResult()
  : mNumberMutations(0)
  {}

void SequencingResult::AddMutation(std::string mutation_id,
                                   unsigned int clone_id,
                                   unsigned int alt_count,
                                   unsigned int depth, 
                                   double ccf)
{
  mNumberMutations++;
  mvMutationId.push_back(mutation_id);
  mvCloneId.push_back(clone_id);
  mvAltCounts.push_back(alt_count);
  mvSequencedDepth.push_back(depth);
  mvCCF.push_back(ccf);
}

std::vector <unsigned int> SequencingResult::AltVector() const {
  return mvAltCounts;
}

std::vector <unsigned int> SequencingResult::CloneIdVector() const {
  return mvCloneId;
}

std::vector <unsigned int> SequencingResult::DepthsVector() const {
  return mvSequencedDepth;
}

std::vector <std::string> SequencingResult::MutationIdVector() const {
  return mvMutationId;
}

std::vector <double> SequencingResult::CCFVector() const {
  return mvCCF;
}

