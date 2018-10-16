//#include "core_genotype.hpp"
#include <fstream>
#include <iostream>

#include <unordered_map>
#include <map>

#include <cmath>
#include <cstdint>

#include <vector>
#include <set>
#include <tuple>

#include <algorithm>
#include <utility>
#include <numeric>

/*! free vs one-sided polyominoes and tile vs orientation determinism */
constexpr uint8_t FREE_POLYOMINO = true ? 2 : 1, DETERMINISM_LEVEL=3;
/*! 
  (true) free polyominoes are not chirally distinct
  (false) one-sided polyominoes are

determinism levels as follows:
    shape       : 1
    tile        : 2
    orientation : 3
*/

/*! phenotype definitions */
struct Phenotype {
  uint8_t dx=1,dy=1;
  std::vector<uint8_t> tiling{1};
};

typedef std::pair<uint8_t,uint16_t> Phenotype_ID;
constexpr Phenotype_ID NULL_pid{0,0};

/*! print phenotype to stdout */
void PrintShape(Phenotype& phen);

/*! phenotype rotations */
void ClockwiseRotation(Phenotype& phen);
void ClockwisePiRotation(Phenotype& phen);
void ChiralFlip(Phenotype& phen);

/*! phenotype comparisons*/
bool ComparePolyominoes(Phenotype& phen1, const Phenotype& phen2);
void MinimizePhenRep(std::vector<uint8_t>& tiling);
void GetMinPhenRepresentation(Phenotype& phen);

namespace model_params
{
  extern bool FIXED_TABLE;
  extern double UND_threshold;
  extern uint16_t phenotype_builds;
}

struct PhenotypeTable {
  std::unordered_map<uint8_t,std::vector<Phenotype> > known_phenotypes;
  
  Phenotype_ID GetPhenotypeID(Phenotype& phen);

  void RelabelPhenotypes(std::vector<Phenotype_ID >& pids);

  /* Count each ID frequency */
  std::map<Phenotype_ID,uint16_t> PhenotypeFrequencies(std::vector<Phenotype_ID >& pids, bool& rare_phenotypes);
        
  void PrintTable(std::ofstream& fout);
private:
  std::unordered_map<uint8_t,std::vector<Phenotype>> undiscovered_phenotypes;
  std::vector<size_t> new_phenotype_xfer;
  std::unordered_map<uint8_t, std::vector<uint16_t>> undiscovered_phenotype_counts;

};