#pragma once
#include <cstdint>
#include <fstream>
#include <iostream>

#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>

//phenotype ID (pid) pair, storing phenotype size and index within that per-size vector
using Phenotype_ID = std::pair<uint8_t,uint16_t>;

//pid printing method
inline std::ostream & operator<<(std::ostream& os,Phenotype_ID& pid) {
  os <<"pID: "<< +pid.first<<", " << +pid.second;
  return os;
}

//standard phenotype IDs (pid) 
constexpr Phenotype_ID NULL_pid{0,0}, UNBOUND_pid{255,0};


struct Phenotype {
  
  //instance variables
  uint8_t dx,dy;
  std::vector<uint8_t> tiling;
  
  //static variables
  inline static bool FREE_POLYOMINO=true;
  inline static uint8_t DETERMINISM_LEVEL=3;

  //constructors
  Phenotype(void) : dx{1}, dy{1}, tiling{1} {}
  Phenotype(uint8_t tdx, uint8_t tdy, const std::vector<uint8_t>& ttiling) : dx{tdx}, dy{tdy}, tiling{ttiling} {}
  
  //utility methods
  bool operator==(const Phenotype& rhs) {return this->dx==rhs.dx && this->dy==rhs.dy && this->tiling==rhs.tiling;}
  inline friend std::ostream& operator<<(std::ostream& out, Phenotype& phen);
    
};

//phenotype printing method
inline std::ostream & operator<<(std::ostream& os,Phenotype& phen) {
  os << "Phenotype dx: " << +phen.dx<<", dy:  " << +phen.dy << "\n";
  for(uint8_t y=0;y<phen.dy;++y) {
    for(uint8_t x=0;x<phen.dx;++x)
      os << +phen.tiling[y*phen.dx+x] << " ";
    os << "\n";
  }
  return os;
}

//rotates polyomino pi/2 clockwise
inline void ClockwiseRotation(Phenotype& phen) {
  std::vector<uint8_t> swapper;
  swapper.reserve(phen.tiling.size());
  
  for(uint8_t column=0;column<phen.dx;++column) 
    for(uint8_t row=phen.dy;row!=0;--row)  
      swapper.emplace_back(phen.tiling[(row-1)*phen.dx+column]);
  
  std::swap(phen.dx,phen.dy);
  phen.tiling=swapper;
}

//reflects polyomino over axis, and reflects tile orientation if required
inline void ChiralFlip(Phenotype& phen) {
  
  //reverse each row to reflect polyomino
  for(auto iter=phen.tiling.begin();iter!=phen.tiling.end();iter+=phen.dx)
    std::reverse(iter,iter+phen.dx);

  //get chiral orientation after flip
  if(Phenotype::DETERMINISM_LEVEL==3)
    for(uint8_t& element : phen.tiling)
      if(element && element%2==0)
	element+=-(element-1)%4+((element-1)%4+2)%4;
}

//relabel the polyomino representation to be minimal
inline void MinimizePhenRep(std::vector<uint8_t>& tiling) {

  //if single tile or only shape-dependent, relabel all as 1s
  if(tiling.size() == 1 || Phenotype::DETERMINISM_LEVEL == 1) {
    for(uint8_t& t : tiling)
      t = t > 0;
    return;
  }

  //offset values to relabel
  for(uint8_t& t:tiling)
    t+=128*(t!=0);

  //relabel tiles to new minimum swap value using 255 as an intermediate
  uint8_t swap_count=1;  
  for(std::vector<uint8_t>::iterator t_iter=std::find_if(tiling.begin(),tiling.end(),[](const int s) { return s>0; });t_iter!=tiling.end();) {
    const uint8_t static_swap=*t_iter;
    for(uint8_t cyclic = 0; cyclic< (Phenotype::DETERMINISM_LEVEL==3?4:1); ++cyclic) {
      const uint8_t pre_swap=Phenotype::DETERMINISM_LEVEL==3 ? (static_swap-(static_swap-1)%4)+((static_swap-1)%4+cyclic)%4 : static_swap;
      std::replace(tiling.begin(),tiling.end(),swap_count,uint8_t(255));
      std::replace(tiling.begin(),tiling.end(),pre_swap,swap_count);
      std::replace(tiling.begin(),tiling.end(),uint8_t(255),pre_swap);
      ++swap_count;
    }  
    t_iter=std::find_if(tiling.begin(),tiling.end(),[swap_count](const int s) { return (s>0 && (Phenotype::DETERMINISM_LEVEL==3?s-(s-1)%4:s) >= swap_count); });
  }
}

//find minimum relabelling over all rotations and chiral flips
inline void GetMinPhenRepresentation(Phenotype& phen) {
  std::vector< std::vector<uint8_t> > min_tilings;

  //always prefer wider than taller, rotate if necessary
  if(phen.dy > phen.dx)
    ClockwiseRotation(phen);

  //iterate through chiral flips
  bool to_flip=true;
  while(true) {
    MinimizePhenRep(phen.tiling);
    min_tilings.emplace_back(phen.tiling);

    //if dy != dx, then only need to test a single rotation by pi
    if(phen.dy != phen.dx) {
      std::reverse(phen.tiling.begin(),phen.tiling.end());
      MinimizePhenRep(phen.tiling);
      min_tilings.emplace_back(phen.tiling);
    }
    //need to test all four unique rotations
    else {
      for(uint8_t rot=0;rot<3;++rot) {
        ClockwiseRotation(phen);
        MinimizePhenRep(phen.tiling);
        min_tilings.emplace_back(phen.tiling);
      }
    }
    
    //flip only if first iteration and free polyomino
    if(to_flip && Phenotype::FREE_POLYOMINO)
      ChiralFlip(phen);
    else
	break;

    to_flip=false;
  }

  //set tiling to smallest representation
  phen.tiling=*std::min_element(min_tilings.begin(),min_tilings.end());
}

//return the minimal phenotype given a vector of assembled tiles
inline Phenotype GetPhenotypeFromGrid(std::vector<int8_t>& placed_tiles) {
  std::vector<int8_t> x_locs, y_locs,tile_vals;
  x_locs.reserve(placed_tiles.size()/3);y_locs.reserve(placed_tiles.size()/3);tile_vals.reserve(placed_tiles.size()/3);

  //seperate out the information in the vector
  for(std::vector<int8_t>::iterator check_iter = placed_tiles.begin();check_iter!=placed_tiles.end();check_iter+=3) {
    x_locs.emplace_back(*check_iter);
    y_locs.emplace_back(*(check_iter+1));
    tile_vals.emplace_back(*(check_iter+2));
  }

  //get polyomino extent
  std::vector<int8_t>::iterator x_left,x_right,y_top,y_bottom;
  std::tie(x_left,x_right)=std::minmax_element(x_locs.begin(),x_locs.end());
  std::tie(y_bottom,y_top)=std::minmax_element(y_locs.begin(),y_locs.end());
  uint8_t dx=*x_right-*x_left+1,dy=*y_top-*y_bottom+1;
  std::vector<uint8_t> spatial_grid(dx*dy);
  
  //assign polyomino tile details based on level of determinism
  for(uint16_t tileIndex=0;tileIndex<x_locs.size();++tileIndex) {
    uint8_t tile_detail=0;
    switch(Phenotype::DETERMINISM_LEVEL) {
    case 1:
      tile_detail=tile_vals[tileIndex] > 0 ? 1 : 0;
      break;
    case 2:
      tile_detail=tile_vals[tileIndex] > 0 ? (tile_vals[tileIndex]-1)/4+1 : 0;
      break;
    case 3:
	[[fallthrough]];
    default:
      tile_detail=tile_vals[tileIndex];
    }
    spatial_grid[(*y_top-y_locs[tileIndex])*dx + (x_locs[tileIndex]-*x_left)]=tile_detail;
  }

  //find the minimum representation for consistency
  Phenotype phen{dx,dy,spatial_grid};
  GetMinPhenRepresentation(phen);
  return phen;
}

//main structure to record information on phenotypes of interest, as well as properties determining assembly
struct PhenotypeTable {

  //if fixed, the table will not add new phenotypes to the table
  bool FIXED_TABLE = false;

  //parameters controlling how many repeated builds are assembled, and the fraction of assembles required to be considered deterministic
  double UND_threshold;
  uint16_t phenotype_builds;

  PhenotypeTable() {
    UND_threshold = 1;
    phenotype_builds = 10;
  }
  PhenotypeTable(uint16_t B) {
    UND_threshold = 1;
    phenotype_builds = B;
  }
  PhenotypeTable(double T) {
    UND_threshold = T;
    phenotype_builds = 10;
  }
  PhenotypeTable(uint16_t B, double T) {
    UND_threshold = T;
    phenotype_builds = B;
  }

  //map with keys given by size and vectors holding the phenotype. The pid can be found as the (map key,index of phenotype in that vector)
  //temporary maps also need to track the count, to later see if common enough to add to table
protected:
  std::unordered_map<uint8_t,std::vector<Phenotype>> undiscovered_phenotypes;
  std::unordered_map<uint8_t, std::vector<uint16_t>> undiscovered_phenotype_counts;
public:
  std::unordered_map<uint8_t,std::vector<Phenotype> > known_phenotypes;
  

  
  inline Phenotype_ID GetPhenotypeID(Phenotype& phen) {
    
    //get phenotype size
    uint8_t phenotype_size=std::count_if(phen.tiling.begin(),phen.tiling.end(),[](const int c){return c != 0;});
    
    //compare against existing table entries, return pid if it exists
    for(size_t phenotype_index=0; phenotype_index != known_phenotypes[phenotype_size].size();++phenotype_index)
      if(phen==known_phenotypes[phenotype_size][phenotype_index]) 
        return Phenotype_ID{phenotype_size,phenotype_index};

    //if not match found and table is fixed, return a default value
    if(FIXED_TABLE)
      return NULL_pid;
  
    //compare against temporary table entries
    uint16_t new_phenotype_index=0;
    while(new_phenotype_index<undiscovered_phenotypes[phenotype_size].size()) { 
      if(phen==undiscovered_phenotypes[phenotype_size][new_phenotype_index])
        goto found_phen;
      ++new_phenotype_index;
    }
    
    //not in temporary either, add to temporary
    undiscovered_phenotypes[phenotype_size].emplace_back(phen);
    undiscovered_phenotype_counts[phenotype_size].emplace_back(0);

    //increment the count for the temporary and return its temporary pid
  found_phen:
    ++undiscovered_phenotype_counts[phenotype_size][new_phenotype_index];
    return Phenotype_ID{phenotype_size,known_phenotypes[phenotype_size].size()+new_phenotype_index+phenotype_builds};
  }

  //helper function to reset temporary trackers
  inline void ClearIncomplete() {
    undiscovered_phenotypes.clear();
    undiscovered_phenotype_counts.clear();
  }

  //relabel temporary pids if sufficiently common to known_phenotypes
  inline void RelabelPIDs(std::vector<Phenotype_ID >& pids,bool clear=false) {
    const uint16_t thresh_val=std::ceil(UND_threshold*phenotype_builds);
    for(auto& kv : undiscovered_phenotype_counts) {
      const size_t table_size=known_phenotypes[kv.first].size(); 
      for(size_t nth=0; nth<kv.second.size(); ++nth)
	//if common enough, change pid and add to main map
        if(kv.second[nth] >= thresh_val) {
          std::replace(pids.begin(),pids.end(),Phenotype_ID{kv.first,table_size+phenotype_builds+nth},Phenotype_ID{kv.first,known_phenotypes[kv.first].size()});
          known_phenotypes[kv.first].emplace_back(undiscovered_phenotypes[kv.first][nth]);
        }
    }
    if(clear)
      ClearIncomplete();
  }

  //relabel templated map, same concept as RelabelPIDs
  template<typename map_val>
  inline void RelabelMaps(std::map<Phenotype_ID, map_val>& map, bool clear=false) {
 const uint16_t thresh_val=std::ceil(UND_threshold*phenotype_builds);
    for(auto& kv : undiscovered_phenotype_counts) {
      const size_t table_size=known_phenotypes[kv.first].size(); 
      size_t back_counts=1;
      for(size_t nth=1; nth<=kv.second.size(); ++nth)
        if(kv.second[kv.second.size()-nth] >= thresh_val) {
	auto nh = map.extract(Phenotype_ID{kv.first,table_size+phenotype_builds-back_counts});
	nh.key() = {kv.first,known_phenotypes[kv.first].size()-back_counts};
	map.insert(std::move(nh));
	++back_counts;
        }
      
    }
    if(clear)
      ClearIncomplete();


  }

  
  //find the frequency of each phenotype in a vector of pids, stripping pids that were below a threshold and adding a rare pid
  inline std::map<Phenotype_ID,uint16_t> PhenotypeFrequencies(std::vector<Phenotype_ID >& pids,const Phenotype_ID RARE_pid=NULL_pid, bool allow_existing=false) { 
    const uint16_t thresh_val=std::ceil(UND_threshold*phenotype_builds);
    std::map<Phenotype_ID, uint16_t> ID_counter;

    //get count of each pid in the vector
    for(auto& pid : pids)
      ++ID_counter[pid];

    //erase pids that were not common enough, or allow them to survive if allow_existing is true
    const size_t all_size= ID_counter.size();
    for(auto map_it = ID_counter.begin(); map_it != ID_counter.end();) {
      if(map_it->second >= thresh_val || (allow_existing && map_it->first.second < known_phenotypes[map_it->first.first].size())) 
        ++map_it;
      else
	map_it=ID_counter.erase(map_it);
        
    }
    if(ID_counter.size()!=all_size)
      ++ID_counter[RARE_pid];

    return ID_counter;    
  }

  //fill a fixed or unfixed table with details from a give file name
  inline void LoadTable(std::string f_name) {
    std::ifstream fin(f_name);
    std::string temp;
    while (std::getline(fin, temp)) {
      std::istringstream buffer(temp);
      std::vector<uint8_t> vec{std::istream_iterator<double>(buffer), std::istream_iterator<double>()};
      
      uint8_t phenotype_size=std::count_if(vec.begin()+2,vec.end(),[](const int c){return c != 0;});
      known_phenotypes[phenotype_size].emplace_back(vec[0],vec[1],std::vector<uint8_t>{vec.begin()+2,vec.end()});
    }
  }

  //print a table to a given file name, with details of pid and phenotype
  inline void PrintTable(std::string f_name) {
    std::ofstream fout(f_name);
    for(auto known_phens : known_phenotypes) {
      uint16_t n_phen=0;
      for(Phenotype known : known_phens.second) {
        fout<<+known_phens.first<<" "<<+n_phen++<<" "<<+known.dx<<" "<<+known.dy<<" ";
        for(uint8_t tile : known.tiling)
          fout<<+tile<<" ";
        fout<<"\n";
      }
    }
  }
};
