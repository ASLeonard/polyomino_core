#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>
#include <tuple>
#include <set>

#include <array>
#include <random>


extern thread_local std::mt19937 RNG_Engine;

using InteractionPair = std::pair<uint8_t,uint8_t>;

struct PotentialTileSites {
  std::vector<std::pair<InteractionPair,std::array<int8_t,3>>> sites; //array with x position, y position, and rotation
  std::vector<double> strengths;
};


template<class Q> //Curiously recurring template pattern
class PolyominoAssembly {
public:
  inline static bool free_seed=true;

  //strip all subunits which cannot interact with initial subunit
  template<typename T, typename A>
  static void StripNoncodingGenotype(std::vector<T,A>& genotype) {
    std::vector<size_t> coding{0},noncoding(genotype.size()/4-1);
    std::iota(noncoding.begin(), noncoding.end(), 1);
  
    for(size_t c_in=0;c_in<coding.size();++c_in)
      for(uint8_t cface=0;cface<4;++cface)
        for(size_t nc_in=0;nc_in<noncoding.size();++nc_in) {
          for(uint8_t ncface=0;ncface<4;++ncface)
            if(Q::InteractionMatrix(genotype[coding[c_in]*4+cface],genotype[noncoding[nc_in]*4+ncface])){
              coding.emplace_back(noncoding[nc_in]);
              noncoding.erase(noncoding.begin()+nc_in--);
              goto newtile;
            }
        newtile: ;
        }

    for(size_t rm=0;rm<noncoding.size();++rm)
      genotype.erase(genotype.begin()+(noncoding[rm]-rm)*4,genotype.begin()+(1+noncoding[rm]-rm)*4);
  }

  //returns set of all interaction,strength tuples in a genotype
  template <typename T, typename A>
  static std::vector<std::pair<InteractionPair,double> > GetActiveInterfaces(const std::vector<T,A>& genotype) {
    std::vector<std::pair<InteractionPair,double> > edge_pairs;
    for(size_t b1=0;b1<genotype.size();++b1)
      for(size_t b2=b1;b2<genotype.size();++b2) 
        if(double B_ij=Q::InteractionMatrix(genotype[b1],genotype[b2]))
          edge_pairs.emplace_back(InteractionPair{b1,b2},B_ij);
    
    return edge_pairs;
  }

  static inline std::vector<int8_t> AssemblePolyomino(const std::vector<std::pair<InteractionPair,double> > edges, std::set<InteractionPair>& interacting_indices) {
    //uint16_t accumulated_time=0;
    if(edges.empty())
      return {0,0,1};
    auto max_subunit = std::max_element(edges.begin(), edges.end(),[](const auto& left, const auto& right){return left.first.second <  right.first.second;})->first.second;
   
    const int8_t seed = 1+Q::free_seed*4*std::uniform_int_distribution<uint8_t>(0,max_subunit/4)(RNG_Engine);
    const size_t UNBOUND_LIMIT= 12*(max_subunit/4+1)*(max_subunit/4+1);

    std::vector<int8_t> placed_tiles{0,0,seed},growing_perimeter;
    std::vector<double> strengths_cdf;
    PotentialTileSites perimeter_sites;
    
    ExtendPerimeter(edges,seed,0,0,placed_tiles,perimeter_sites);

    while(!perimeter_sites.strengths.empty()) {
      //select new site proportional to binding strength 
      strengths_cdf.resize(perimeter_sites.strengths.size());
      std::partial_sum(perimeter_sites.strengths.begin(), perimeter_sites.strengths.end(), strengths_cdf.begin());
      std::uniform_real_distribution<double> random_interval(0,strengths_cdf.back());
      size_t selected_choice=static_cast<size_t>(std::lower_bound(strengths_cdf.begin(),strengths_cdf.end(),random_interval(RNG_Engine))-strengths_cdf.begin());
      
      auto chosen_site=perimeter_sites.sites[selected_choice];
      /*
        if(true) {
        accumulated_time+=1+std::geometric_distribution<uint16_t>(perimeter_sites.strengths[selected_choice])(RNG_Engine);
        if(accumulated_time>2)
	break;
        }
      */
      //place new tile 
      placed_tiles.insert(placed_tiles.end(),chosen_site.second.begin(),chosen_site.second.end());
      interacting_indices.insert(chosen_site.first);
      if(placed_tiles.size()>UNBOUND_LIMIT)
        return {};
      
      //remove all further options in same tile location
      auto [f_x, f_y, f_t] = chosen_site.second;
      for(size_t cut_index=0;cut_index<perimeter_sites.strengths.size();) {
        if(f_x==perimeter_sites.sites[cut_index].second[0] && f_y==perimeter_sites.sites[cut_index].second[1]) {
          perimeter_sites.strengths.erase(perimeter_sites.strengths.begin()+cut_index);
          perimeter_sites.sites.erase(perimeter_sites.sites.begin()+cut_index);
        }
        else
          ++cut_index;
      }
      ExtendPerimeter(edges,f_t,f_x,f_y,placed_tiles,perimeter_sites);
    }
    return placed_tiles;
  }
 
  static inline void ExtendPerimeter(const std::vector<std::pair<InteractionPair,double> >& edges,uint8_t tile_detail, int8_t x,int8_t y, std::vector<int8_t>& placed_tiles,PotentialTileSites& perimeter_sites) {
    int8_t dx=0,dy=0,tile=(tile_detail-1)/4,theta=(tile_detail-1)%4;
    for(uint8_t f=0;f<4;++f) {
      switch(f) {
      case 0:dx=0;dy=1;break;
      case 1:dx=1;dy=0;break;
      case 2:dx=0;dy=-1;break;
      case 3:dx=-1;dy=0;break;
      }
      /*! site already occupied, move on */
      for(std::vector<int8_t>::reverse_iterator tile_info=placed_tiles.rbegin();tile_info!=placed_tiles.rend();tile_info+=3) 
        if((x+dx)==*(tile_info+2) && (y+dy)==*(tile_info+1)) 
          goto nextloop;

      {//empty scope to prevent cross-initilisation errors from goto
        uint8_t g_index=static_cast<uint8_t>(tile*4+(f-theta+4)%4);
        std::vector<std::pair<InteractionPair,double> >::const_iterator iter = edges.begin();
        while ((iter = std::find_if(iter, edges.end(),[&g_index](const auto& edge){ return (edge.first.first == g_index || edge.first.second == g_index);})) != edges.end()) {
          int8_t base = iter->first.first==g_index ? iter->first.second : iter->first.first;
          perimeter_sites.sites.emplace_back(iter->first,std::array<int8_t,3>{static_cast<int8_t>(x+dx),static_cast<int8_t>(y+dy),static_cast<int8_t>(base-base%4+((f+2)%4-base%4+4)%4+1)});
          perimeter_sites.strengths.emplace_back(iter->second);
          ++iter;
        }
      }
      /*! find all above threshold bindings and add to new sites */
    nextloop: ; //continue if site already occupied
    }   
  }
  
};
