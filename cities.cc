/*
 * Implementation for Travelling-Salesperson Cities class and utilities
 */

#include <algorithm>
#include <fstream>
#include <numeric>
#include <random>
#include <chrono>
#include <cassert>

#include "cities.hh"


//////////////////////////////////////////////////////////////////////////////
Cities::Cities(const std::string& filename)
 : cities_()
{
  std::ifstream inp(filename);
  if (inp.is_open()) {
    inp >> *this;
  }
}

//////////////////////////////////////////////////////////////////////////////
// Given a permutation, return a new Cities object where the order of the
// cities reflects the original order of this class after reordering with
// the given ordering. So for example, the ordering { 1, 0 } simply swaps
// the first two elements (coordinates) in the new Cities object.
Cities
Cities::reorder(const permutation_t& ordering) const
{
  Cities ret;
  for (auto i : ordering) {
    ret.cities_.push_back(cities_[i]);
  }
  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// For a given permutation of the cities in this object,
// compute how long (distance) it would take to traverse all the cities in the
// order of the permutation, and then returning to the first city.
// The distance between any two cities is computed as the Euclidean 
// distance on a plane between their coordinates.
double
Cities::total_path_distance(const permutation_t& ord) const
{
  double sum = 0.;
  for (unsigned i = 1; i <= ord.size(); ++i) {
    const unsigned j = i % ord.size();
    const auto dx = cities_[ord[j]].first - cities_[ord[i-1]].first;
    const auto dy = cities_[ord[j]].second - cities_[ord[i-1]].second;
    sum += std::hypot(dx, dy);
  }
  return sum;
}


//////////////////////////////////////////////////////////////////////////////
std::istream& operator>>(std::istream& is, Cities& cities)
{
  cities.cities_ = Cities::cities_t();
  Cities::coord_t city;

  while (is >> city.first) {
    is >> city.second;
    cities.cities_.push_back(city);
  }

  return is;
}

//////////////////////////////////////////////////////////////////////////////
std::ostream& operator<<(std::ostream& os, const Cities& cities)
{
  for (const auto& c : cities.cities_) {
    os << c.first << "\t" << c.second << "\n";
  }

  return os;
}

//////////////////////////////////////////////////////////////////////////////
Cities::permutation_t
random_permutation(unsigned len)
{
  assert(int(len) > 0);
  Cities::permutation_t out;
  int end = len-1;
  //following 2 lines borrowed from a StackOverflow user: see README of hw 09
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(seed);
  std::uniform_int_distribution<unsigned int> distribution(0, end);
  while(out.size() != len){
    auto num = distribution(generator);
    if(count(out.begin(), out.end(), num) == 0){
      out.push_back(num);
    }
  }
  return out;
}


