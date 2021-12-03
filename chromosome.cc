/*
 * Implementation for Chromosome class
 */

#include <algorithm>
#include <cassert>
#include <random>
#include <chrono>
#include <random>

#include "chromosome.hh"

using namespace std;

Cities::permutation_t range(int num); //helper function to find the range of a number

//////////////////////////////////////////////////////////////////////////////
// Generate a completely random permutation from a list of cities
Chromosome::Chromosome(const Cities* cities_ptr)
  : cities_ptr_(cities_ptr),
    order_(random_permutation(cities_ptr->size()))
{
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  generator_ = generator;
  assert(is_valid());
}

//////////////////////////////////////////////////////////////////////////////
// Clean up as necessary
Chromosome::~Chromosome()
{
  assert(is_valid());
}

//////////////////////////////////////////////////////////////////////////////
// Perform a single mutation on this chromosome
void
Chromosome::mutate()
{
  assert(is_valid());

  unsigned pos1 = generator_() % order_.size();
  unsigned pos2 = generator_() % order_.size();
  swap (order_[pos1], order_[pos2]);
}

//////////////////////////////////////////////////////////////////////////////
// Return a pair of offsprings by recombining with another chromosome
// Note: this method allocates memory for the new offsprings
pair<Chromosome*, Chromosome*>
Chromosome::recombine(const Chromosome* other)
{
  assert(is_valid());
  assert(other->is_valid());

  auto child_1 = create_crossover_child(this, other,
                                        generator_() % order_.size(),
                                        generator_() % order_.size());
  auto child_2 = create_crossover_child(other, this,
                                        generator_() % order_.size(),
                                        generator_() % order_.size());
  pair<Chromosome*, Chromosome*> offspring (child_1, child_2);
  return offspring;
}

//////////////////////////////////////////////////////////////////////////////
// For an ordered set of parents, return a child using the ordered crossover.
// The child will have the same values as p1 in the range [b,e),
// and all the other values in the same order as in p2.
Chromosome*
Chromosome::create_crossover_child(const Chromosome* p1, const Chromosome* p2,
                                   unsigned b, unsigned e) const
{
  Chromosome* child = p1->clone(); //this is where new memory is allocated

  // We iterate over both parents separately, copying from parent1 if the
  // value is within [b,e) and from parent2 otherwise
  unsigned i = 0, j = 0;

  for ( ; i < p1->order_.size() && j < p2->order_.size(); ++i) {
    if (i >= b && i < e) {
      child->order_[i] = p1->order_[i];
    }
    else { // Increment j as long as its value is in the [b,e) range of p1
      while (p1->is_in_range(p2->order_[j], b, e)) {
        ++j;
      }
      assert(j < p2->order_.size());
      child->order_[i] = p2->order_[j];
      j++;
    }
  }

  assert(child->is_valid());
  return child;
}

// Return a positive fitness value, with higher numbers representing
// fitter solutions (shorter total-city traversal path).
double
Chromosome::get_fitness() const
{
  double dist = calculate_total_distance();
  assert(dist > 0);
  return(1000 / dist); //smaller the distance is, the larger its fitness will be
}

// A Chromosome is valid if it has no repeated values in its permutation,
// as well as no indices above the range (length) of the chromosome.
// We implement this check with a sort, which is a bit inefficient, but simple
bool
Chromosome::is_valid() const
{
  auto copy = order_;
  sort(copy.begin(), copy.end());
  return (copy == range(order_.size()));
}

// Find whether a certain value appears in a given range of the chromosome.
// Returns true if value is within the specified the range specified
// [begin, end) and false otherwise.
bool
Chromosome::is_in_range(unsigned value, unsigned begin, unsigned end) const
{
  assert(end < order_.size());
  for(unsigned i = begin; i < end; ++i){
    if(order_[i] == value) { return true; }
  }
  return false;
}

Cities::permutation_t range(int num) //helper function
{
  assert(num >= 0);
  Cities::permutation_t v;
  for(int i = 0; i < num; ++i){ v.push_back(i); }
  return v;
}
