/*
 * Declarations for Deme class to evolve a genetic algorithm for the
 * travelling-salesperson problem.  A deme is a population of individuals.
 */

#include "chromosome.hh"
#include "deme.hh"
#include <algorithm>
//#include <numeric>
#include <chrono>
#include <cassert>

using vec_size_t = std::vector<Chromosome*>::size_type;

Chromosome* mut_decider(double rand, Chromosome* chromosome_ptr, double mut_rate);
double frac(std::default_random_engine generator){ return double(generator()) / double(generator.max()); }

//Deme Members

// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
// Also seeds the generator idk if this is important but i added it
Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate){
	for(unsigned i = 0; i < pop_size; ++i){
		Chromosome* chromosome_ptr = new Chromosome(cities_ptr);
		pop_.push_back(chromosome_ptr);
		mut_rate_ = mut_rate;
	}
	std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
	generator_ = generator;
}

// Clean up as necessary
Deme::~Deme(){
	for(Chromosome* chromosome_ptr : pop_){
		delete chromosome_ptr;
	}
}

// Evolve a single generation of new chromosomes, as follows:
// We select pop_size/2 pairs of chromosomes (using the select() method below).
// Each chromosome in the pair can be randomly selected for mutation, with
// probability mut_rate, in which case it calls the chromosome mutate() method.
// Then, the pair is recombined once (using the recombine() method) to generate
// a new pair of chromosomes, which are stored in the Deme.
// After we've generated pop_size new chromosomes, we delete all the old ones.
void Deme::compute_next_generation(){
	std::vector< std::pair<Chromosome*,Chromosome*> > chrom_pairs;
	for(vec_size_t i = 0; i < pop_.size() / 2; i++){
		chrom_pairs.push_back(std::pair<Chromosome*, Chromosome*>(select_parent(), select_parent()));
	}
	for(std::pair<Chromosome*, Chromosome*> pair : chrom_pairs){
		mut_decider(frac(generator_), pair.first, mut_rate_);
		mut_decider(frac(generator_), pair.second, mut_rate_);
		pair = pair.first->recombine(pair.second);
		pop_.push_back(pair.first);
		pop_.push_back(pair.second);
	}
	int half_pop_size = pop_.size() / 2;
	for(int i = 0; i < half_pop_size; i++){
		delete pop_[0];
		pop_.erase(pop_.begin());
	}
}

Chromosome* mut_decider(double rand, Chromosome* chromosome_ptr, double mut_rate){
	if(rand < mut_rate){
		chromosome_ptr->mutate();
	}
	return chromosome_ptr;
}

// Return a copy of the chromosome with the highest fitness.
// ^ This is literally impossible, the chromosome.hh implementation we were given deletes the copy and assignment constructor. I'm going with the moodle instructions, which say to return a pointer to the best chromosome.
const Chromosome* Deme::get_best() const{
	Chromosome* best = pop_[0];
	if(pop_.size() > 1){
		double best_fitness = best->get_fitness();
		for(vec_size_t i = 1; i < pop_.size(); ++i){
			if(pop_[i]->get_fitness() > best_fitness){
				best = pop_[i];
				best_fitness = best->get_fitness();
			}
		}
	}
	assert(best);
	return best;
}

// Randomly select a chromosome in the population based on fitness and
// return a pointer to that chromosome.
Chromosome* Deme::select_parent(){
	unsigned fit_sum = 0;
	for(Chromosome* chrom : pop_){
		fit_sum += chrom->get_fitness();
	}
	unsigned long selector = generator_() & fit_sum;
	unsigned long fit_check = 0;
	int i = 0;
	while(fit_check < selector){
		fit_check += pop_[i]->get_fitness();
		++i;
	}
	++i;
	return pop_[i];
}

/*
	std::vector<vec_size_t> fps_table(pop_.size());	//table for fitness proportionate selection

	fps_table[0] = pop_[0]->get_fitness();	//initialize first value as fitness of first cromosome
	for(vec_size_t i = 1; i < pop_.size(); ++i){
		fps_table[i] = fps_table[i-1] + pop_[i]->get_fitness();	//each chromosome "takes up" the amount of "space" proportionate to their fitness
	}
	unsigned gend = generator_(); 						std::cout << "gend: " << gend << std::endl;
														std::cout << "t size: " << fps_table.size() << std::endl;
	unsigned max = fps_table[fps_table.size() - 1];		std::cout << "max: " << max << std::endl;
	//vec_size_t rand = (unsigned(generator_()) % unsigned(fps_table[fps_table.size() - 1])); //number between 0 and the sum of all fitnesses, i.e. the value of the last element
	vec_size_t rand(gend % max);

	int i = 0;
	while(fps_table[i] < rand){	//iterates through the table until the random value is within the "space" of the chromosome's fitness on the table
		++i;
	}
	assert(pop_[i]);
	return pop_[i];
}*/
/*
using namespace std;
Chromosome* Deme::select_parent()
{
  double sum = 0, partial_sum = 0;
  for(auto chromo: pop_){
    sum += chromo->get_fitness();
  }
  uniform_int_distribution<int> distr(0, sum);

  int benchmark = distr(generator_);

  for(auto chromo: pop_)
  {
    assert(chromo);
    double fitness = chromo->get_fitness();
    if(partial_sum + fitness < benchmark){
		partial_sum += fitness; }
    else{
      auto order = chromo->get_ordering();

      return chromo;
    }
  }
}
*/