/*
 * Declarations for Deme class to evolve a genetic algorithm for the
 * travelling-salesperson problem.  A deme is a population of individuals.
 */

#include "chromosome.hh"
#include "deme.hh"
#include <algorithm>
//#include <numeric>

// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].

Chromosome* mut_decider(std::default_random_engine& generator, Chromosome* chromosome_ptr, double mut_rate);

Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate){
	for(int i = 0; i < int(pop_size); ++i){
		Chromosome* chromosome_ptr = new Chromosome(cities_ptr);
		pop_.push_back(chromosome_ptr);
		//pop_.push_back(std::make_shared<Chromosome>(cities_ptr));
		mut_rate_ = mut_rate;
	}
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
	for(int i = 0; i <= int(pop_.size()) / 2; ++i){
		Chromosome* parent_1 = mut_decider(generator_, select_parent(), mut_rate_);
		Chromosome* parent_2 = mut_decider(generator_, select_parent(), mut_rate_);

		std::pair<Chromosome*, Chromosome*> children = parent_1->recombine(parent_2);
		delete parent_1;
		delete parent_2;
		parent_1 = children.first;
		parent_2 = children.second;
	}
}

Chromosome* mut_decider(std::default_random_engine& generator, Chromosome* chromosome_ptr, double mut_rate){
	if(generator.max() / generator() < mut_rate){
		chromosome_ptr->mutate();
	}
	return chromosome_ptr;
}

// Return a copy of the chromosome with the highest fitness.
// ^ This is literally impossible, the chromosome.hh implementation we were given deletes the copy and assignment constructor. I'm going with the moodle instructions, which say to return a pointer to the best chromosome.
const Chromosome* Deme::get_best() const{
	//return std::max_element(/*construct an array of the get_fitness of eawch element of pop_*/)
	//not sure if this ^ would be faster or slower; I'll compare later if I have time

	Chromosome* ret = pop_[0];
	//for(chromosome_ptr_t chromosome_ptr : pop_){ //more spicy code i can't use
	if(pop_.size() > 1){
		double current_fitness = ret->get_fitness();	//saves on computing the fitness of the current best option every time they're compared later
		for(int i = 1; i < int(pop_.size()); ++i){	//is it faster to assign pop_.size() to a variable so i don't have to access it like this? i checked the implementation for vectors and it's calculated by subtraction of the begin and end iterators
			if(pop_[i]->get_fitness() > current_fitness){
				ret = pop_[i];
				current_fitness = ret->get_fitness();
			}
		}
	}
	return ret;
}

// Randomly select a chromosome in the population based on fitness and
// return a pointer to that chromosome.
Chromosome* Deme::select_parent(){
	std::vector<int> fps_table(pop_.size());	//table for fitness proportionate selection
	fps_table[0] = pop_[0]->get_fitness();	//initialize first value as fitness of first cromosome
	for(int i = 1; i < int(pop_.size()); ++i){
		fps_table[i] = fps_table[i-1] + pop_[i]->get_fitness();	//each chromosome "takes up" the amount of "space" proportionate to their fitness
	}

	int rand = generator_() % *(fps_table.end() - 1); //number between 0 and the sum of all fitnesses, i.e. the value of the last element	//the right side of the % operator is a little cursed

	int i = 0;	// the place of the chromosome in which
	while(fps_table[i] < rand){	//iterates through the table until the random value is within the "space" of the chromosome's fitness on the table
		++i;
	}
	++i;	// sets i to the first chromosome with a fitness space larger than the random number
	return pop_[i];
}
