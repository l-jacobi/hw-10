/*
 * Declarations for Deme class to evolve a genetic algorithm for the
 * travelling-salesperson problem.  A deme is a population of individuals.
 */

#include "chromosome.hh"
#include "deme.hh"
#include <algorithm>
//#include <numeric>

using vec_size_t = std::vector<Chromosome*>::size_type;

Chromosome* mut_decider(std::default_random_engine& generator, Chromosome* chromosome_ptr, double mut_rate){
	if(generator.max() / generator() < mut_rate){
		chromosome_ptr->mutate();
	}
	return chromosome_ptr;
}

//Deme Members

// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate){	std::cout << std::endl << "deme created" << std::endl;
	for(unsigned i = 0; i < pop_size; ++i){
		Chromosome* chromosome_ptr = new Chromosome(cities_ptr);
		pop_.push_back(chromosome_ptr);
		mut_rate_ = mut_rate;
	}
}

// Clean up as necessary
Deme::~Deme(){	std::cout << std::endl << "deme destroyed" << std::endl;
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
void Deme::compute_next_generation(){	std::cout << std::endl << "computing next generation" << std::endl;
	for(vec_size_t i = 0; i <= pop_.size() / 2; ++i){
		Chromosome* parent_1 = mut_decider(generator_, select_parent(), mut_rate_);
		Chromosome* parent_2 = mut_decider(generator_, select_parent(), mut_rate_);

		std::pair<Chromosome*, Chromosome*> children = parent_1->recombine(parent_2);
		delete parent_1;
		delete parent_2;
		parent_1 = children.first;
		parent_2 = children.second;
	}
}

// Return a copy of the chromosome with the highest fitness.
// ^ This is literally impossible, the chromosome.hh implementation we were given deletes the copy and assignment constructor. I'm going with the moodle instructions, which say to return a pointer to the best chromosome.
const Chromosome* Deme::get_best() const{ std::cout << std::endl << "getting best" << std::endl;
	//return std::max_element(/*construct an array of the get_fitness of eawch element of pop_*/)
	//not sure if this ^ would be faster or slower; I'll compare later if I have time

	Chromosome* best = pop_[0];
	//for(chromosome_ptr_t chromosome_ptr : pop_){ //more spicy code i can't use
	if(pop_.size() > 1){
		double best_fitness = best->get_fitness();	//saves on computing the fitness of the current best option every time they're compared later
		for(vec_size_t i = 1; i < pop_.size(); ++i){	//is it faster to assign pop_.size() to a variable so i don't have to access it like this? i checked the implementation for vectors and it's calculated by subtraction of the begin and end iterators
			if(pop_[i]->get_fitness() > best_fitness){
				best = pop_[i];
				best_fitness = best->get_fitness();
			}
		}
	}
	return best;
}

// Randomly select a chromosome in the population based on fitness and
// return a pointer to that chromosome.
Chromosome* Deme::select_parent(){	std::cout << std::endl << "selecting parent" << std::endl;
	std::vector<int> fps_table(pop_.size());	//table for fitness proportionate selection
	fps_table[0] = pop_[0]->get_fitness();	//initialize first value as fitness of first cromosome
	for(vec_size_t i = 1; i < pop_.size(); ++i){
		std::cout << "fitness: " << pop_[i]->get_fitness() << ", putting into table as ";
		fps_table[i] = fps_table[i-1] + pop_[i]->get_fitness();	//each chromosome "takes up" the amount of "space" proportionate to their fitness
		std::cout << fps_table[i] << std::endl;
	}

	//int rand = generator_() % *(fps_table.end() - 1); //number between 0 and the sum of all fitnesses, i.e. the value of the last element	//the right side of the % operator is a little cursed
	//bad code ^

	//int rand = generator_() % fps_table[fps_table.size() - 1]; //number between 0 and the sum of all fitnesses, i.e. the value of the last element
	int rand = generator_();							std::cout << "rand: " << rand << std::endl;
	auto max_size = fps_table[fps_table.size() - 1];	std::cout << "max_size: " << max_size << std::endl;
	rand = rand % max_size;

	int i = 0;	// the place of the chromosome in which
	while(fps_table[i] < rand){	//iterates through the table until the random value is within the "space" of the chromosome's fitness on the table
		++i;
	}
	++i;	// sets i to the first chromosome with a fitness space larger than the random number
	return pop_[i];
}
