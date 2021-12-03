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
double frac(std::default_random_engine generator){
	return double(generator()) / double(generator.max());
}

//Deme Members

// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
// Also seeds the generator idk if this is important but i added it
Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate){	std::cout << std::endl << "####deme created" << std::endl;
	for(unsigned i = 0; i < pop_size; ++i){
		Chromosome* chromosome_ptr = new Chromosome(cities_ptr);
		pop_.push_back(chromosome_ptr);
		mut_rate_ = mut_rate;
	}
	std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
	generator_ = generator;
}

// Clean up as necessary
Deme::~Deme(){	std::cout << std::endl << "####deme destroyed" << std::endl;
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
void Deme::compute_next_generation(){	std::cout << std::endl << "####computing next generation" << std::endl;
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



/*
		Chromosome* parent_1 = mut_decider(frac(generator_), select_parent(), mut_rate_);
		Chromosome* parent_2 = mut_decider(frac(generator_), select_parent(), mut_rate_);

		//parent_1->recombine(parent_2);
		std::pair<Chromosome*, Chromosome*> children = parent_1->recombine(parent_2);
		std::cout << "parent 1: " << parent_1 << ", parent 2: " << parent_2 << std::endl;
		delete parent_1;
		if(parent_1 != parent_2){delete parent_2;}
		std::cout << "parent 1: " << parent_1 << ", parent 2: " << parent_2 << std::endl;
		parent_1 = children.first;
		parent_2 = children.second;
		std::cout << "parent 1: " << parent_1 << ", parent 2: " << parent_2 << std::endl;
		assert(parent_1 && parent_2);
*/
}

Chromosome* mut_decider(double rand, Chromosome* chromosome_ptr, double mut_rate){ std::cout << std::endl << "####mut_deciding" << std::endl;
	if(true /*rand < mut_rate*/){
//		std::cout << "rand: " << rand << ", mut_rate: " << mut_rate << std::endl;
		chromosome_ptr->mutate();
		std::cout << "mutating went ok :)" << std::endl;
	}else{
		std::cout << "not mutating" << std::endl;
	}
	return chromosome_ptr;
}

// Return a copy of the chromosome with the highest fitness.
// ^ This is literally impossible, the chromosome.hh implementation we were given deletes the copy and assignment constructor. I'm going with the moodle instructions, which say to return a pointer to the best chromosome.
const Chromosome* Deme::get_best() const{ std::cout << std::endl << "####getting best" << std::endl;
	Chromosome* best = pop_[0];
	if(pop_.size() > 1){
		double best_fitness = best->get_fitness();
		std::cout << "pop_ size " << pop_.size() << std::endl;
		for(vec_size_t i = 1; i < pop_.size(); ++i){
			std::cout << "checking " << i << std::endl;
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
Chromosome* Deme::select_parent(){	std::cout << std::endl << "####selecting parent" << std::endl;
	std::vector<vec_size_t> fps_table(pop_.size());	//table for fitness proportionate selection

	std::cout << "fitness: " << pop_[0]->get_fitness();
	fps_table[0] = pop_[0]->get_fitness();	//initialize first value as fitness of first cromosome
	std::cout << ", put into table as " << fps_table[0] << std::endl;
	for(vec_size_t i = 1; i < pop_.size(); ++i){
		std::cout << "fitness: " << pop_[i]->get_fitness();
		fps_table[i] = fps_table[i-1] + pop_[i]->get_fitness();	//each chromosome "takes up" the amount of "space" proportionate to their fitness
		std::cout << ", put into table as " << fps_table[i] << std::endl;
	}

	vec_size_t rand = (generator_() % fps_table[fps_table.size() - 1]); //number between 0 and the sum of all fitnesses, i.e. the value of the last element
	std::cout << "selector: " << rand << std::endl;

	int i = 0;
	while(fps_table[i] < rand){	//iterates through the table until the random value is within the "space" of the chromosome's fitness on the table
//		std::cout << "i: " << i << ", table value: " << fps_table[i] << std::endl;
		++i;
	}
	std::cout << "today's lucky chromosome is " << i << " with fps table value " << fps_table[i] << std::endl;
//	if(i == int(fps_table.size() - 1)){ std::cout << "#####################congratulations bitch! you fixed the problem##########" << std::endl; }
	assert(pop_[i]);
	return pop_[i];
}
