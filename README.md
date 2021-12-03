# hw-10

Rebecca: Chromosome class
==============

The compiler did not agree with the generator_ member that was included in the skeleton code Eitan provided, so I decided to comment it out and I reused the random int generator code that I implemented in random_permutation in HW9. So the Chromosome constructor only establishes order_ and cities_ptr_. Since the Chromosome does not have any members that allocate new memory, the destructor function just needed to reassert that the Chromosome was valid. Then for mutate(), I defined two positions in order_ using randomly-generated ints to indicate their places, and simply swapped the two values. recombine() does a similar technique, although it is cumbersome, it basically just generates random positions for each call of create_crossover_child, and passes This and the pointer to the other Chromosome as the first and second parameters. get_fitness() returns 1000 divided by the total path distance of the current ordering, since the smaller this distance is, the larger the double produced. is_valid takes advantage of the fact that order_ must be a permutation of the range from 0 to the size of the Cities, so it just sorts a copy of order_ and compares that to the range, which is calculated using a helper function. Finally, is_in_range checks every value of order_ between begin and end, returning True if it finds a match with its target value, and False otherwise.

For part 3, a population size of 50 and mutation rate of 0.75 seemed to be a good compromise, generally producing a result between 18-19 in about 35 seconds.
Note: modifying tsp.cc for this part was a lot more confusing as we would have had to find ways for our original code to agree with what Eitan had written. So, I just made randomized.tsv using our code from HW9, since I knew that would work, and then moved it into the HW10 folder to make the comparison gif.

El: Deme class
==============

# deme.hh

unchanged, other than formatting

# deme.cc

### frac
Generates a random decimal number between 0 and 1 by dividing a random number from the generator by the generator's maximum value

### constructor

`Deme` populates `pop_` with new chromosomes it allocates memory for and stores the mutation rate it's given into `mut_rate_`.

It also seeds the random number generator with a number bansed on the amount of time since the epoch, something Rebecca found for the last homework and which I thought would be more "random" than c's `rand()`.

### destructor

This just goes through `pop_` and deallocates the memory for each chromosome. Valgrind shows no memory leaks, so it looks like this works well.

### compute_next_generation()

This creates a vector of pairs chosen from `pop_` using `select_parent()` (see below). It may be marginally faster to initialize the vector as size `pop_size() / 2`, but I didn't have time to implement that.

Then for each pair in the vector, it uses `mut_decider` to decider whether to mutate, and mutate if decided, each chromosome in the pair. It then recombines the pair and stores that into a new pair, the chromosomes of which are stored into `pop_`. It then goes through the old chromosomes of `pop_`, deallocates their memory, and removes them from the vector.

### mut_decider

This decides to mutate a chromosome it is passed by poitner based on if a random number generated by `frac()` (see above) is lower than the mutation rate is passed, and mutates it if it decides to.

### get_best()

This iterates through the `pop_` and if the fitness of the chromosome it's checking is greater than the previous chromosome with best fitness, stores it as the new best fitness. It then returns the chromosome with the greatest fitness.

### select_parent()

This iterates through all of `pop_` and adds all of the fitnesses in it to store in `fit_sum`. A random number between 0 and a maximum value, which is what is stored in `fit_sum`, is generated and stored in `selector`. It then goes through `pop_` again with an interator `i`, and adds the fitnesses once again until the sum of fitnesses is greater than `selector`. Once it is greater, this means that `selector` falls within the fitness "range" of the next chromosome on `pop_`, so `i` is increased once more and the function returns `pop_[i]`.
