# hw-10

Rebecca: Chromosome class
==============

The compiler did not agree with the generator_ member that was included in the skeleton code Eitan provided, so I decided to comment it out and I reused the random int generator code that I implemented in random_permutation in HW9. So the Chromosome constructor only establishes order_ and cities_ptr_. Since the Chromosome does not have any members that allocate new memory, the destructor function just needed to reassert that the Chromosome was valid. Then for mutate(), I defined two positions in order_ using randomly-generated ints to indicate their places, and simply swapped the two values. recombine() does a similar technique, although it is cumbersome, it basically just generates random positions for each call of create_crossover_child, and passes This and the pointer to the other Chromosome as the first and second parameters. get_fitness() returns 1 divided by the total path distance of the current ordering, since the smaller this distance is, the larger the double produced. is_valid takes advantage of the fact that order_ must be a permutation of the range from 0 to the size of the Cities, so it just sorts a copy of order_ and compares that to the range, which is calculated using a helper function. Finally, is_in_range checks every value of order_ between begin and end, returning True if it finds a match with its target value, and False otherwise.

El: Deme class
==============

I could've manually managed memory in the `deme` constructor/destructor, but valgrind never likes that and so I decided to just start by changing the `pop_` member to be a vector of shared pointers rather than pointers. If this comes back to bite me later, oof. I don't think it should though. Also this might be helpful if during one of the methods i delete members from the population prematurely, though idk.

