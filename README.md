# hw-10

Rebecca: Chromosome class
==============

:(

El: Deme class
==============

I could've manually managed memory in the `deme` constructor/destructor, but valgrind never likes that and so I decided to just start by changing the `pop_` member to be a vector of shared pointers rather than pointers. If this comes back to bite me later, oof. I don't think it should though. Also this might be helpful if during one of the methods i delete members from the population prematurely, though idk.

