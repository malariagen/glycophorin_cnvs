# simulate_uc_events

The program in th


## Compilation

`simulate_uc_events` requires a recent version of gcc that has C++11 support.  Compilation can be acheived with:
```
g++ -std=c++11 simulate_uc_events simulate_uc_events.cpp
```

## Usage

Running `simulate_uc_events` with no arguments gives a brief usage summary.  The general usage is:
```
simulate_uc_events <target profile> <number of generations>
```

`simulate_uc_events` works by attempting to generate a *target* copy number profile using crossover events that occur at predetermined positions (i.e. between predetermined segments), starting from a *reference* chromosome.  The reference chromosome is assumed to contain each segment exactly once; it is always of the form `012345...`.  (Each element in this profile is a single character, and the order is that specified by ASCII.  We recommend using the numeric digits unless more segments are required.)  `simulate_uc_events` works iteratively.  In each iteration, `simulate_uc_events` considers all possible ways to generate new chromosomes from the chromosomes already created.  Then, after a specified number of iterations, `simulate_uc_events` outputs all the chromosomes that match a specified copy number profile.

## Example

For example, suppose we are interested in chromosomes that (relative to a reference sequence) have two duplicated segments separated by a non-duplicated segment.  In the syntax of `simulate_uc_events`, these could be encoded as the target profile:
```
0112334
```

This profile encodes a chromosome with 2 tandem copies of segment 1, 2 of segment 3, and one copy of the segment in between.  `simulate_uc_events` then tries to generate this profile by crossover events, starting with the reference sequence `01234`.  For example, to use `simulate_uc_events` to generate possible histories of the above profile in two iterations:

```
simulate_uc_events "0 1 1 2 3 3 4" 2
```
The output is:
```
# simulate_uc_events
# Starting sequence is: 01234
# Target is: 0112334
# Attempting to make target in 2 generations.
Computing CNVS in generation 1...
After generation 1: 13 haplotypes.
Computing CNVS in generation 2...
After generation 2: 355 haplotypes.
# After 3 generations, a total of 355 CNVs are possible.
#?After 3 generations, a total of 54 histories generate the target copy number profile.
index	result	generation	left	right	break	copy_number_profile_match	exact_match
0	0112334	2	01123|4	0112|34	3|3	1	1
1	0112334	2	01123|4	01212|34	3|3	1	1
2	0112334	2	01123|4	0122|34	3|3	1	1
(etc.)

```

This output tells us a few things

* It's possible to make 355 distinct chromosomes from 2 unequal crossover events given 3 possible breakpoints.
* 54 possible one- or two-event histories generate the target *copy number profile* (i.e. give rise to chromosomes with the same number of each segment as the target profile).
* Of these, 22 histories generate the target profile exactly (as indicated in the `exact_match` column).

The output additionally includes events that appear in possible histories of the matching profiles.

For example, the output tells us that one way to make the target profile `0112334` is via the history

| event | left sequence | right sequence | result       |
| --- | ------------- | -------------- | -------------|
| 1 | `01|234` | `0|1234` | `011234` |
| 2 | `01123|4` | `0112|34` | `0112334` |

