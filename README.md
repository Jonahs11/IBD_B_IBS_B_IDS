This is the repo for the Personal Genomics (CSE 284) project of Michael Iter, Dan Musachio, and Jonah Silverman.

We implemented a method to detect Identical By Descent (IBD) segments using unphased genotypes between two individuals.

We are calling our tool IBDbIBSbIDS, as we are using runs of Identical by State (IBS) alleles to infer when the chromosomal segment being considered shares the same ancestral origin (IBD), and we are IDS (Iter, Dan, Silverman).

Our algorithm is inspired by and will benchmark against the tool TRUFFLE[^1].

[^1]: Dimitromanolakis A, Paterson AD, Sun L. Fast and Accurate Shared Segment Detection and Relatedness Estimation in Un-phased Genetic Data via TRUFFLE. Am J Hum Genet. 2019 Jul 3;105(1):78-88. doi: 10.1016/j.ajhg.2019.05.007. Epub 2019 Jun 6. PMID: 31178127; PMCID: PMC6612710.
