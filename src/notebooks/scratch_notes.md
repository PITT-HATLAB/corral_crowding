Outline of paper:

0. Motivation, hardware crowding + connectivity (SWAPs, high-fidelity 2Q gates, etc) -> transition into (1) SNAIL-coupling hardware architecture co-design. this paper is spiritual successor to corral paper, so we can talk about the short-comings of that paper here, ie how do we know what is the right size corral to use (previously we guesstimated)? More- in the first paper we showed that dense connectivity is better at reducing SWAPs (let's try again with better SWAP routing algorithms) but more importantly, we didn't account for the fact that we can't keep gate fidelities perfectly ideal. Now we have better analysis to make sure we are making appropriate tradeoff considerations)
1. Derive physics-informed frequency constraints
   - expand terms, variation of RWA (keep some of the slow non-zero terms)
   - stark shift, kerr shift, coupling theory, etc
2. Frequency allocation problem (classical CS + regularly appears in engineering real-world network e.g.), constraint satisfaction, graph coloring (?), linear programming, etc
3. spectator errors as coherent errors (errors that propagate outside the module can be cancelled out by the other module's SNAIL with a compensation pulse :))
4. bring full circle to topology/connectivity problem
   - feasibile architecture designs with reasonable fab requirements that maintain high fidelity 2Q gates
   - how can we push to increasingly high connected systems (or do we even want to?)
   - can look at these thresholds by incorporting incoherent lifetime errors, now with tradeoffs between slower gates (ie DRAG style compensation or more plainly by narrowing frequency spectrum)

________

Assumptions, limitations of our model:
- that we flux tune snail to have no fourth-order nonlinearity
- pulse shapes have non-zero width spectral components
- ideal transmon nonlineality
- fab precison of target couplings,frequencies
- we are keeping some slow-rotating terms but we are not keeping all terms
- *for now*, I am assuming stark-shift is small, because I am having trouble defining what I expect an updated definition of $\omega_p$ should be, if I know pump will shift around qubit frequencies. I need further discussion on this point. But for our initial model of intrinsic error budgeting we can assume we don't correct for this.


for the nearest-neighbor model:
- assume that couplings between snails is small, this helps a ton because if each snail only see its own modules qubit then we have more straightforward hybridization model
- again for simplicity we assume that all qubit-snail couplings are identical, lambda_i = lambda
- BIG if multiple snails are pumped at same time then we have multiple pump displacments and things get confusing, for our initial model we are going to assume the neighboring snail drives are small and are only acting as small coherent error correction pulses; but for a practical architectue we would want to verify that constraints allow for simulatenous gates
