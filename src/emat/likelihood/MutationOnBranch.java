package emat.likelihood;

import java.util.Comparator;

/** represents a mutation for a given branch and site **/
class MutationOnBranch implements Comparable<MutationOnBranch> {
	
	public MutationOnBranch(float brancheFraction, int stateTransition) {
		this.brancheFraction = brancheFraction;
		this.stateTransition = stateTransition;
	}

	/** 
	 * branchFraction is fraction of branch above node with number nodeNr 
	 * that defines the location of a mutation
	 * **/
	float brancheFraction;
	
	/**
	 * stateTransition = fromState * stateCount + toState
	 * **/
	int stateTransition;


	@Override
	public int compareTo(MutationOnBranch o) {
		if (brancheFraction < o.brancheFraction) {
			return -1;
		} else if (brancheFraction > o.brancheFraction) {
			return 1;
		}
		return 0;
	}
}