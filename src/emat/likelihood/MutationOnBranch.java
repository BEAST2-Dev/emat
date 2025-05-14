package emat.likelihood;

/** represents a mutation for a given branch and site **/
public class MutationOnBranch implements Comparable<MutationOnBranch> {
	
	public MutationOnBranch(float brancheFraction, int stateTransition) {
		this.brancheFraction = brancheFraction;
		this.stateTransition = stateTransition;
	}

	/** 
	 * branchFraction is fraction of branch above node with number nodeNr 
	 * that defines the location of a mutation
	 * **/
	protected float brancheFraction;
	
	/**
	 * stateTransition = fromState * stateCount + toState
	 * **/
	protected int stateTransition;


	public float getBrancheFraction() {
		return brancheFraction;
	}


	public void setBrancheFraction(float brancheFraction) {
		this.brancheFraction = brancheFraction;
	}


	public int getStateTransition() {
		return stateTransition;
	}


	public void setStateTransition(int stateTransition) {
		this.stateTransition = stateTransition;
	}


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