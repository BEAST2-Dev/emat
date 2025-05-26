package emat.likelihood;

import java.text.DecimalFormat;

/** represents a mutation for a given branch and site **/
public class MutationOnBranch implements Comparable<MutationOnBranch> {
	
	public MutationOnBranch(int nodeNr, double brancheFraction, int stateTransition, int siteNr) {
		this.nodeNr = nodeNr;
		this.brancheFraction = brancheFraction;
		this.stateTransition = stateTransition;
		this.siteNr = siteNr;
	}

	/** 
	 * branchFraction is fraction of branch above node with number nodeNr 
	 * that defines the location of a mutation
	 * **/
	protected double brancheFraction;
	
	/**
	 * stateTransition = fromState * stateCount + toState
	 * **/
	protected int stateTransition;
	
	protected int siteNr;
	protected int nodeNr;


	public double getBrancheFraction() {
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
	
	
	private static DecimalFormat f = new DecimalFormat("#.###");
	
	@Override
	public String toString() {
		return nodeNr + "@" + siteNr + ":" + stateTransition/4 + "=>" + stateTransition%4 +"(" + f.format(brancheFraction) +")";
	}
}