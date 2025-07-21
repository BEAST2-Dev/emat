package emat.likelihood;

import java.text.DecimalFormat;

/** represents a mutation for a given branch and site **/
public class MutationOnBranch implements Comparable<MutationOnBranch> {
	static String chars = "ACGTXXXXXXXXXXXXXXXX";
	
	public MutationOnBranch(int nodeNr, double brancheFraction, int fromState, int toState, int siteNr) {
		this.nodeNr = nodeNr;
		this.brancheFraction = brancheFraction;
		this.fromState = fromState;
		this.toState = toState;
		this.siteNr = siteNr;
		if (brancheFraction > 1) {
			throw new IllegalArgumentException("branch fraction must be less than 1");
		}
	}

	/** 
	 * branchFraction is fraction of branch above node with number nodeNr 
	 * that defines the location of a mutation
	 * **/
	protected double brancheFraction;
	
	/**
	 * stateTransition = fromState * stateCount + toState
	 * **/
	protected int fromState, toState;
	
	protected int siteNr;
	protected int nodeNr;


	public double getBrancheFraction() {
		return brancheFraction;
	}


	public void setBrancheFraction(double brancheFraction) {
		if (brancheFraction > 1) {
			throw new IllegalArgumentException("branch fraction must be less than 1");
		}
		this.brancheFraction = brancheFraction;
	}


	public int getFromState() {
		return fromState;
	}


	public void setFromState(int fromState) {
		this.fromState = fromState;
	}

	public int getToState() {
		return toState;
	}
	
	public int siteNr() {
		return siteNr;
	}
	
	public double brancheFraction() {
		return brancheFraction;
	}


	public void setToState(int toState) {
		this.toState = toState;
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
		return chars.charAt(fromState) + "" + siteNr + chars.charAt(toState) +"(" + f.format(brancheFraction) +")";
	}


	public void setNodeNr(int nr) {
		this.nodeNr = nr;
	}
}