package emat.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.evolution.tree.Node;
import emat.likelihood.MutationOnBranch;
import emat.likelihood.MutationState;

public class MutationOperatorUtil {

	
	/** create mutation list that contains mutations above a node as well as the ones above its parent **/
	public static List<MutationOnBranch> amalgamateBranchWithParentBranch(MutationState state, Node sibling) {
		int siblingNr = sibling.getNr();
		Node parent = sibling.getParent();
		List<MutationOnBranch> siblingMutations = state.getMutationList(siblingNr);
		double siblingLength = sibling.getLength();
		List<MutationOnBranch> parentMutations = state.getMutationList(parent.getNr());
		double parentLength = parent.getLength();
		List<MutationOnBranch> newSiblingMutations = new ArrayList<>();
		double f = siblingLength/(siblingLength + parentLength);
		for (MutationOnBranch m : siblingMutations) {
			newSiblingMutations.add(new MutationOnBranch(siblingNr, f * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
		}
		for (MutationOnBranch m : parentMutations) {
			newSiblingMutations.add(new MutationOnBranch(siblingNr, f + (1-f) * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
		}
		return newSiblingMutations;
	}
}
