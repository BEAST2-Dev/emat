package emat.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import emat.likelihood.Edit;
import emat.likelihood.EditableNode;
import emat.likelihood.EditableTree;
import emat.likelihood.MutationOnBranch;
import emat.likelihood.MutationState;
import emat.likelihood.MutationStateTreeLikelihood;
import emat.stochasticmapping.StochasticMapping;
import emat.stochasticmapping.TimeStateInterval;
import emat.stochasticmapping.UniformisationStochasticMapping;

@Description("Subtree-prune-regraft base operator for EditableTrees")
public class SPR extends MutationOnNodeResampler {
    final public Input<EditableTree> treeInput = new Input<>("tree", "beast.tree on which this operation is performed", Validate.REQUIRED);

    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
	
	protected EditableTree tree;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		super.initAndValidate();
	}


	@Override
	public double proposal() {
		// select two leaf nodes
		EditableNode n1 = (EditableNode) tree.getNode(Randomizer.nextInt(tree.getLeafNodeCount()));
		EditableNode n2 = n1;
		while (n2 == n1) {
			n2 = (EditableNode) tree.getNode(Randomizer.nextInt(tree.getLeafNodeCount()));
		}
		subtreePruneRegraft(n1, n2, Randomizer.nextDouble() * n2.getLength());
		return 0;
	}
	
	/**
	 * Perform an SPR move
	 * @param subtree = MRCA of the subtree to be removed
	 * @param targetBranch = MRCA above which the subtree will be grafted 
	 * @param newHeight = height at which the subtree will be grafted
	 */
	protected void subtreePruneRegraft(EditableNode subtree, EditableNode targetBranch, double newHeight) {
		
		Node parent = subtree.getParent();
		Node sibling = getOtherChild(parent, subtree);
		
		Edit e = tree.doSPR(subtree.getNr(), targetBranch.getNr(), newHeight);
		
		// resample mutations on branches above and below subtree
		resample(subtree.getParent());
		
		// resample mutations on original branch between parent of subtree and its sibling
		List<MutationOnBranch> branchMutations = new ArrayList<>();
		int siblingNr = sibling.getNr();
		int [] states = state.getNodeSequence(siblingNr);
		int [] parentstates = state.getNodeSequence(parent.getNr());
		double [] weightsN = setUpWeights(sibling.getLength() * clockModel.getRateForBranch(sibling));
		for (int i = 0; i < states.length; i++) {
			double [] p = new double[M_MAX_JUMPS];
			for (int r = 0; r < M_MAX_JUMPS; r++) {
				p[r] = weightsN[r] * qUnifPowers.get(r)[states[i]][parentstates[i]];
			}			
			int N = Randomizer.randomChoicePDF(p);
			generatePath(siblingNr, i, parentstates[i], states[i], N, branchMutations);
		}
		state.setBranchMutations(siblingNr, branchMutations);

		
		
//		StochasticMapping mapping = new UniformisationStochasticMapping();
//		mapping.setRatematrix(substModel.getRateMatrix());
//
//		// sibling gets mutations from parent branch
//		List<MutationOnBranch> siblingMutations = state.getMutationList(e.siblingNr());
//		List<MutationOnBranch> parentMutations = state.getMutationList(e.parentNr());
//		List<MutationOnBranch> newSiblingMutations = new ArrayList<>();
//		newSiblingMutations.addAll(siblingMutations);
//		newSiblingMutations.addAll(parentMutations);
//		state.setBranchMutations(e.siblingNr(), newSiblingMutations);
//
//		// target mutations get split at newHeight:
//		List<MutationOnBranch> targetMutations = state.getMutationList(e.targetNr());
//		List<MutationOnBranch> newTargetMutations = new ArrayList<>();
//		List<MutationOnBranch> newParentMutations = new ArrayList<>();
//		
//		// states will contain states for parent node after evolving through mutations from target node
//		// this will be used as endpoint for simulating mutations going from node to parent
//		int [] states = state.getNodeSequence(e.targetNr());
//		double parentBranchFraction = (newHeight - tree.getNode(e.targetNr()).getHeight()) / tree.getNode(e.targetNr()).getLength();
//		int i = 0;
//		while (i < targetMutations.size() && targetMutations.get(i).brancheFraction() < parentBranchFraction) {
//			MutationOnBranch m = targetMutations.get(i);
//			newTargetMutations.add(m);
//			states[m.siteNr()] = m.getFromState();
//			i++;
//		}
//		while (i < targetMutations.size()) {
//			newParentMutations.add(targetMutations.get(i++));
//		}
//		state.setBranchMutations(e.targetNr(), newTargetMutations);
//		state.setBranchMutations(e.parentNr(), newParentMutations);
//		
//		// sample new path for node
//		Node node = tree.getNode(e.nodeNr());
//		double length = node.getLength() * clockModel.getRateForBranch(node);
//		List<TimeStateInterval> path = mapping.generatePath(state.getNodeSequence(e.nodeNr()), states, length);
//		List<MutationOnBranch> nodeMutations = new ArrayList<>(); 
//		for (int j = 0; j < path.size() - 1; j++) {
//			TimeStateInterval interval = path.get(j);
//			if (path.get(i).endTime() < length) {
//				nodeMutations.add(new MutationOnBranch(e.nodeNr(), interval.endTime()/length, interval.state(), path.get(j+1).state(), interval.site()));
//			}
//		}
//		state.setBranchMutations(e.nodeNr(), nodeMutations);
		
	}
	
    /**
     * @param parent the parent
     * @param child  the child that you want the sister of
     * @return the other child of the given parent.
     */
    protected Node getOtherChild(final Node parent, final Node child) {
        if (parent.getLeft().getNr() == child.getNr()) {
            return parent.getRight();
        } else {
            return parent.getLeft();
        }
    }

	
}
