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
	
	final static private boolean debug = true;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		super.initAndValidate();
	}


	@Override
	public double proposal() {
		double logHR = 0;
		
		// select two leaf nodes
		EditableNode n1 = (EditableNode) tree.getNode(FastRandomiser.nextInt(tree.getLeafNodeCount()));
		while (n1.getParent().isRoot()) {
			n1 = (EditableNode) tree.getNode(FastRandomiser.nextInt(tree.getLeafNodeCount()));
		}
		
		Node parent = n1.getParent();
		double height = parent.getHeight();
		
		List<EditableNode> eligableNodes = new ArrayList<>();
		for (Node n : tree.getNodesAsArray()) {
			if (!n.isRoot() && n.getHeight() < height && n.getParent().getHeight() > height) {
				if (n != n1) {
					eligableNodes.add((EditableNode) n);
				}
			}
		}
		
		EditableNode n2 = eligableNodes.get(FastRandomiser.nextInt(eligableNodes.size()));
		
//		EditableNode n2 = n1;
//		while (n2 == n1 || n2 == getOtherChild(n1.getParent(), n1)) {
//			n2 = (EditableNode) tree.getNode(FastRandomiser.nextInt(tree.getLeafNodeCount()));
//		}
//		height = n2.getHeight() + FastRandomiser.nextDouble() * n2.getLength();
//		subtreePruneRegraft(n1, n2, n2.getHeight() + FastRandomiser.nextDouble() * n2.getLength());

		logHR += subtreePruneRegraft(n1, n2, height);
		return logHR;
	}
	
	/**
	 * Perform an SPR move
	 * @param subtree = MRCA of the subtree to be removed
	 * @param targetBranch = MRCA above which the subtree will be grafted 
	 * @param newHeight = height at which the subtree will be grafted
	 */
	protected double subtreePruneRegraft(EditableNode subtree, EditableNode targetBranch, double newHeight) {
		
		Node parent = subtree.getParent();
		Node sibling = getOtherChild(parent, subtree);
		int siblingNr = sibling.getNr();
		int parentNr = parent.getNr();
		
		
		double logHR = 0;
        int mutationCount = state.getMutationList(subtree.getNr()).size();
        logHR = -mutationCount * Math.log(newHeight - subtree.getHeight());

		
		
		// amalgamate mutations above sibling
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

        logHR += -siblingMutations.size() * Math.log(f);
        logHR += -parentMutations.size() * Math.log(1-f);
		
		
		
		// split target branch mutations
		int targetNr = targetBranch.getNr();
		List<MutationOnBranch> targetMutations = state.getMutationList(targetNr);
		double targetLength = targetBranch.getLength();
		List<MutationOnBranch> newTargetMutations = new ArrayList<>();
		List<MutationOnBranch> newParentMutations = new ArrayList<>();
		f = targetLength / (newHeight - targetBranch.getHeight());
		
		double f2 = (targetBranch.getParent().getHeight() - targetBranch.getHeight()) / (targetBranch.getParent().getHeight() - newHeight); 
		double threshold = 1 / f;
		int [] states0 = state.getNodeSequence(targetNr);
		int [] states = state.getNodeSequenceForUpdate(parentNr);
		int [] parentNodeStates = states.clone();
		System.arraycopy(states0, 0, states, 0, states.length);
		int [] statesOrig = null;
		if (debug) {
			statesOrig = new int[states.length];
			System.arraycopy(states0, 0, statesOrig, 0, states.length);
		}
		for (MutationOnBranch m : targetMutations) {
			if (m.getBrancheFraction() < threshold) {
				newTargetMutations.add(new MutationOnBranch(targetNr, f * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
				states[m.siteNr()] = m.getToState();
			} else {
				newParentMutations.add(new MutationOnBranch(parentNr, f2 * (m.brancheFraction() - threshold), m.getFromState(), m.getToState(), m.siteNr()));
			}
			if (debug) {
				statesOrig[m.siteNr()] = m.getToState();
			}
		}

        logHR += -newTargetMutations.size() * Math.log(f);
        logHR += -newParentMutations.size() * Math.log(f2);
		
		if (debug) {
			int [] gpstates = state.getNodeSequenceForUpdate(targetBranch.getParent().getNr());
			
			for (int i = 0; i < states.length; i++) {
				if (statesOrig[i] != gpstates[i]) {
					System.err.println("Something wrong with the state reconstruction at site "  + i);
				}
			}
		}		
		
		// resample mutations on branch above subtree (only if necessary)
		int nodeNr = subtree.getNr();
		int [] nodeStates = state.getNodeSequence(nodeNr);

		// TODO: only do this when substModel changes
		substModel.setupRateMatrix();
		setRatematrix(substModel.getRateMatrix());

		List<MutationOnBranch> nodeMutations = new ArrayList<>();
		double totalTime = (newHeight - subtree.getHeight()) * clockModel.getRateForBranch(subtree);
		double [] weightsN = setUpWeights(totalTime);
		double [] p = new double[M_MAX_JUMPS];

//		for (int i = 0; i < states.length; i++) {
//			int nodeState = nodeStates[i];
//			for (int r = 0; r < M_MAX_JUMPS; r++) {
//				p[r] = weightsN[r] * qUnifPowers.get(r)[states[i]][nodeState];
//			}			
//			int N = FastRandomiser.randomChoicePDF(p);
//			generatePath(nodeNr, i, states[i], nodeState, N, nodeMutations);
//		}

		boolean [] needsResampling = new boolean[states.length];
		for (int i = 0; i < states.length; i++) {
			if (parentNodeStates[i] != states[i]) {
				needsResampling[i] = true;
			}
		}
		// keep mutations on sites that do not differ
		List<MutationOnBranch> currentNodeMutations = state.getMutationList(nodeNr);
		for (MutationOnBranch m : currentNodeMutations) {
			if (!needsResampling[m.siteNr()]) {
				nodeMutations.add(m);
			}
		}
		// resample sites that do differ
		for (int i = 0; i < states.length; i++) {
			if (needsResampling[i]) {
				int nodeState = nodeStates[i]; 
				for (int r = 0; r < M_MAX_JUMPS; r++) {
					p[r] = weightsN[r] * qUnifPowers.get(r)[states[i]][nodeState];
				}			
				int N = FastRandomiser.randomChoicePDF(p);
				generatePath(nodeNr, i, states[i], nodeState, N, nodeMutations);
			}
		}

		
		Edit e = tree.doSPR(subtree.getNr(), targetBranch.getNr(), newHeight);
		
		state.setBranchMutations(siblingNr, newSiblingMutations);
		state.setBranchMutations(nodeNr, nodeMutations);
		state.setBranchMutations(targetNr, newTargetMutations);
		state.setBranchMutations(parentNr, newParentMutations);

		
        logHR += nodeMutations.size() * Math.log(newHeight - subtree.getHeight());
        
        //return 0;
        return logHR;
		
		
		// resample mutations on branches above and below subtree
//		substModel.setupRateMatrix();
//		setRatematrix(substModel.getRateMatrix());
//		resample(subtree.getParent());
		
		// resample mutations on original branch between parent of subtree and its sibling
//		List<MutationOnBranch> branchMutations = new ArrayList<>();
//		int siblingNr = sibling.getNr();
//		int [] states = state.getNodeSequence(siblingNr);
//		int [] parentstates = state.getNodeSequence(sibling.getParent().getNr());
//		double [] weightsN = setUpWeights(sibling.getLength() * clockModel.getRateForBranch(sibling));
//		for (int i = 0; i < states.length; i++) {
//			double [] p = new double[M_MAX_JUMPS];
//			for (int r = 0; r < M_MAX_JUMPS; r++) {
//				p[r] = weightsN[r] * qUnifPowers.get(r)[parentstates[i]][states[i]];
//			}			
//			int N = FastRandomiser.randomChoicePDF(p);
//			generatePath(siblingNr, i, parentstates[i], states[i], N, branchMutations);
//		}
//		state.setBranchMutations(siblingNr, branchMutations);
		
		
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
