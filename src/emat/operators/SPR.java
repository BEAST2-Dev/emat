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
public class SPR extends EditableTreeOperator {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree", Validate.REQUIRED);
	final public Input<MutationStateTreeLikelihood> likelihoodInput = new Input<>("likelihood" ,"tree likelihood from which substitution model and branch rate models are obtained for stochastic mapping", Validate.REQUIRED);

    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
	
	protected EditableTree tree;
	protected MutationState state;
	protected BranchRateModel.Base clockModel;
	private GeneralSubstitutionModel substModel;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		state = stateInput.get();
		substModel = (GeneralSubstitutionModel) likelihoodInput.get().getSubstModel();
		clockModel = likelihoodInput.get().branchRateModelInput.get();
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
		
		Edit e = tree.doSPR(subtree.getNr(), targetBranch.getNr(), newHeight);
		
		StochasticMapping mapping = new UniformisationStochasticMapping();
		mapping.setRatematrix(substModel.getRateMatrix());

		// sibling gets mutations from parent branch
		List<MutationOnBranch> siblingMutations = state.getMutationList(e.siblingNr());
		List<MutationOnBranch> parentMutations = state.getMutationList(e.parentNr());
		List<MutationOnBranch> newSiblingMutations = new ArrayList<>();
		newSiblingMutations.addAll(siblingMutations);
		newSiblingMutations.addAll(parentMutations);
		state.setBranchMutations(e.siblingNr(), newSiblingMutations);

		// target mutations get split at newHeight:
		List<MutationOnBranch> targetMutations = state.getMutationList(e.targetNr());
		List<MutationOnBranch> newTargetMutations = new ArrayList<>();
		List<MutationOnBranch> newParentMutations = new ArrayList<>();
		
		// states will contain states for parent node after evolving through mutations from target node
		// this will be used as endpoint for simulating mutations going from node to parent
		int [] states = state.getNodeSequence(e.targetNr());
		double parentBranchFraction = (newHeight - tree.getNode(e.targetNr()).getHeight()) / tree.getNode(e.targetNr()).getLength();
		int i = 0;
		while (i < targetMutations.size() && targetMutations.get(i).brancheFraction() < parentBranchFraction) {
			MutationOnBranch m = targetMutations.get(i);
			newTargetMutations.add(m);
			states[m.siteNr()] = m.getFromState();
			i++;
		}
		while (i < targetMutations.size()) {
			newParentMutations.add(targetMutations.get(i++));
		}
		state.setBranchMutations(e.targetNr(), newTargetMutations);
		state.setBranchMutations(e.parentNr(), newParentMutations);
		
		// sample new path for node
		Node node = tree.getNode(e.nodeNr());
		double length = node.getLength() * clockModel.getRateForBranch(node);
		List<TimeStateInterval> path = mapping.generatePath(state.getNodeSequence(e.nodeNr()), states, length);
		List<MutationOnBranch> nodeMutations = new ArrayList<>(); 
		for (int j = 0; j < path.size() - 1; j++) {
			TimeStateInterval interval = path.get(j);
			if (path.get(i).endTime() < length) {
				nodeMutations.add(new MutationOnBranch(e.nodeNr(), interval.endTime()/length, interval.state(), path.get(j+1).state(), interval.site()));
			}
		}
		state.setBranchMutations(e.nodeNr(), nodeMutations);
	}
	
}
