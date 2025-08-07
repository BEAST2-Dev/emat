package emat.operators;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import emat.likelihood.MutationOnBranch;
import emat.likelihood.MutationState;
import emat.likelihood.MutationStateTreeLikelihood;
import emat.substitutionmodel.EmatSubstitutionModel;

@Description("Gibbs operator that resamples mutations on a node and surrounding branches")
public class MutationOnNodeResampler extends Operator {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree",
			Validate.REQUIRED);
	final public Input<MutationStateTreeLikelihood> likelihoodInput = new Input<>("likelihood",
			"tree likelihood from which substitution model and branch rate models are obtained for stochastic mapping",
			Validate.REQUIRED);

	protected MutationState state;
	protected EmatSubstitutionModel substModel;
	protected BranchRateModel clockModel;
	protected int stateCount;

	@Override
	public void initAndValidate() {
		state = stateInput.get();
		stateCount = state.getStateCount();

		substModel = likelihoodInput.get().substModelInput.get();
		clockModel = likelihoodInput.get().branchRateModelInput.get();

	}

	@Override
	public double proposal() {
		// randomly select internal node
		TreeInterface tree = state.treeInput.get();
		int nodeNr = tree.getLeafNodeCount() + FastRandomiser.nextInt(tree.getInternalNodeCount());
		Node node = tree.getNode(nodeNr);

		if (node.isRoot()) {
			resampleRoot(node, EmatSubstitutionModel.M_MAX_JUMPS);
		} else {
			resample(node, EmatSubstitutionModel.M_MAX_JUMPS);
		}

		return Double.POSITIVE_INFINITY;
	}


	private void resampleRoot(Node root, final int M_MAX_JUMPS) {
		List<MutationOnBranch> branchMutationsLeft = new ArrayList<>();
		List<MutationOnBranch> branchMutationsRight = new ArrayList<>();
		int[] states = state.getNodeSequenceForUpdate(root.getNr());
		
		MutationOperatorUtil.resampleRoot(root, M_MAX_JUMPS,
				branchMutationsLeft,
				branchMutationsRight,
				states,
				state,
				substModel,
				clockModel);
		
		state.setBranchMutations(root.getLeft().getNr(), branchMutationsLeft);
		state.setBranchMutations(root.getRight().getNr(), branchMutationsRight);
	}


	protected void resample(Node node, final int M_MAX_JUMPS) {
		int nodeNr = node.getNr();
		List<MutationOnBranch> branchMutations = new ArrayList<>();
		List<MutationOnBranch> branchMutationsLeft = new ArrayList<>();
		List<MutationOnBranch> branchMutationsRight = new ArrayList<>();
		int [] nodeSequence = state.getNodeSequenceForUpdate(nodeNr);
		
		MutationOperatorUtil.resample(node, M_MAX_JUMPS, 
				branchMutations, 
				branchMutationsLeft, 
				branchMutationsRight,
				nodeSequence,
				state, substModel, clockModel);

		state.setBranchMutations(nodeNr, branchMutations);
		state.setBranchMutations(node.getLeft().getNr(), branchMutationsLeft);
		state.setBranchMutations(node.getRight().getNr(), branchMutationsRight);
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
