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

	private double[] weightsN, leftWeightsN, rightWeightsN;

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
		Node left = root.getLeft();
		Node right = root.getRight();
		int[] leftStates = state.getNodeSequence(left.getNr());
		int[] rightStates = state.getNodeSequence(right.getNr());
		
		double totalTimeLeft = left.getLength() * clockModel.getRateForBranch(left);
		double totalTimeRight = right.getLength() * clockModel.getRateForBranch(right);

		double totalTime = totalTimeLeft + totalTimeRight;
		weightsN = MutationOperatorUtil.setUpWeights(totalTime, substModel.getLambdaMax(), M_MAX_JUMPS);
		
		List<MutationOnBranch> branchMutations = new ArrayList<>();

		List<double[][]> qUnifPowers = substModel.getQUnifPowers();
		for (int i = 0; i < leftStates.length; i++) {
			double [] p = new double[M_MAX_JUMPS];
			for (int r = 0; r < M_MAX_JUMPS; r++) {
				p[r] = weightsN[r] * qUnifPowers.get(r)[leftStates[i]][rightStates[i]];
			}			
			int N = FastRandomiser.randomChoicePDF(p);
			MutationOperatorUtil.generatePath(left.getNr(), i, rightStates[i], leftStates[i], N, qUnifPowers ,branchMutations);
		}

		Collections.sort(branchMutations);
		
		int[] states = state.getNodeSequenceForUpdate(root.getNr());
		System.arraycopy(leftStates, 0, states, 0, states.length);
		
		double leftFraction = totalTimeLeft / totalTime;
		double rightFraction = 1 - leftFraction;
		List<MutationOnBranch> branchMutationsLeft = new ArrayList<>();
		List<MutationOnBranch> branchMutationsRight = new ArrayList<>();
		for (MutationOnBranch m : branchMutations) {
			if (m.brancheFraction() < leftFraction) {
				m.setBrancheFraction(m.getBrancheFraction() / leftFraction);
				branchMutationsLeft.add(m);
				states[m.siteNr()] = m.getToState();
			} else {
				branchMutationsRight.add(
						new MutationOnBranch(right.getNr(), (1-m.getBrancheFraction()) / rightFraction, m.getToState(), m.getFromState(), m.siteNr()));
			}
		}
		
		
		state.setBranchMutations(left.getNr(), branchMutationsLeft);
		state.setBranchMutations(right.getNr(), branchMutationsRight);
		
		
	}


	protected void resample(Node node, final int M_MAX_JUMPS) {
		int nodeNr = node.getNr();

		int[] states = state.getNodeSequence(node.getParent().getNr());
		int[] leftStates = state.getNodeSequence(node.getLeft().getNr());
		int[] rightStates = state.getNodeSequence(node.getRight().getNr());
		int [] nodeSequence = state.getNodeSequenceForUpdate(nodeNr);

		double totalTime = node.getLength() * clockModel.getRateForBranch(node);
		double totalTimeLeft = node.getLeft().getLength() * clockModel.getRateForBranch(node.getLeft());
		double totalTimeRight = node.getRight().getLength() * clockModel.getRateForBranch(node.getRight());

		// TODO: cach weights per branch & update only when evolutionary distance = (branch lengths * clock rate) changes
		double lambdaMax = substModel.getLambdaMax();
		List<double[][]> qUnifPowers = substModel.getQUnifPowers();
		weightsN = MutationOperatorUtil.setUpWeights(totalTime, lambdaMax, M_MAX_JUMPS);
		leftWeightsN = MutationOperatorUtil.setUpWeights(totalTimeLeft, lambdaMax, M_MAX_JUMPS);
		rightWeightsN = MutationOperatorUtil.setUpWeights(totalTimeRight, lambdaMax, M_MAX_JUMPS);
		
		int Relevant_M_MAX_JUMPS = 3;
		
		List<MutationOnBranch> branchMutations = new ArrayList<>();
		List<MutationOnBranch> branchMutationsLeft = new ArrayList<>();
		List<MutationOnBranch> branchMutationsRight = new ArrayList<>();
		
		double [] pNodeState = new double[stateCount];

		for (int i = 0; i < states.length; i++) {
			for (int nodeState = 0; nodeState < stateCount; nodeState++) {
				double p = 0;
				for (int r = 0; r < Relevant_M_MAX_JUMPS; r++) {
					double p0 = weightsN[r] * qUnifPowers.get(r)[states[i]][nodeState];
					for (int s = 0; s < Relevant_M_MAX_JUMPS; s++) {
						double p1 = leftWeightsN[s] * qUnifPowers.get(s)[nodeState][leftStates[i]] * p0;
						for (int t = 0; t < Relevant_M_MAX_JUMPS; t++) {
							p +=  p1 * rightWeightsN[t] * qUnifPowers.get(t)[nodeState][rightStates[i]];
						}
					}
				}
				pNodeState[nodeState] = p;
			}
			int nodeState = FastRandomiser.randomChoicePDF(pNodeState);
			nodeSequence[i] = nodeState;
			
			double [] p = new double[M_MAX_JUMPS];
			for (int r = 0; r < M_MAX_JUMPS; r++) {
				p[r] = weightsN[r] * qUnifPowers.get(r)[states[i]][nodeState];
			}			
			int N = FastRandomiser.randomChoicePDF(p);
			MutationOperatorUtil.generatePath(nodeNr, i, states[i], nodeState, N, qUnifPowers, branchMutations);
			
			for (int r = 0; r < M_MAX_JUMPS; r++) {
				p[r] = leftWeightsN[r] * qUnifPowers.get(r)[nodeState][leftStates[i]];
			}			
			int Nleft = FastRandomiser.randomChoicePDF(p);
			MutationOperatorUtil.generatePath(node.getLeft().getNr(), i, nodeState, leftStates[i], Nleft, qUnifPowers, branchMutationsLeft);

			for (int r = 0; r < M_MAX_JUMPS; r++) {
				p[r] = rightWeightsN[r] * qUnifPowers.get(r)[nodeState][rightStates[i]];
			}			
			int Nright = FastRandomiser.randomChoicePDF(p);
			MutationOperatorUtil.generatePath(node.getRight().getNr(), i, nodeState, rightStates[i], Nright, qUnifPowers, branchMutationsRight);
		}

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
