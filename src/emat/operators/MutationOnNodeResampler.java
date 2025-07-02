package emat.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.CombinatoricsUtils;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import emat.likelihood.MutationOnBranch;
import emat.likelihood.MutationState;
import emat.likelihood.MutationStateTreeLikelihood;
import emat.stochasticmapping.UniformisationStochasticMapping;

@Description("Gibbs operator that resamples mutations on a node and surrounding branches")
public class MutationOnNodeResampler extends Operator {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree",
			Validate.REQUIRED);
	final public Input<MutationStateTreeLikelihood> likelihoodInput = new Input<>("likelihood",
			"tree likelihood from which substitution model and branch rate models are obtained for stochastic mapping",
			Validate.REQUIRED);

	private MutationState state;
	private GeneralSubstitutionModel substModel;
	private BranchRateModel clockModel;
	private int stateCount;

	private int M_MAX_JUMPS = 20;
	
	private double lambdaMax;
	private double[][] qUnif;
	private List<double[][]> qUnifPowers;

	private double[] weightsN, leftWeightsN, rightWeightsN;

	@Override
	public void initAndValidate() {
		state = stateInput.get();
		stateCount = state.getStateCount();

		substModel = (GeneralSubstitutionModel) likelihoodInput.get().getSubstModel();
		clockModel = likelihoodInput.get().branchRateModelInput.get();

	}

	@Override
	public double proposal() {
		TreeInterface tree = state.treeInput.get();
		setRatematrix(substModel.getRateMatrix());

		// randomly select internal node
		int nodeNr = tree.getLeafNodeCount() + Randomizer.nextInt(tree.getInternalNodeCount()-1);
		Node node = tree.getNode(nodeNr);
		int[] states = state.getNodeSequence(node.getParent().getNr());
		int[] leftStates = state.getNodeSequence(node.getLeft().getNr());
		int[] rightStates = state.getNodeSequence(node.getRight().getNr());

		state.replaceNodeStates(nodeNr);
		int [] nodeSequence = state.getNodeSequence(nodeNr);

		double totalTime = node.getLength() * clockModel.getRateForBranch(node);
		double totalTimeLeft = node.getLeft().getLength() * clockModel.getRateForBranch(node.getLeft());
		double totalTimeRight = node.getRight().getLength() * clockModel.getRateForBranch(node.getRight());

		// TODO: cach weights per branch & update only when evolutionary distance = (branch lengths * clock rate) changes
		weightsN = setUpWeights(totalTime);
		leftWeightsN = setUpWeights(totalTimeLeft);
		rightWeightsN = setUpWeights(totalTimeRight);
		
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
			int nodeState = Randomizer.randomChoicePDF(pNodeState);
			nodeSequence[i] = nodeState;
			
			double [] p = new double[M_MAX_JUMPS];
			for (int r = 0; r < M_MAX_JUMPS; r++) {
				p[r] = weightsN[r] * qUnifPowers.get(r)[states[i]][nodeState];
			}			
			int N = Randomizer.randomChoicePDF(p);
			generatePath(nodeNr, i, states[i], nodeState, N, branchMutations);
			
			for (int r = 0; r < M_MAX_JUMPS; r++) {
				p[r] = leftWeightsN[r] * qUnifPowers.get(r)[nodeState][leftStates[i]];
			}			
			int Nleft = Randomizer.randomChoicePDF(p);
			generatePath(node.getLeft().getNr(), i, nodeState, leftStates[i], Nleft, branchMutationsLeft);

			for (int r = 0; r < M_MAX_JUMPS; r++) {
				p[r] = rightWeightsN[r] * qUnifPowers.get(r)[nodeState][rightStates[i]];
			}			
			int Nright = Randomizer.randomChoicePDF(p);
			generatePath(node.getRight().getNr(), i, nodeState, rightStates[i], Nright, branchMutationsRight);
		}

		state.setBranchMutations(nodeNr, branchMutations);
		state.setBranchMutations(node.getLeft().getNr(), branchMutationsLeft);
		state.setBranchMutations(node.getRight().getNr(), branchMutationsRight);

		return Double.POSITIVE_INFINITY;
	}


	private void generatePath(int nodeNr, int siteNr, int endState, int startState, int N,
			List<MutationOnBranch> mutations) {
        // --- Step 2: Sample the Jump Times (tau_1, ..., tau_N) ---
        double[] jumpTimes = new double[N];
        for (int i = 0; i < N; i++) {
            jumpTimes[i] = Randomizer.nextDouble();
        }
        Arrays.sort(jumpTimes);
        
        // --- Step 3: Sample the Sequence of States (S_0, ..., S_N) ---
        int[] stateSequence = new int[N + 1];
        stateSequence[0] = startState;
        stateSequence[N] = endState;

        for (int k = 1; k < N; k++) { // For S_k (state AFTER k-th jump)
            int prevState = stateSequence[k - 1];
            double[] transitionProbsToNextCandidates = new double[stateCount];
            double sumProbs = 0.0;

            int remainingJumps = N - k;
            double[][] qUnifPowRemaining = qUnifPowers.get(remainingJumps); // Use precomputed

            // TODO: precalculate transitionProbsToNextCandidates[remainingJumps][prev][end][states] 
            for (int nextCandidateState = 0; nextCandidateState < stateCount; nextCandidateState++) {
                double prob = qUnif[prevState][nextCandidateState] * qUnifPowRemaining[nextCandidateState][endState];
                transitionProbsToNextCandidates[nextCandidateState] = prob;
                sumProbs += prob;
            }

            if (sumProbs <= 1e-100) { // Should not happen if P_ij(t) > 0 and N is sampled correctly
                System.err.println("Error: Stuck during state sequence sampling at k=" + k +". Sum of probabilities is zero. This indicates an issue.");
                // This might happen if M_MAX_JUMPS was too small and N was sampled near the edge
                // or if numerical precision led to an inconsistent state.
                // One recovery could be to force to endState if k=N, but that's a hack.
                 if (k == N) {
                    stateSequence[k] = endState; // Force if last jump
                    System.err.println("Forcing state to endState as k=N.");
                 } else {
                    // A more robust solution might involve re-sampling N or throwing an error
                    System.err.println("Cannot determine next state. Aborting path generation.");
                    return;
                 }
            } else {
                stateSequence[k] = Randomizer.randomChoicePDF(transitionProbsToNextCandidates); 
            }
        }
        
        // Final check: S_N should be endState
        if (N > 0 && stateSequence[N] != endState) {
             // This can happen due to numerical precision or M_MAX_JUMPS issues
            System.err.println("Warning: Sampled stateSequence[N]=" + stateSequence[N] + " != endState=" + endState + ". Forcing S_N = endState.");
            stateSequence[N] = endState;
        } else if (N == 0 && startState != endState) {
             System.err.println("Error: N=0 sampled but startState != endState. Path impossible.");
             return;
        }


        // --- Step 4: Construct the Stochastic Map (and filter fictitious jumps) ---
        int currentActualState = startState;

        if (N == 0) {
            if (startState == endState) {
                 return;
            } else {
                // This case should have been caught by p_ij_t_denominator being zero
                System.err.println("Error: N=0 but startState != endState. Path construction failed.");
                return;
            }
        } else {
            for (int k = 0; k < N; k++) { // Iterate through N jumps, creating N+1 intervals
                double jumpOccursAt = jumpTimes[k];
                int stateBeforeThisJump = stateSequence[k]; // S_k
                int stateAfterThisJump = stateSequence[k+1]; // S_{k+1}

                if (stateBeforeThisJump != currentActualState) { // Should not happen if logic is correct
                     System.err.println("State tracking error before jump " + k);
                }

                if (stateBeforeThisJump != stateAfterThisJump) {
                	mutations.add(new MutationOnBranch(nodeNr, jumpOccursAt, stateBeforeThisJump, stateAfterThisJump, siteNr));
                }
                currentActualState = stateAfterThisJump;
            }
        }
	}

	public void setRatematrix(double[][] rateMatrixR) {
		int numStates = rateMatrixR.length;

		// --- Step 0: Precomputation & Initialization ---
		lambdaMax = 0.0;
		for (int i = 0; i < numStates; i++) {
			if (-rateMatrixR[i][i] > lambdaMax) {
				lambdaMax = -rateMatrixR[i][i];
			}
		}

		qUnif = getQUnif(rateMatrixR, lambdaMax);

		// Precompute powers of Q_unif to avoid re-computation
		qUnifPowers = new ArrayList<>();
		qUnifPowers.add(UniformisationStochasticMapping.identity(numStates)); // Q_unif^0

		for (int n = 0; n <= M_MAX_JUMPS; n++) {
			if (n > 0) {
				qUnifPowers.add(UniformisationStochasticMapping.multiply(qUnifPowers.get(n - 1), qUnif));
			}
		}
	}

	private double[] setUpWeights(double totalTime) {
		double[] weightsN = new double[M_MAX_JUMPS + 1];

		for (int n = 0; n <= M_MAX_JUMPS; n++) {

			double poissonTerm;
			if (lambdaMax * totalTime == 0 && n == 0) { // Handle 0^0 case for Poisson
				poissonTerm = Math.exp(-lambdaMax * totalTime);
			} else if (lambdaMax * totalTime == 0 && n > 0) {
				poissonTerm = 0;
			} else {
				poissonTerm = Math.exp(-lambdaMax * totalTime) * Math.pow(lambdaMax * totalTime, n)
						/ CombinatoricsUtils.factorial(n);
			}
			weightsN[n] = poissonTerm;
		}

		return weightsN;
	}

	
    /** return matrix I+R/lambdaMax **/
    double [][] getQUnif(double [][] rateMatrixR, double lambdaMax) {
//		return add(
//		        identity(numStates),
//		        multiplyByScalar(rateMatrixR, 1.0 / lambdaMax)
        int rows = rateMatrixR.length;
        int cols = rateMatrixR[0].length;
        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = rateMatrixR[i][j] /lambdaMax;
            }
            result[i][i] += 1.0;
        }
        return result;
    }
}
