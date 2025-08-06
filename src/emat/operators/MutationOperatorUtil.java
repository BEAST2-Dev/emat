package emat.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.CombinatoricsUtils;

import beast.base.evolution.branchratemodel.BranchRateModel;
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
	
	
	public static List<MutationOnBranch> resample(
			Node node,
			MutationState state,
			BranchRateModel clockModel,
			List<double[][]> qUnifPowers,
			double lambdaMax, int M_MAX_JUMPS) {
		int[] states = state.getNodeSequence(node.getNr());
		int[] parentStates = state.getNodeSequence(node.getParent().getNr());
		
		double totalTime = node.getLength() * clockModel.getRateForBranch(node);

		double [] weightsN = MutationOperatorUtil.setUpWeights(totalTime, lambdaMax, M_MAX_JUMPS);
		
		List<MutationOnBranch> branchMutations = new ArrayList<>();

		for (int i = 0; i < states.length; i++) {
			double [] p = new double[M_MAX_JUMPS];
			for (int r = 0; r < M_MAX_JUMPS; r++) {
				p[r] = weightsN[r] * qUnifPowers.get(r)[states[i]][parentStates[i]];
			}			
			int N = FastRandomiser.randomChoicePDF(p);
			MutationOperatorUtil.generatePath(node.getNr(), i, parentStates[i], states[i], N, qUnifPowers ,branchMutations);
		}
		return branchMutations;
	}

	
	
	static public void generatePath(int nodeNr, int siteNr, int endState, int startState, int N,
			List<double[][]> qUnifPowers,
			List<MutationOnBranch> mutations) {
        // --- Step 2: Sample the Jump Times (tau_1, ..., tau_N) ---
        double[] jumpTimes = new double[N];
        for (int i = 0; i < N; i++) {
            jumpTimes[i] = FastRandomiser.nextDouble();
        }
        Arrays.sort(jumpTimes);
        
        // --- Step 3: Sample the Sequence of States (S_0, ..., S_N) ---
        int[] stateSequence = new int[N + 1];
        stateSequence[0] = startState;
        stateSequence[N] = endState;
        
        double [][] qUnif = qUnifPowers.get(1);
        int stateCount = qUnif.length;

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
                stateSequence[k] = FastRandomiser.randomChoicePDF(transitionProbsToNextCandidates); 
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
	
	public static double[] setUpWeights(double totalTime, double lambdaMax, int M_MAX_JUMPS) {
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


	public static int [] splitBranchMutations(MutationState state, Node targetBranch,
			List<MutationOnBranch> newTargetMutations, List<MutationOnBranch> newParentMutations, double newHeight) {
		// split target branch mutations
		int targetNr = targetBranch.getNr();
		int parentNr = targetBranch.getParent().getNr();
		List<MutationOnBranch> targetMutations = state.getMutationList(targetNr);
		double targetLength = targetBranch.getLength();
		double f = targetLength / (newHeight - targetBranch.getHeight());
		
		double f2 = (targetBranch.getParent().getHeight() - targetBranch.getHeight()) / (targetBranch.getParent().getHeight() - newHeight); 
		double threshold = 1 / f;
		int [] states0 = state.getNodeSequence(targetNr);
		int [] states = state.getNodeSequenceForUpdate(parentNr);
		int [] parentNodeStates = states.clone();
		System.arraycopy(states0, 0, states, 0, states.length);
		int [] statesOrig = null;

		for (MutationOnBranch m : targetMutations) {
			if (m.getBrancheFraction() < threshold) {
				newTargetMutations.add(new MutationOnBranch(targetNr, f * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
				states[m.siteNr()] = m.getToState();
			} else {
				newParentMutations.add(new MutationOnBranch(parentNr, f2 * (m.brancheFraction() - threshold), m.getFromState(), m.getToState(), m.siteNr()));
			}
		}
		return states;
	}

}
