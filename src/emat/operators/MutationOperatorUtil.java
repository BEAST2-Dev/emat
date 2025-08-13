package emat.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.util.CombinatoricsUtils;

import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import emat.likelihood.MutationOnBranch;
import emat.likelihood.MutationState;
import emat.substitutionmodel.EmatSubstitutionModel;

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

	private static void resample(Node node, Node virtualparent, int M_MAX_JUMPS,
			List<MutationOnBranch> branchMutations,
			List<MutationOnBranch> branchMutationsLeft,
			List<MutationOnBranch> branchMutationsRight,
			int [] nodeSequence,
			MutationState state,
			EmatSubstitutionModel substModel,
			double totalTime, double totalTimeLeft, double totalTimeRight
			) {
		
		int nodeNr = node.getNr();
		int stateCount = substModel.getStateCount();
		
		int[] states = state.getNodeSequence(virtualparent.getNr());
		int[] leftStates = state.getNodeSequence(node.getLeft().getNr());
		int[] rightStates = state.getNodeSequence(node.getRight().getNr());
		//int [] nodeSequence = state.getNodeSequenceForUpdate(nodeNr);


		int Relevant_M_MAX_JUMPS = 3;
		M_MAX_JUMPS = Relevant_M_MAX_JUMPS;
		
		
		// TODO: cache weights per branch & update only when evolutionary distance = (branch lengths * clock rate) changes
		double lambdaMax = substModel.getLambdaMax();
		List<double[][]> qUnifPowers = substModel.getQUnifPowers();
		double [] weightsN = MutationOperatorUtil.setUpWeights(totalTime, lambdaMax, M_MAX_JUMPS);
		double [] leftWeightsN = MutationOperatorUtil.setUpWeights(totalTimeLeft, lambdaMax, M_MAX_JUMPS);
		double [] rightWeightsN = MutationOperatorUtil.setUpWeights(totalTimeRight, lambdaMax, M_MAX_JUMPS);
		
		
		
		double [] pNodeState = new double[stateCount];

		double [][][] pNode = new double[stateCount][stateCount][M_MAX_JUMPS];
		double [][][] pLeft = new double[stateCount][stateCount][M_MAX_JUMPS];
		double [][][] pRight = new double[stateCount][stateCount][M_MAX_JUMPS];
		for (int src = 0; src < stateCount; src++) {
			for (int trgt = 0; trgt < stateCount; trgt++) {
				for (int r = 0; r < M_MAX_JUMPS; r++) {
					pNode[src][trgt][r] = weightsN[r] * qUnifPowers.get(r)[src][trgt];
					pLeft[src][trgt][r] = leftWeightsN[r] * qUnifPowers.get(r)[src][trgt];
					pRight[src][trgt][r] = rightWeightsN[r] * qUnifPowers.get(r)[src][trgt];
				}			
			}
		}

		int [] skipCount = new int[stateCount];
		int [] skipCountLeft = new int[stateCount];
		int [] skipCountRight = new int[stateCount];
		for (int i = 0; i < stateCount; i++) {
			skipCount[i] = drawSkipCount(pNode[i][i]);
			skipCountLeft[i] = drawSkipCount(pLeft[i][i]);
			skipCountRight[i] = drawSkipCount(pRight[i][i]);
		}
		
		
		for (int i = 0; i < states.length; i++) {
			for (int nodeState = 0; nodeState < stateCount; nodeState++) {
				double p = 0;
				for (int r = 0; r < Relevant_M_MAX_JUMPS; r++) {
					double p0 = weightsN[r] * qUnifPowers.get(r)[states[i]][nodeState];
					for (int s = 0; s < Relevant_M_MAX_JUMPS; s++) {
						double p1 = leftWeightsN[s] * qUnifPowers.get(s)[nodeState][leftStates[i]] * p0;
						for (int t = 0; t < Relevant_M_MAX_JUMPS; t++) {
							p += p1 * rightWeightsN[t] * qUnifPowers.get(t)[nodeState][rightStates[i]];
						}
					}
				}
				pNodeState[nodeState] = p;
			}
			int nodeState = FastRandomiser.randomChoicePDF(pNodeState);
			nodeSequence[i] = nodeState;
			
		}
		
		generatePaths(nodeSequence, states, pNode, branchMutations, nodeNr, skipCount, qUnifPowers);
		generatePaths(leftStates, nodeSequence, pLeft, branchMutationsLeft, node.getLeft().getNr(), skipCountLeft, qUnifPowers);
		generatePaths(rightStates, nodeSequence, pRight, branchMutationsRight, node.getRight().getNr(), skipCountRight, qUnifPowers);
	}
	
	
	private static void generatePaths(int[] nodeSequence, int[] states, double[][][] pNode,
			List<MutationOnBranch> branchMutations,
			int nodeNr, int [] skipCount, List<double[][]> qUnifPowers) {

//		for (int i = 0; i < nodeSequence.length; i++) {
//			int nodeState = nodeSequence[i];
//			int N = FastRandomiser.randomChoicePDF(pNode[states[i]][nodeState]);
//			generatePath(nodeNr, i, states[i], nodeState, N, qUnifPowers, branchMutations);
//		}
		
		for (int i = 0; i < nodeSequence.length; i++) {
			int nodeState = nodeSequence[i];
			if (states[i] == nodeState) {
				if (skipCount[nodeState] > 0) {
					skipCount[nodeState]--;
				} else {
					int N = FastRandomiser.randomChoicePDF(pNode[states[i]][nodeState]);
					generatePath(nodeNr, i, states[i], nodeState, N, qUnifPowers, branchMutations);
					skipCount[nodeState] = drawSkipCount(pNode[nodeState][nodeState]);
				}
			} else {
				int N = FastRandomiser.randomChoicePDF(pNode[states[i]][nodeState]);
				generatePath(nodeNr, i, states[i], nodeState, N, qUnifPowers, branchMutations);
			}
		}
		
	}


	// draw number of sites to skip when parent and node sites are the same
	private static int drawSkipCount(double[] pNode) {
		double p = pNode[0];
		double sum = 0;
		for (double d : pNode) {
			sum += d;
		}
		p /= sum;
		double u = FastRandomiser.nextDouble();
		int skip = (int) Math.floor(Math.log(u) / Math.log(p));
		if (skip < 0) {
			// numerical instability
			skip = Integer.MAX_VALUE;
		}
		return skip;
	}


	protected static void resample(Node node, final int M_MAX_JUMPS,
		List<MutationOnBranch> branchMutations,
		List<MutationOnBranch> branchMutationsLeft,
		List<MutationOnBranch> branchMutationsRight,
		int [] nodeSequence,
		MutationState state,
		EmatSubstitutionModel substModel,
		BranchRateModel clockModel
		) {

		double totalTime = node.getLength() * clockModel.getRateForBranch(node);
		double totalTimeLeft = node.getLeft().getLength() * clockModel.getRateForBranch(node.getLeft());
		double totalTimeRight = node.getRight().getLength() * clockModel.getRateForBranch(node.getRight());
		
		resample(node, node.getParent(), M_MAX_JUMPS, branchMutations, branchMutationsLeft, branchMutationsRight, nodeSequence, state, substModel, totalTime, totalTimeLeft, totalTimeRight);
	}

	protected static void resampleRoot(Node root, final int M_MAX_JUMPS,
			List<MutationOnBranch> branchMutationsLeft,
			List<MutationOnBranch> branchMutationsRight,
			int [] states,
			MutationState state,
			EmatSubstitutionModel substModel,
			BranchRateModel clockModel) {

		Node left = root.getLeft();
		Node right = root.getRight();
		int[] leftStates = state.getNodeSequence(left.getNr());
		int[] rightStates = state.getNodeSequence(right.getNr());
		
		double totalTimeLeft = left.getLength() * clockModel.getRateForBranch(left);
		double totalTimeRight = right.getLength() * clockModel.getRateForBranch(right);

		double totalTime = totalTimeLeft + totalTimeRight;
		double [] weightsN = MutationOperatorUtil.setUpWeights(totalTime, substModel.getLambdaMax(), M_MAX_JUMPS);
		
		int stateCount = substModel.getStateCount();
		List<double[][]> qUnifPowers = substModel.getQUnifPowers();

		double [][][] pNode = new double[stateCount][stateCount][M_MAX_JUMPS];
		for (int src = 0; src < stateCount; src++) {
			for (int trgt = 0; trgt < stateCount; trgt++) {
				for (int r = 0; r < M_MAX_JUMPS; r++) {
					pNode[src][trgt][r] = weightsN[r] * qUnifPowers.get(r)[src][trgt];
				}			
			}
		}

		int [] skipCount = new int[stateCount];
		for (int i = 0; i < stateCount; i++) {
			skipCount[i] = drawSkipCount(pNode[i][i]);
		}
		
		List<MutationOnBranch> branchMutations = new ArrayList<>();
		generatePaths(leftStates, rightStates, pNode, branchMutations, left.getNr(), skipCount, qUnifPowers);
		
//		for (int i = 0; i < leftStates.length; i++) {
//			double [] p = pNode[leftStates[i]][rightStates[i]];
//			int N = FastRandomiser.randomChoicePDF(p);
//			MutationOperatorUtil.generatePath(left.getNr(), i, rightStates[i], leftStates[i], N, qUnifPowers ,branchMutations);
//		}

		Collections.sort(branchMutations);
		
		System.arraycopy(leftStates, 0, states, 0, states.length);
		
		double leftFraction = totalTimeLeft / totalTime;
		double rightFraction = 1 - leftFraction;
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
	}

	protected static void resampleBelowRoot(Node node, final int M_MAX_JUMPS,
			List<MutationOnBranch> branchMutationsLeftOfRoot,
			List<MutationOnBranch> branchMutationsRightOfRoot,
			List<MutationOnBranch> branchMutationsLeft,
			List<MutationOnBranch> branchMutationsRight,
			int [] rootSequence,
			int [] nodeSequence,
			MutationState state,
			EmatSubstitutionModel substModel,
			BranchRateModel clockModel
			) {

		    // first, resample node as if the branch above goes all the way around the root to the sibling
			Node root = node.getParent();
			if (!root.isRoot()) {
				throw new IllegalArgumentException("Node must have root as parent");
			}
			Node leftOfRoot = root.getLeft();
			Node rightOfRoot = root.getRight();
			double timeLeftOfRoot = leftOfRoot.getLength() * clockModel.getRateForBranch(leftOfRoot);
			double timeRightOfRoot = rightOfRoot.getLength() * clockModel.getRateForBranch(rightOfRoot); 
			double totalTime = timeLeftOfRoot + timeRightOfRoot;
			double totalTimeLeft = node.getLeft().getLength() * clockModel.getRateForBranch(node.getLeft());
			double totalTimeRight = node.getRight().getLength() * clockModel.getRateForBranch(node.getRight());

			List<MutationOnBranch> branchMutations = new ArrayList<>();
			Node sibling = leftOfRoot == node ? rightOfRoot : leftOfRoot;
			resample(node, sibling, M_MAX_JUMPS, branchMutations, branchMutationsLeft, branchMutationsRight, nodeSequence, state, substModel, totalTime, totalTimeLeft, totalTimeRight);

			// next, redistribute the mutations, starting at nodeSequence going to the sibling
			double nodeFraction = (leftOfRoot == node ? timeLeftOfRoot : timeRightOfRoot) / totalTime;
			double siblingFraction = 1 - nodeFraction;
			List<MutationOnBranch> nodeMutations = null;
			List<MutationOnBranch> siblingMutations = null;
			if (leftOfRoot == node) {
				nodeMutations = branchMutationsLeftOfRoot;
				siblingMutations = branchMutationsRightOfRoot;
			} else {
				nodeMutations = branchMutationsRightOfRoot;
				siblingMutations = branchMutationsLeftOfRoot;
			}

			System.arraycopy(nodeSequence, 0, rootSequence, 0, nodeSequence.length);
			for (MutationOnBranch m : branchMutations) {
				if (m.brancheFraction() < nodeFraction) {
					m.setBrancheFraction(m.getBrancheFraction() / nodeFraction);
					nodeMutations.add(m);
					rootSequence[m.siteNr()] = m.getToState();
				} else {
					siblingMutations.add(
							new MutationOnBranch(sibling.getNr(), (1-m.getBrancheFraction()) / siblingFraction, m.getToState(), m.getFromState(), m.siteNr()));
				}
			}
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

		if (lambdaMax * totalTime == 0) {
			weightsN[0] = 1;
			return weightsN;
		}
		
		for (int n = 0; n <= M_MAX_JUMPS; n++) {

			double poissonTerm;
//			if (lambdaMax * totalTime == 0 && n == 0) { // Handle 0^0 case for Poisson
//				poissonTerm = Math.exp(-lambdaMax * totalTime);
//			} else if (lambdaMax * totalTime == 0 && n > 0) {
//				poissonTerm = 0;
//			} else {
				poissonTerm = Math.exp(-lambdaMax * totalTime) * Math.pow(lambdaMax * totalTime, n)
						/ CombinatoricsUtils.factorial(n);
//			}
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
