package emat.operators;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.util.CombinatoricsUtils;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import emat.likelihood.MutationOnBranch;
import emat.likelihood.MutationState;
import emat.likelihood.MutationStateTreeLikelihood;
import emat.substitutionmodel.EmatSubstitutionModel;

@Description("For testing purposes only. "
		+ "Operator that resamples sequence on a node using the substitution and clock model and "
		+ "resamples mutations on surrounding branches. "
		+ "Note this requires exponentiating the substitution matrix.")
public class SimpleMutationOnNodeResampler extends Operator {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree",
			Validate.REQUIRED);
	final public Input<MutationStateTreeLikelihood> likelihoodInput = new Input<>("likelihood",
			"tree likelihood from which substitution model and branch rate models are obtained for stochastic mapping",
			Validate.REQUIRED);

	protected MutationState state;
	protected GeneralSubstitutionModel substModel;
	protected BranchRateModel clockModel;
	protected int stateCount;

	protected int M_MAX_JUMPS = 20;
	
	protected double lambdaMax;
	protected double[][] qUnif;
	protected List<double[][]> qUnifPowers;

	private double[] weightsN, leftWeightsN, rightWeightsN;

    private double [][] pMatrix;
    private double [][] qMatrix;
    private double gamma;

	
	@Override
	public void initAndValidate() {
		state = stateInput.get();
		stateCount = state.getStateCount();

		substModel = (GeneralSubstitutionModel) ((SiteModel)likelihoodInput.get().siteModelInput.get()).substModelInput.get();
		clockModel = likelihoodInput.get().branchRateModelInput.get();

		qMatrix = substModel.getRateMatrix();
        this.gamma = calculateUniformisationRate(qMatrix);
        this.pMatrix = createDtmsTransitionMatrix(qMatrix, gamma);
	}

	public static double calculateUniformisationRate(double [][] qMatrix) {
        double maxRate = 0.0;
        for (int i = 0; i < qMatrix.length; i++) {
            maxRate = Math.max(maxRate, -qMatrix[i][i]);
        }
        return maxRate;
    }

    public static double[][] createDtmsTransitionMatrix(double [][] qMatrix, double gamma) {
    	int stateCount = qMatrix.length;
        double[][] p = new double[stateCount][stateCount];
        if (gamma == 0) {
            for(int i = 0; i < stateCount; i++) p[i][i] = 1.0;
            return p;
        }
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                if (i == j) {
                    p[i][j] = 1.0 + qMatrix[i][j] / gamma;
                } else {
                    p[i][j] = qMatrix[i][j] / gamma;
                }
            }
        }
        return p;
    }

	public static void generatePath(int nodeNr, int siteNr, int endState, int startState, double branchLength,
			List<MutationOnBranch> mutations, double [][]pMatrix, double gamma) { 

        // Rejection sampling loop
        while (true) {
            // 1. Simulate the number of potential events (N) from a Poisson distribution
            double lambda = gamma * branchLength;
            int numEvents = FastRandomiser.drawFromPoisson(lambda);
            
            if (numEvents == 0) {
            	if (startState == endState) {
                	return;
            	} else {
            		numEvents = 1;
            	}
            }
            
            // 2. Simulate the path of N steps using the P matrix
            List<Integer> path = new ArrayList<>();
            int currentState = startState;
            path.add(currentState);

            for (int i = 0; i < numEvents; i++) {
                // Get the probability distribution for the current state
                double[] probabilities = pMatrix[currentState];
                // Draw the next state from this distribution
                int prevState = currentState;
                currentState = FastRandomiser.drawFromCategorical(probabilities);
                if (currentState != prevState) {
                	path.add(currentState);
                }
            }

            // 3. Accept or reject the simulated path
            if (currentState == endState) {
                // Path is accepted, draw times, add to mutations and return
            	numEvents = path.size() - 1;
            	if (numEvents == 0) {
            		return;
            	}
            	double [] times = new double[numEvents];
            	for (int i = 0; i < numEvents; i++) {
            		times[i] = FastRandomiser.nextDouble();
            	}
            	Arrays.sort(times);
            	for (int i = 0; i < numEvents; i++) {
            		mutations.add(new MutationOnBranch(nodeNr, times[i], path.get(i), path.get(i+1), siteNr));
            	}            	
                return;
            }
            // If not, the loop continues and we try again.
        }
    }

	@Override
	public double proposal() {
		// randomly select internal node
		TreeInterface tree = state.treeInput.get();
		int nodeNr = tree.getLeafNodeCount() + FastRandomiser.nextInt(tree.getInternalNodeCount());
		Node node = tree.getNode(nodeNr);

		if (node.isRoot()) {
			resampleRoot(node);
		} else {
			resample(node);
		}
		

		return Double.POSITIVE_INFINITY;
	}


	private void resampleRoot(Node root) {
		if (true) {
//			return;
		}
		Node left = root.getLeft();
		Node right = root.getRight();
		int[] leftStates = state.getNodeSequence(left.getNr());
		int[] rightStates = state.getNodeSequence(right.getNr());
		
		
		double [] Pleft, Pright;
		Pleft = new double[stateCount * stateCount];
		Pright = new double[stateCount * stateCount];

		double [] freqs = substModel.getFrequencies();
		
		substModel.getTransitionProbabilities(root.getLeft(), root.getLeft().getHeight(), root.getHeight(), clockModel.getRateForBranch(root.getLeft()), Pleft);
		substModel.getTransitionProbabilities(root.getRight(), root.getRight().getHeight(), root.getHeight(), clockModel.getRateForBranch(root.getRight()), Pright);		
		
		List<MutationOnBranch> branchMutationsLeft = new ArrayList<>();
		List<MutationOnBranch> branchMutationsRight = new ArrayList<>();

		double totalTimeLeft = root.getLeft().getLength() * clockModel.getRateForBranch(root.getLeft());
		double totalTimeRight = root.getRight().getLength() * clockModel.getRateForBranch(root.getRight());
		int [] nodeSequence = state.getNodeSequenceForUpdate(root.getNr());
		double [] pNodeState = new double[stateCount];

		for (int i = 0; i < nodeSequence.length; i++) {
			for (int nodeState = 0; nodeState < stateCount; nodeState++) {
				pNodeState[nodeState] = freqs[nodeState] * 
						Pleft[nodeState * stateCount + leftStates[i]] * 
						Pright[nodeState * stateCount + rightStates[i]];
			}
			int nodeState = FastRandomiser.randomChoicePDF(pNodeState);
			nodeSequence[i] = nodeState;

			generatePath(root.getLeft().getNr(), i, nodeState, leftStates[i], totalTimeLeft, branchMutationsLeft, pMatrix, gamma);
			generatePath(root.getRight().getNr(), i, nodeState, rightStates[i], totalTimeRight, branchMutationsRight, pMatrix, gamma);
	}

		Collections.sort(branchMutationsLeft);
		Collections.sort(branchMutationsRight);
		
		
		state.setBranchMutations(left.getNr(), branchMutationsLeft);
		state.setBranchMutations(right.getNr(), branchMutationsRight);
		
		
	}

	protected void resample(Node node) {
		int nodeNr = node.getNr();
		substModel.setupRateMatrix();
		setRatematrix(substModel.getRateMatrix());

		int[] states = state.getNodeSequence(node.getParent().getNr());
		int[] leftStates = state.getNodeSequence(node.getLeft().getNr());
		int[] rightStates = state.getNodeSequence(node.getRight().getNr());
		int [] nodeSequence = state.getNodeSequenceForUpdate(nodeNr);

		double totalTime = node.getLength() * clockModel.getRateForBranch(node);
		double totalTimeLeft = node.getLeft().getLength() * clockModel.getRateForBranch(node.getLeft());
		double totalTimeRight = node.getRight().getLength() * clockModel.getRateForBranch(node.getRight());

		// TODO: cach weights per branch & update only when evolutionary distance = (branch lengths * clock rate) changes
		weightsN = setUpWeights(totalTime);
		leftWeightsN = setUpWeights(totalTimeLeft);
		rightWeightsN = setUpWeights(totalTimeRight);
		
		
		List<MutationOnBranch> branchMutations = new ArrayList<>();
		List<MutationOnBranch> branchMutationsLeft = new ArrayList<>();
		List<MutationOnBranch> branchMutationsRight = new ArrayList<>();
		
		double [] P, Pleft, Pright;
		P = new double[stateCount * stateCount];
		Pleft = new double[stateCount * stateCount];
		Pright = new double[stateCount * stateCount];
		substModel.getTransitionProbabilities(node, node.getHeight(), node.getParent().getHeight(), clockModel.getRateForBranch(node), P);
		substModel.getTransitionProbabilities(node.getLeft(), node.getLeft().getHeight(), node.getHeight(), clockModel.getRateForBranch(node.getLeft()), Pleft);
		substModel.getTransitionProbabilities(node.getRight(), node.getRight().getHeight(), node.getHeight(), clockModel.getRateForBranch(node.getRight()), Pright);
				
		double [] pNodeState = new double[stateCount];

		for (int i = 0; i < states.length; i++) {
			for (int nodeState = 0; nodeState < stateCount; nodeState++) {
				pNodeState[nodeState] = P[states[i] * stateCount + nodeState] * 
						Pleft[nodeState * stateCount + leftStates[i]] * 
						Pright[nodeState * stateCount + rightStates[i]];
			}
			int nodeState = FastRandomiser.randomChoicePDF(pNodeState);
			nodeSequence[i] = nodeState;
			
			generatePath(nodeNr, i, states[i], nodeState, totalTime, branchMutations, pMatrix, gamma);
			generatePath(node.getLeft().getNr(), i, nodeState, leftStates[i], totalTimeLeft, branchMutationsLeft, pMatrix, gamma);
			generatePath(node.getRight().getNr(), i, nodeState, rightStates[i], totalTimeRight, branchMutationsRight, pMatrix, gamma);
		}

		state.setBranchMutations(nodeNr, branchMutations);
		state.setBranchMutations(node.getLeft().getNr(), branchMutationsLeft);
		state.setBranchMutations(node.getRight().getNr(), branchMutationsRight);
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
		qUnifPowers.add(EmatSubstitutionModel.identity(numStates)); // Q_unif^0

		for (int n = 0; n <= M_MAX_JUMPS; n++) {
			if (n > 0) {
				qUnifPowers.add(EmatSubstitutionModel.multiply(qUnifPowers.get(n - 1), qUnif));
			}
		}
	}

	protected double[] setUpWeights(double totalTime) {
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
