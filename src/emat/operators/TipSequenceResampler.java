package emat.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import emat.likelihood.MutationOnBranch;
import emat.likelihood.MutationState;
import emat.likelihood.MutationStateTreeLikelihood;
import emat.substitutionmodel.EmatSubstitutionModel;

@Description("Resample any missing data on tips")
public class TipSequenceResampler extends Operator {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree",
			Validate.REQUIRED);
	final public Input<MutationStateTreeLikelihood> likelihoodInput = new Input<>("likelihood",
			"tree likelihood from which substitution model and branch rate models are obtained for stochastic mapping",
			Validate.REQUIRED);
	
	protected MutationState state;
	protected EmatSubstitutionModel substModel;
	protected BranchRateModel clockModel;
	protected int stateCount, taxonCount, siteCount;
	private TreeInterface tree;

    private double [][] pMatrix;
    private double [][] qMatrix;
    private double gamma;

    private Set<Integer> [] missing;

	@Override
	public void initAndValidate() {
		state = stateInput.get();
		tree = state.treeInput.get();
		stateCount = state.getStateCount();

		substModel = likelihoodInput.get().substModelInput.get();
		clockModel = likelihoodInput.get().branchRateModelInput.get();

		qMatrix = substModel.getRateMatrix();
        this.gamma = SimpleMutationOnNodeResampler.calculateUniformisationRate(qMatrix);
        this.pMatrix = SimpleMutationOnNodeResampler.createDtmsTransitionMatrix(qMatrix, gamma);

        // collect missing data information
        Alignment data = state.dataInput.get();
        DataType datatype = data.getDataType();
        taxonCount = data.getTaxonCount();
        siteCount = data.getSiteCount();
        missing = new Set[taxonCount];
        for (int i = 0; i < taxonCount; i++) {
        	missing[i] = new HashSet<>();
        }
    	for (int j = 0; j < siteCount; j++) {
    		int [] pattern = data.getPattern(data.getPatternIndex(j));
    		for (int k = 0; k < taxonCount; k++) {
    			if (datatype.isAmbiguousCode(pattern[k])) {
    				missing[k].add(j);
    			}
    		}
    	}
	}

	@Override
	public double proposal() {
		
		//for (int k = 0; k < taxonCount; k++) {
		int k = FastRandomiser.nextInt(taxonCount); {
			if (missing[k].size() > 0) {
				Set<Integer> sites = missing[k];
				
				// clear mutations ending in a sampled site for this taxon
				List<MutationOnBranch> mutations = new ArrayList<>();
				mutations.addAll(state.getMutationList(k));
				for (int i = mutations.size()-1; i >= 0; i--) {
					MutationOnBranch m = mutations.get(i);
					if (sites.contains(m.siteNr())) {
						mutations.remove(i);
					}
				}
			
				// sample new trajectory for the missing sites starting at the (known) parent site
				Node node = tree.getNode(k);
//System.err.println("tip sampling taxon " + k + " " + node.getID());
				int [] parentSequence = state.getNodeSequence(node.getParent().getNr());
				int [] sequence = state.getNodeSequenceForUpdate(k);
				double distance = node.getLength() * clockModel.getRateForBranch(node); 
				for (int site : missing[k]) {
					sequence[site] = generatePath(k, site, parentSequence[site], distance, mutations, pMatrix, gamma);
				}
				
				// set branch mutations and update statistics in state
				state.setBranchMutations(k, mutations);
			}
		}
		return Double.POSITIVE_INFINITY;
	}
	
	/** generate path for a given end state (at top of branch)
	 * @param nodeNr
	 * @param siteNr
	 * @param endState
	 * @param distance
	 * @param mutations
	 * @param pMatrix
	 * @param gamma
	 * @return start state (at bottom of branch)
	 */
	public static int generatePath(int nodeNr, int siteNr, int endState, double distance,
		List<MutationOnBranch> mutations, double [][]pMatrix, double gamma) { 

        // 1. Simulate the number of potential events (N) from a Poisson distribution
        double lambda = gamma * distance;
        int numEvents = FastRandomiser.drawFromPoisson(lambda);
        
        if (numEvents == 0) {
        	return endState;
        }

        // 2. Simulate the path of N steps using the P matrix
        List<Integer> path = new ArrayList<>();
        int currentState = endState;
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

        // Draw times, add to mutations and return
    	numEvents = path.size() - 1;
    	if (numEvents == 0) {
    		return endState;
    	}
    	double [] times = new double[numEvents];
    	for (int i = 0; i < numEvents; i++) {
    		times[i] = FastRandomiser.nextDouble();
    	}
    	Arrays.sort(times);
    	for (int i = 0; i < numEvents; i++) {
    		mutations.add(new MutationOnBranch(nodeNr, times[numEvents-1-i], path.get(i+1), path.get(i), siteNr));
    	}            	
        return currentState;
    }

}