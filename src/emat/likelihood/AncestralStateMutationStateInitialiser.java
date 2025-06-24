package emat.likelihood;



import java.util.List;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.inference.CalculationNode;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;
import emat.stochasticmapping.StochasticMapping;
import emat.stochasticmapping.TimeStateInterval;
import emat.stochasticmapping.UniformisationStochasticMapping;

@Description("Mutation State initialiser by treelikelihood based ancestral state reconstruction")
public class AncestralStateMutationStateInitialiser extends TreeLikelihood implements StateNodeInitialiser {
    public static final String STATES_KEY = "states";
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree to be initialised", Validate.REQUIRED);

    public Input<Boolean> useMAPInput = new Input<Boolean>("useMAP", "whether to use maximum aposteriori assignments or sample", false);
    public Input<Boolean> sampleTipsInput = new Input<Boolean>("sampleTips", "if tips have missing data/ambigous values sample them for logging (default true)", true);
    
	/** and node number associated with parameter **/
	int[] leafNr;

	int traitDimension;

    int patternCount;
    int stateCount;

    int[][] tipStates; // used to store tip states
    
    
    protected DataType dataType;
    private int[][] reconstructedStates;
    private boolean useMAP = false;
    
    @Override
    public void initAndValidate() {
    	if (dataInput.get().getSiteCount() == 0) {
    		return;
    	}
    	
    	
    	String sJavaOnly = null;
		sJavaOnly = System.getProperty("java.only");
		System.setProperty("java.only", "" + true);
    	super.initAndValidate();
    	if (sJavaOnly != null) {
    		System.setProperty("java.only", sJavaOnly);
    	} else {
    		System.clearProperty("java.only");
    	}
    	
        TreeInterface treeModel = treeInput.get();
        patternCount = dataInput.get().getPatternCount();
        dataType = dataInput.get().getDataType();
        stateCount = dataType.getStateCount();

        reconstructedStates = new int[treeModel.getNodeCount()][patternCount];

        this.useMAP = useMAPInput.get();

        int tipCount = treeModel.getLeafNodeCount();
        tipStates = new int[tipCount][];

        Alignment data = dataInput.get();
        for (Node node : treeInput.get().getExternalNodes()) {
            String taxon = node.getID();
            int taxonIndex = data.getTaxonIndex(taxon);
            if (taxonIndex == -1) {
            	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                    taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
                }
                if (taxonIndex == -1) {
                	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
                }
            }
            tipStates[node.getNr()] = new int[patternCount];
            if (!m_useAmbiguities.get()) {
            	likelihoodCore.getNodeStates(node.getNr(), tipStates[node.getNr()]);
            } else {
            	int [] states = tipStates[node.getNr()];
	            for (int i = 0; i < patternCount; i++) {
	                int code = data.getPattern(taxonIndex, i);
	                int[] statesForCode = data.getDataType().getStatesForCode(code);
	                if (statesForCode.length==1)
	                    states[i] = statesForCode[0];
	                else
	                    states[i] = code; // Causes ambiguous states to be ignored.
	            }

            }
    	}
        
        if (m_siteModel.getCategoryCount() > 1)
            throw new RuntimeException("Reconstruction not implemented for multiple categories yet.");
        
        
        // stuff for dealing with ambiguities in tips
		traitDimension = tipStates[0].length;

    }

    
	@Override
	public void initStateNodes() {
		
		optimiseClock();
		
		calculateLogP();
        // redraw states and return joint density of drawn states
        redrawAncestralStates();

		
		MutationState state = stateInput.get();
		Alignment data = dataInput.get();
		TreeInterface tree = treeInput.get();
		
		
		
		StochasticMapping mapper = new UniformisationStochasticMapping();
		((GeneralSubstitutionModel)substitutionModel).setupRateMatrix();
		double [][] rateMatrixR = ((GeneralSubstitutionModel)substitutionModel).getRateMatrix();
		
		// initialise internal node sequences 
		// + set mutations in the middle (if any)
		for (int nodeNr = 0; nodeNr < treeInput.get().getNodeCount() - 1; nodeNr++) {
			Node node = tree.getNode(nodeNr);
			double length = node.getLength();
			int [] nodePatternStates = reconstructedStates[nodeNr];
			int [] parentPatternStates = reconstructedStates[tree.getNode(nodeNr).getParent().getNr()];
			int [] nodeSequence = new int[data.getSiteCount()];
			int [] parentSequence = new int[data.getSiteCount()];
			for (int siteNr = 0; siteNr < data.getSiteCount(); siteNr++) {
				int k = data.getPatternIndex(siteNr);
				nodeSequence[siteNr] = nodePatternStates[k];
				parentSequence[siteNr] = parentPatternStates[k];
			}
			List<TimeStateInterval> path = mapper.generatePath(rateMatrixR, nodeSequence, parentSequence, length);
			for (int i = 0; i < path.size() - 1; i++) {
				if (path.get(i).endTime()<length) {
					state.addMutation(path.get(i).site(), nodeNr, path.get(i).endTime()/length, path.get(i).state(), path.get(i+1).state());
				}
			}
			state.setNodeSequence(nodeNr, nodeSequence);
		}
		
		// set root state
		int [] nodePatternStates = reconstructedStates[tree.getRoot().getNr()];
		int [] rootState = new int[data.getSiteCount()];
		for (int siteNr = 0; siteNr < data.getSiteCount(); siteNr++) {
			int k = data.getPatternIndex(siteNr);
			rootState[siteNr] = nodePatternStates[k];
		}
		state.setNodeSequence(tree.getNodeCount()-1, rootState);

		
		// self initialise MutationState 
		state.reinitialise();
		
		
		// save memory
		likelihoodCore = null;
		reconstructedStates = null;
		tipStates = null;
	}
    

    public void redrawAncestralStates() {
        logP = 0;
        TreeInterface tree = treeInput.get();
        traverseSample(tree, tree.getRoot(), null);
    }
    
    @Override
    public double calculateLogP() {
        logP = super.calculateLogP();
        return logP;
    }

    // optimise likelihood over mean clock rate for this tree
    private double optimiseClock() {
    	    	
		Function clockModel = branchRateModel.meanRateInput.get();
		if (!(clockModel instanceof RealParameter r)) {
			Log.warning("Warning: Cannot optimise clock model, since mean rate is not a RealParameter");
			Log.warning("Perhaps explicitly specify the StrictClockModel");
			return 0;
		}
		
		
		
        UnivariateFunction f = new UnivariateFunction() {
			@Override
			public double value(double x0) {
				double x = Math.exp(x0);
				r.setValue(x);
				((StrictClockModel)branchRateModel).requiresRecalculation();
				double v = calculateLogP();
				// System.err.println(x0 + " " + x + " " + v);
				return v;
			}
        	
        };
		
        double max = 2;
        double min = -15;
        BrentOptimizer optimizer = new BrentOptimizer(1e-4, 1e-5);
        UnivariatePointValuePair p =
		optimizer.optimize(new MaxEval(200),
                new UnivariateObjectiveFunction(f),
                GoalType.MAXIMIZE,
                new InitialGuess(new double [] {1.0}),
                new SearchInterval(min, max));
        double optimal = Math.exp(p.getPoint());
        
        Log.warning("Setting mean clock rate to " + optimal);
		r.setValue(optimal);
		((StrictClockModel)branchRateModel).requiresRecalculation();
		return p.getValue();
	}

	private int drawChoice(double[] measure) {
        if (useMAP) {
            double max = measure[0];
            int choice = 0;
            for (int i = 1; i < measure.length; i++) {
                if ((measure[i] - max)/(measure[i] + max) > 1e-10) {
                    max = measure[i];
                    choice = i;
                }
            }
            return choice;
        } else {
            return Randomizer.randomChoicePDF(measure);
        }
    }

    /**
     * Traverse (pre-order) the tree sampling the internal node states.
     *
     * @param tree        - TreeModel on which to perform sampling
     * @param node        - current node
     * @param parentState - character state of the parent node to 'node'
     */
    public void traverseSample(TreeInterface tree, Node node, int[] parentState) {

        int nodeNum = node.getNr();

        Node parent = node.getParent();

        // This function assumes that all partial likelihoods have already been calculated
        // If the node is internal, then sample its state given the state of its parent (pre-order traversal).

        double[] conditionalProbabilities = new double[stateCount];
        int[] state = new int[patternCount];

        if (!node.isLeaf()) {

            if (parent == null) {

                double[] rootPartials = m_fRootPartials;

                double[] rootFrequencies = substitutionModel.getFrequencies();
                if (rootFrequenciesInput.get() != null) {
                    rootFrequencies = rootFrequenciesInput.get().getFreqs();
                }

                // This is the root node
                for (int j = 0; j < patternCount; j++) {
                	System.arraycopy(rootPartials, j * stateCount, conditionalProbabilities, 0, stateCount);

                    for (int i = 0; i < stateCount; i++) {
                        conditionalProbabilities[i] *= rootFrequencies[i];
                    }
                    try {
                        state[j] = drawChoice(conditionalProbabilities);
                    } catch (Error e) {
                        System.err.println(e.toString());
                        System.err.println("Please report error to Marc");
                        state[j] = 0;
                    }
                    reconstructedStates[nodeNum][j] = state[j];

                    //System.out.println("Pr(j) = " + rootFrequencies[state[j]]);
                    logP += Math.log(rootFrequencies[state[j]]);
                }

            } else {

                // This is an internal node, but not the root
                double[] partialLikelihood = new double[stateCount * patternCount];
                likelihoodCore.getNodePartials(node.getNr(), partialLikelihood);


                for (int j = 0; j < patternCount; j++) {

                    int parentIndex = parentState[j] * stateCount;
                    int childIndex = j * stateCount;

                    for (int i = 0; i < stateCount; i++) {
                        conditionalProbabilities[i] = partialLikelihood[childIndex + i] * probabilities[parentIndex + i];
                    }

                    state[j] = drawChoice(conditionalProbabilities);
                    reconstructedStates[nodeNum][j] = state[j];
                    double contrib = probabilities[parentIndex + state[j]];
                    logP += Math.log(contrib);
                }
            }

            // Traverse down the two child nodes
            Node child1 = node.getChild(0);
            traverseSample(tree, child1, state);

            Node child2 = node.getChild(1);
            traverseSample(tree, child2, state);
        } else {

            // This is an external leaf
        	System.arraycopy(tipStates[nodeNum], 0, reconstructedStates[nodeNum], 0, reconstructedStates[nodeNum].length);
        	
        	if (sampleTipsInput.get()) {
	            // Check for ambiguity codes and sample them
	            for (int j = 0; j < patternCount; j++) {
	
	                final int thisState = reconstructedStates[nodeNum][j];
	                final int parentIndex = parentState[j] * stateCount;
	                likelihoodCore.getNodeMatrix(nodeNum, 0, probabilities);
	                if (dataType.isAmbiguousCode(thisState)) {
		                    
	                    boolean [] stateSet = dataType.getStateSet(thisState);
	                    for (int i = 0; i < stateCount; i++) {
	                        conditionalProbabilities[i] =  stateSet[i] ? probabilities[parentIndex + i] : 0;
	                    }
	                    
	                    reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);
	                }
	
	                double contrib = probabilities[parentIndex + reconstructedStates[nodeNum][j]];
	                logP += Math.log(contrib);
	            }
        	}
        	
        }
    }

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(stateInput.get());
	}
}
