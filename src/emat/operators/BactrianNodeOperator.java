package emat.operators;

import java.text.DecimalFormat;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import emat.likelihood.EditableTree;
import emat.likelihood.MutationState;

@Description("Randomly selects true internal tree node (i.e. not the root) and move node height uniformly in interval " +
        "restricted by the nodes parent and children.")
public class BactrianNodeOperator extends EditableTreeOperator {
    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: larger means more bold proposals", 0.1);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());
    // root and internal nodes use different ways to sample
    // so auto-optimisation for the root differs from that of the internal nodes
    // therefore they should be separated out in different operators 
    final public Input<Boolean> rootOnlyInput = new Input<>("rootOnly", "flag to indicate that only the root should be operated on", false);

    final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree",
			Validate.REQUIRED);

    protected KernelDistribution kernelDistribution;


    protected double scaleFactor;
    private MutationState state;

	// empty constructor to facilitate construction by XML + initAndValidate
	public BactrianNodeOperator() {
	}
	
	public BactrianNodeOperator(Tree tree) {
	    try {
	        initByName(treeInput.getName(), tree);
	    } catch (Exception e) {
	        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
	        throw new RuntimeException("Failed to construct Uniform Tree Operator.");
	    }
	}
	
	@Override
	public void initAndValidate() {
    	kernelDistribution = kernelDistributionInput.get();
	    scaleFactor = scaleFactorInput.get();
	    state = stateInput.get();
	}

    /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	EditableTree tree = (EditableTree) InputUtil.get(treeInput, this);

        // randomly select internal node
        int nodeCount = tree.getNodeCount();
        
        // Abort if no non-root internal nodes
        if (tree.getInternalNodeCount()==1)
            return Double.NEGATIVE_INFINITY;
        
        if (rootOnlyInput.get()) {
        	Node node = tree.getRoot();
            double scale = kernelDistribution.getScaler(node.getNr(), node.getHeight(), getCoercableParameterValue());
            final double height = node.getHeight();
            final double newHeight = node.getHeight() * scale;

            if (newHeight < Math.max(node.getLeft().getHeight(), node.getRight().getHeight())) {
                return Double.NEGATIVE_INFINITY;
            }
            tree.setHeight(node.getNr(), newHeight);

            double logHR = Math.log(scale);
            // Since we also scale mutations in child branches, these need to be
            // taken in account for the HR as well.
            int mutationCountLeft = state.getMutationList(node.getLeft().getNr()).size();
            double leftHeight = node.getLeft().getHeight();
            logHR += mutationCountLeft * Math.log((newHeight - leftHeight)/(height - leftHeight));

            int mutationCountRight = state.getMutationList(node.getRight().getNr()).size();
            double rightHeight = node.getRight().getHeight();
            logHR += mutationCountRight * Math.log((newHeight - rightHeight)/(height - rightHeight));
            
            return logHR;
        }
        Node node;
        do {
            int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isLeaf() || node.isRoot());
        
        double upper = node.getParent().getHeight();
        double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        
        double scale = kernelDistribution.getScaler(0, Double.NaN, scaleFactor);

        // transform value
        double value = node.getHeight();
        double y = (upper - value) / (value - lower);
        y *= scale;
        double newValue = (upper + lower * y) / (y + 1.0);
        
        if (newValue < lower || newValue > upper) {
        	return Double.NEGATIVE_INFINITY;
        	//throw new RuntimeException("programmer error: new value proposed outside range");
        }
        
        tree.setHeight(node.getNr(), newValue);

        double logHR = Math.log(scale) + 2.0 * Math.log((newValue - lower)/(value - lower));

        // Since we also scale mutations in surrounding branches, these need to be
        // taken in account for the HR as well.
        int mutationCount = state.getMutationList(node.getNr()).size();
        double parentHeight = node.getParent().getHeight();
        logHR += mutationCount * Math.log((parentHeight - newValue)/(parentHeight - value));
        
        int mutationCountLeft = state.getMutationList(node.getLeft().getNr()).size();
        double leftHeight = node.getLeft().getHeight();
        logHR += mutationCountLeft * Math.log((newValue - leftHeight)/(value - leftHeight));

        int mutationCountRight = state.getMutationList(node.getRight().getNr()).size();
        double rightHeight = node.getRight().getHeight();
        logHR += mutationCountRight * Math.log((newValue - rightHeight)/(value - rightHeight));
         
        return logHR;
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(double value) {
    	scaleFactor = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
    	if (optimiseInput.get()) {
	        double delta = calcDelta(logAlpha);
	        double scaleFactor = getCoercableParameterValue();
	        delta += Math.log(scaleFactor);
	        scaleFactor = Math.exp(delta);
	        setCoercableParameterValue(scaleFactor);
    	}
    }
    
    @Override
    public double getTargetAcceptanceProbability() {
    	return 0.3;
    }
    
    @Override
    public String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = getCoercableParameterValue() * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10 || prob > 0.40) {
            return "Try setting scale factor to about " + formatter.format(newWindowSize);
        } else return "";
    }

}
