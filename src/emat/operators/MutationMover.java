package emat.operators;

import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.operator.kernel.KernelDistribution;
import emat.likelihood.*;

@Description("Randomly picks a mutation and moves it up or down on the same branch")
public class MutationMover extends Operator {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree", Validate.REQUIRED);
	final public Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());

    final public Input<Double> deltaInput = new Input<>("delta", "Magnitude of change for two randomly picked values.", 0.25);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

	protected MutationState mutationState;
	protected KernelDistribution kernelDistribution;
	protected double delta;
	protected boolean optimise;
	
	@Override
	public void initAndValidate() {
		mutationState = stateInput.get();
		kernelDistribution = kernelDistributionInput.get();
		delta = deltaInput.get();
		optimise = optimiseInput.get();
	}

	@Override
	public double proposal() {
		int [] ids = new int[2];
		MutationOnBranch mutation = mutationState.getRandomMutation(ids);
		int siteNr = ids[0]; 
		int nodeNr = ids[1]; 
		double oldBranchFraction = mutation.getBrancheFraction();
		double newBranchFraction = oldBranchFraction + (float) kernelDistribution.getRandomDelta(0, 0, delta);

		if (newBranchFraction < 0 || newBranchFraction > 1) {
			// invalid branch fraction
			return Double.NEGATIVE_INFINITY;
		}
		
		mutationState.moveBranchFraction(siteNr, nodeNr, mutation, oldBranchFraction, newBranchFraction);
		
		return 0;
	}


	@Override
    public double getCoercableParameterValue() {
        return delta;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        delta = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    @Override
    public void optimize(final double logAlpha) {
        // must be overridden by operator implementation to have an effect
        if (optimise) {
            double _delta = calcDelta(logAlpha);
            _delta += Math.log(delta);
            delta = Math.exp(_delta);
        }

    }

    @Override
    public final String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double newDelta = delta * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else return "";
    }
    
	
    @Override
    public double getTargetAcceptanceProbability() {
    	return 0.3;
    }

}
