package emat.likelihood;


import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;

@Description("Likelihood for a explicit mutation annotated tree")
public class MutationStateTreeLikelihood extends GenericTreeLikelihood {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree", Validate.REQUIRED);

	private TreeInterface tree;
	private MutationState state;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		tree = treeInput.get();
		state = stateInput.get();
	}
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}

}
