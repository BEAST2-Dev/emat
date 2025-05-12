package emat.likelihood;


import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;

@Description("Likelihood for a explicit mutation annotated tree")
public class EMATTreeLikelihood extends GenericTreeLikelihood {
	

	private TreeInterface tree;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		tree = treeInput.get();
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
