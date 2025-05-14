package emat.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;
import emat.likelihood.MutationOnBranch;

@Description("Slide mutations between branches")
public class MutationSlider extends MutationMover {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    private TreeInterface tree;
    
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		tree = treeInput.get();
	}

	@Override
	public double proposal() {
		int attempts = 0;
		double logHR = Double.NEGATIVE_INFINITY;

		while (Double.isInfinite(logHR) && attempts < 1000) {
			if (Randomizer.nextInt(3) == 0) {
				logHR = moveUp();
			} else {
				logHR = moveDown();
			}
			attempts++;
		}
		return logHR;
	}

	private double moveUp() {
		int [] ids = new int[2];
		MutationOnBranch mutation = mutationState.getRandomMutation(ids);
		int siteNr = ids[0]; 
		int nodeNr = ids[1]; 

		Node node = tree.getNode(nodeNr);
		Node parent = node.getParent();
		if (parent.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}
		Node other = parent.getLeft() == node ? parent.getRight() : parent.getLeft();
		Node grandParent = parent.getParent();
		int nodeChar = mutationState.getCharAt(nodeNr, siteNr);
		int parentChar = mutationState.getCharAt(parent.getNr(), siteNr);
		int grandParentChar = mutationState.getCharAt(grandParent.getNr(), siteNr);
		int otherChar = mutationState.getCharAt(other.getNr(), siteNr);
		
		
		return 0;
	}

	private double moveDown() {
		int [] ids = new int[2];
		MutationOnBranch mutation = mutationState.getRandomMutation(ids);
		int siteNr = ids[0]; 
		int nodeNr = ids[1]; 

		Node node = tree.getNode(nodeNr);
		if (node.isLeaf()) {
			return Double.NEGATIVE_INFINITY;
		}

		
		return 0;
	}
}
