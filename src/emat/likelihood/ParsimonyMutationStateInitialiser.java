package emat.likelihood;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beastlabs.parsimony.ParsimonyCriterion;
import beastlabs.parsimony.FitchParsimonyFactory;

@Description("Initialises a mutation state based on parsimony")
public class ParsimonyMutationStateInitialiser extends CalculationNode implements StateNodeInitialiser {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree to be initialised", Validate.REQUIRED);
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void initStateNodes() {
		MutationState state = stateInput.get();
		Alignment data = state.dataInput.get();
		int stateCount = data.getMaxStateCount();
		
		ParsimonyCriterion parsimony = FitchParsimonyFactory.newFitchParsimony(data, false);
		TreeInterface tree = treeInput.get();
		parsimony.getScore((Tree) tree);
		
		// add mutations
		for (int nodeNr = 0; nodeNr < tree.getNodeCount()-1; nodeNr++) {
			int [] nodePatternStates = parsimony.getStates((Tree)tree, tree.getNode(nodeNr));
			int [] parentPatternStates = parsimony.getStates((Tree) tree, (tree.getNode(nodeNr).getParent()));
			int [] nodeSequence = new int[data.getSiteCount()];
			for (int siteNr = 0; siteNr < data.getSiteCount(); siteNr++) {
				int k = data.getPatternIndex(siteNr);
				nodeSequence[siteNr] = nodePatternStates[k];
				if (nodePatternStates[k] != parentPatternStates[k]) {
					state.addMutation(siteNr, nodeNr, 0.5f, parentPatternStates[k], nodePatternStates[k]);
				}
			}
			state.setNodeSequence(nodeNr, nodeSequence);
		}
		
		// set root state
		int [] nodePatternStates = parsimony.getStates((Tree)tree, tree.getRoot());
		int [] rootState = new int[data.getSiteCount()];
		for (int siteNr = 0; siteNr < data.getSiteCount(); siteNr++) {
			int k = data.getPatternIndex(siteNr);
			rootState[siteNr] = nodePatternStates[k];
		}
		state.setNodeSequence(tree.getNodeCount()-1, rootState);
		
		state.reinitialise();
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(stateInput.get());
	}

}
