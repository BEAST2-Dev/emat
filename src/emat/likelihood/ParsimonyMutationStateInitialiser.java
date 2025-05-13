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
import beastlabs.parsimony.FitchParsimony;

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
		
		FitchParsimony parsimony = new FitchParsimony(data, false);
		TreeInterface tree = treeInput.get();
		parsimony.getScore((Tree) tree);
		
		// add mutations
		for (int nodeNr = 0; nodeNr < tree.getNodeCount()-1; nodeNr++) {
			int [] nodePatternStates = parsimony.getPatternStates(nodeNr);
			int [] parentPatternStates = parsimony.getPatternStates(tree.getNode(nodeNr).getParent().getNr());
			for (int siteNr = 0; siteNr < data.getSiteCount(); siteNr++) {
				int k = data.getPatternIndex(siteNr);
				if (nodePatternStates[k] != parentPatternStates[k]) {
					int stateTransition = parentPatternStates[k] * stateCount + nodePatternStates[k];
					state.addMutation0(siteNr, nodeNr, 0.5f, stateTransition);
				}
			}
		}
		
		// set root state
		int [] nodePatternStates = parsimony.getPatternStates(tree.getNodeCount()-1);
		int [] rootState = new int[data.getSiteCount()];
		for (int siteNr = 0; siteNr < data.getSiteCount(); siteNr++) {
			int k = data.getPatternIndex(siteNr);
			rootState[siteNr] = nodePatternStates[k];
		}
		state.setRootState(rootState);
		
		state.calcLengths();
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(stateInput.get());
	}

}
