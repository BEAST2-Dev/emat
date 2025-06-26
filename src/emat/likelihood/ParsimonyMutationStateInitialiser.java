package emat.likelihood;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beastlabs.parsimony.ParsimonyCriterion;
import beastlabs.parsimony.FitchParsimony;
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
		//parsimony.getScore((Tree) tree);
		
		// set root state
		int [] nodePatternStates = parsimony.getStates((Tree)tree, tree.getRoot());
		int [] rootState = new int[data.getSiteCount()];
		for (int siteNr = 0; siteNr < data.getSiteCount(); siteNr++) {
			int k = data.getPatternIndex(siteNr);
			rootState[siteNr] = nodePatternStates[k];
			if (rootState[siteNr] >= stateCount) {
				// replace ambiguous character with majority
				// should never get here if parsimony assignment is correct
				rootState[siteNr] = majorityState(data.getPattern(k), stateCount);
				
			}
		}
		
		
		
		state.setNodeSequence(tree.getNodeCount()-1, rootState);

		// add mutations
		traverse(tree.getRoot(), parsimony, tree, data, state);
		
		state.reinitialise();
	}

	private int majorityState(int[] pattern, int stateCount) {
		int [] count = new int[stateCount];
		for (int i : pattern) {
			if (i < stateCount) {
				count[i]++;
			}
		}
		int iMax = 0;
		int maxCount = count[0];
		for (int i = 1; i < stateCount; i++) {
			if (count[i] > maxCount) {
				iMax = i;
				maxCount = count[i];
			}
		}
		return iMax;
	}

	private void traverse(Node node, ParsimonyCriterion parsimony, TreeInterface tree, Alignment data, MutationState state) {
		if (!node.isRoot()) {
			int nodeNr = node.getNr();
			int stateCount = data.getMaxStateCount();
			
			int [] nodePatternStates = parsimony.getStates((Tree)tree, tree.getNode(nodeNr));
			int [] parentSites = state.getNodeSequence(node.getParent().getNr()); 
			int [] nodeSequence = new int[data.getSiteCount()];
			for (int siteNr = 0; siteNr < data.getSiteCount(); siteNr++) {
				int k = data.getPatternIndex(siteNr);
				nodeSequence[siteNr] = nodePatternStates[k];
				if (nodePatternStates[k] != parentSites[siteNr]) {
					if (nodePatternStates[k] < stateCount) {
						state.addMutation(siteNr, nodeNr, 0.5f, nodePatternStates[k], parentSites[siteNr]);
					} else {
						// fill in missing data, gaps and ambiguous characters by parent state
						if (parentSites[siteNr] >= stateCount) {
							int h = 3;
							h++;
						}
						nodeSequence[siteNr] = parentSites[siteNr];
					}
				}
			}
			state.setNodeSequence(nodeNr, nodeSequence);
		}
		
		if (!node.isLeaf()) {
			traverse(node.getLeft(), parsimony, tree, data, state);
			traverse(node.getRight(), parsimony, tree, data, state);
		}
		
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(stateInput.get());
	}

}
