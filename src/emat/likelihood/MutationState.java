package emat.likelihood;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import beast.base.util.Randomizer;

@Description("State node that represents a full mutation history on a tree")
public class MutationState extends StateNode {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

	private TreeInterface tree;
	private Alignment data;

	enum EditType {add, delete, replace}
	
	/** list of edits, used to restore the state after a rejected proposal **/
	protected List<Edit> editList = new ArrayList<>();

	
	/** for every site, map nodeNr onto list of mutations **/
	private Map<Integer,List<MutationOnBranch>> [] mutations;
	private int mutationCount;
	private int [] rootState;
	
	int nodeCount;
	int siteCount;
	int stateCount;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		data = dataInput.get();
		
		nodeCount = tree.getNodeCount();
		siteCount = data.getSiteCount();
		stateCount = data.getMaxStateCount();

		mutations = new Map[siteCount];
		for (int i = 0; i < mutations.length; i++) {
			mutations[i] = new HashMap<>();
		}
		
		mutationCount = 0;	
		rootState = new int[siteCount];
	}

	protected void setRootState(int [] rootState) {
		System.arraycopy(rootState, 0, this.rootState, 0, rootState.length);
	}
	
	/** operations on a MutationState: add, delete, replace **/
	
	public void addMutation(int siteNr, int nodeNr, float brancheFraction, int stateTransition) {
		MutationOnBranch mutation = addMutation0(siteNr, nodeNr, brancheFraction, stateTransition);
		editList.add(new Edit(EditType.add, siteNr, nodeNr, mutation));
	}
	
	protected MutationOnBranch addMutation0(int siteNr, int nodeNr, float brancheFraction, int stateTransition) {
		MutationOnBranch mutation = new MutationOnBranch(brancheFraction, stateTransition);
		List<MutationOnBranch> list = mutations[siteNr].get(nodeNr);
		if (list == null) {
			list = new ArrayList<>();
			list.add(mutation);
			mutations[siteNr].put(nodeNr, list);
		} else {
			list.add(mutation);
		}
		mutationCount++;
		return mutation;
	}
	
	public void deleteMutation(int siteNr, int nodeNr, MutationOnBranch mutation) {
		deleteMutation0(siteNr, nodeNr, mutation);
		editList.add(new Edit(EditType.delete, siteNr, nodeNr, mutation));
	}
	
	protected void deleteMutation0(int siteNr, int nodeNr, MutationOnBranch mutation) {
		List<MutationOnBranch> list = mutations[siteNr].get(nodeNr);
		if (list == null) {
			throw new IllegalArgumentException("Could not find mutation");
		}
		list.remove(mutation);
		// maybe save memory?
		// if (list.size() == 0) {
		//  	mutations[siteNr].remove(nodeNr);
		// }
		mutationCount--;
	}

	public void replaceMutation(int siteNr, int nodeNr, MutationOnBranch oldMutation, MutationOnBranch newMutation) {
		replaceMutation0(siteNr, nodeNr, oldMutation, newMutation);
		editList.add(new Edit(EditType.replace, siteNr, nodeNr, oldMutation, newMutation));
	}
	
	protected void replaceMutation0(int siteNr, int nodeNr, MutationOnBranch oldMutation, MutationOnBranch newMutation) {
		List<MutationOnBranch> list = mutations[siteNr].get(nodeNr);
		if (list == null) {
			throw new IllegalArgumentException("Could not find mutation");
		}
		list.set(list.indexOf(oldMutation), newMutation);
	}
	
	
	/**
	 * Uniformly select an existing mutation
	 * @param ids[0] = siteNr, ids[1] = nodeNr
	 * @return MutationOnBranch at site siteNr for node nodeNr
	 */
	public MutationOnBranch getRandomMutation(int [] ids) {
		int i = Randomizer.nextInt(mutationCount);
		int j = 0;
		while (i > 0) {
			for (int nodeNr : mutations[j].keySet()) {
				List<MutationOnBranch> list = mutations[j].get(nodeNr);
				if (i < list.size()) {
					ids[0] = j;
					ids[1] = nodeNr;
					return list.get(i);
				}
			}
		}
		return null;
	}
	
	
	
	@Override
	public void init(PrintStream out) {
		// TODO Auto-generated method stub
	}

	@Override
	public void log(long sample, PrintStream out) {
		// TODO Auto-generated method stub

	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub

	}

	@Override
	public int getDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getArrayValue(int dim) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setEverythingDirty(boolean isDirty) {
		// TODO Auto-generated method stub

	}

	@Override
	public StateNode copy() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void assignTo(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void assignFrom(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void assignFromFragile(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void fromXML(org.w3c.dom.Node node) {
		// TODO Auto-generated method stub

	}

	@Override
	public int scale(double scale) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected void store() {
		editList.clear();
	}

	@Override
	public void restore() {
		for (Edit e : editList) {
			e.undo(this);
		}
		editList.clear();
	}

	
	
	
	private double [][] stateLengths;
	private int [] mutationCounts;
	
	protected void calcLengths() {
		if (stateLengths == null) {
			stateLengths = new double[siteCount][stateCount];
			mutationCounts = new int[stateCount * stateCount];
		}
		collectStateLengths(tree.getRoot());
	}

	private int [] collectStateLengths(Node node) {
		int [] states = null;
		if (node.isLeaf()) {
            states = new int[siteCount];
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (int i = 0; i < siteCount; i++) {
            	int k = data.getPatternIndex(i);
                int code = data.getPattern(taxonIndex, k);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                if (statesForCode.length == 1) {
                    states[i] = statesForCode[0];
                } else {
                    states[i] = code; // Causes ambiguous states to be ignored.
                }
            }
		} else {
			int [] state0 = collectStateLengths(node.getLeft());
			states = collectStateLengths(node.getRight());
			// sanity check: make sure state and state0 are identical
		}
		
		// update stateLengths and apply mutations
		double len = node.getLength();
		int nodeNr = node.getNr();
		for (int i = 0; i < siteCount; i++) {
			Map<Integer, List<MutationOnBranch>> map = mutations[i];
			List<MutationOnBranch> list = map.get(nodeNr);
			if (list == null | list.size() == 0) {
				stateLengths[i][states[i]] += len;
			} else {
				Collections.sort(list);
				int k = states[i];
				double offset = 0;
				for (MutationOnBranch m : list) {
					stateLengths[i][k] += len * (m.brancheFraction - offset);
					offset = m.brancheFraction;
					k = m.stateTransition / stateCount;
					mutationCounts[m.stateTransition]++;
				}
				stateLengths[i][k] += len * (1.0 - offset);
				states[i] = k;
			}
		}
		
		return states;
	}

	/**
    *
    * @param taxon the taxon name as a string
    * @param data the alignment
    * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
    *         or -1 if the taxon is not in the alignment.
    */
   private int getTaxonIndex(String taxon, Alignment data) {
       int taxonIndex = data.getTaxonIndex(taxon);
       if (taxonIndex == -1) {
       	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
               taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
           }
           if (taxonIndex == -1) {
           	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
           }
       }
       return taxonIndex;
	}
}
