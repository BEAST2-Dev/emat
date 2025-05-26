package emat.likelihood;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.BEASTInterface;
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

	/** list of edits, used to restore the state after a rejected proposal **/
	protected List<Edit> editList;

	
	/** for every branch, track list of mutations **/
	private List<MutationOnBranch> [] branchMutations;
	
	private int mutationCount;
	private int [][] nodeSequence;
	
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

		branchMutations = new List[nodeCount];
		for (int i = 0; i < nodeCount; i++) {
			branchMutations[i] = new ArrayList<>();
		}
		
		mutationCount = 0;	
		nodeSequence = new int[nodeCount][];
		
		
		for (BEASTInterface o : getOutputs()) {
			if (o instanceof EditList list) {
				editList = list.list;
			}
		}
		if (editList == null) {
			throw new RuntimeException("Expected an " + EditList.class.getName() + " to be present");
		}
	}

	protected void setNodeSequence(int nodeNr, int [] nodeSequence) {
		this.nodeSequence[nodeNr] = nodeSequence;
	}
	
	public int getCharAt(int nodeNr, int siteNr) {
		return nodeSequence[nodeNr][siteNr];
	}

	/** operations on a MutationState: add, delete, replace **/
	public void addMutation(int siteNr, int nodeNr, float brancheFraction, int stateTransition) {
		startEditing(null);
		MutationOnBranch mutation = addMutation0(siteNr, nodeNr, brancheFraction, stateTransition);
		editList.add(new Edit(EditType.addMutation, siteNr, nodeNr, mutation));
	}
	
	protected MutationOnBranch addMutation0(int siteNr, int nodeNr, float brancheFraction, int stateTransition) {
		MutationOnBranch mutation = new MutationOnBranch(nodeNr, brancheFraction, stateTransition, siteNr);
		List<MutationOnBranch> list = branchMutations[nodeNr];
		int i = 0;
		while (i < list.size() && brancheFraction > list.get(i).brancheFraction) {
			i++;
		}
		list.add(i, mutation);
		mutationCount++;
		return mutation;
	}
	
	public void deleteMutation(int siteNr, int nodeNr, MutationOnBranch mutation) {
		startEditing(null);
		deleteMutation0(siteNr, nodeNr, mutation);
		editList.add(new Edit(EditType.deleteMutation, siteNr, nodeNr, mutation));
	}
	
	protected void deleteMutation0(int siteNr, int nodeNr, MutationOnBranch mutation) {
		List<MutationOnBranch> list = branchMutations[nodeNr];
		if (list == null) {
			throw new IllegalArgumentException("Could not find mutation");
		}
		list.remove(mutation);
		mutationCount--;
	}

	public void replaceMutation(int siteNr, int nodeNr, MutationOnBranch oldMutation, MutationOnBranch newMutation) {
		startEditing(null);
		replaceMutation0(siteNr, nodeNr, oldMutation, newMutation);
		editList.add(new Edit(EditType.replaceMutation, siteNr, nodeNr, oldMutation, newMutation));
	}
	
	protected void replaceMutation0(int siteNr, int nodeNr, MutationOnBranch oldMutation, MutationOnBranch newMutation) {
		List<MutationOnBranch> list = branchMutations[nodeNr];
		if (list == null) {
			throw new IllegalArgumentException("Could not find mutation");
		}
		list.set(list.indexOf(oldMutation), newMutation);
	}
	
	public void moveBranchFraction(int siteNr, int nodeNr, MutationOnBranch mutation, float oldBranchFraction, float newBranchFraction) {
		startEditing(null);
		moveBranchFraction0(mutation, siteNr, nodeNr, newBranchFraction);
		editList.add(new Edit(EditType.moveBranchFraction, siteNr, nodeNr, mutation, oldBranchFraction, newBranchFraction));
	}
	
	protected void moveBranchFraction0(MutationOnBranch mutation, int siteNr, int nodeNr, float branchFraction) {
		mutation.brancheFraction = branchFraction;

		// ensure list remains in order of branch fractions
		Collections.sort(branchMutations[nodeNr]);
	}
	
	
	/**
	 * Uniformly select an existing mutation
	 * @param ids[0] = siteNr, ids[1] = nodeNr
	 * @return MutationOnBranch at site siteNr for node nodeNr
	 */
	public MutationOnBranch getRandomMutation(int [] ids) {
		int i = Randomizer.nextInt(mutationCount);
		int j = 0;
		while (i >= 0) {
			List<MutationOnBranch> list = branchMutations[j];
			if (i < list.size()) {
				ids[0] = list.get(i).siteNr;
				ids[1] = j;
				return list.get(i);
			} else {
				i -= list.size();
			}
			j++;
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


	/** number of sites having a particular root state **/
	private int [] rootStateFreqs;

	public int[] getRootStateFreqs() {
		return rootStateFreqs;
	}

	/** 
	 * recalculate stateLengths,  totalStateLengths and mutationCounts from scratch
	 * instead of incrementally
	 **/
	protected void calcLengths() {
		if (rootStateFreqs == null) {
			rootStateFreqs = new int[stateCount];
		} else {
			Arrays.fill(rootStateFreqs, 0);
		}

		int [] rootStates = collectStateLengths(tree.getRoot());
		for (int i : rootStates) {
			rootStateFreqs[i]++;
		}
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
		
		// apply mutations
		int nodeNr = node.getNr();
		List<MutationOnBranch> list = branchMutations[nodeNr];
		Collections.sort(list);
		for (MutationOnBranch m : list) {
			states[m.siteNr] = m.stateTransition/stateCount;
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

	public List<MutationOnBranch> getMutationList(int nodeNr) {
		return branchMutations[nodeNr];
	}



}
