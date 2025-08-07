package emat.likelihood;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import beast.base.util.Randomizer;

@Description("State node that represents a full mutation history on a tree")
public class MutationState extends StateNode {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

    final static boolean debug = true;
    
	private TreeInterface tree;
	private Alignment data;

	/** list of edits, used to restore the state after a rejected proposal **/
	protected EditList editList;

	
	/** for every branch, track list of mutations **/
	private List<MutationOnBranch> [] branchMutations;
	
	private int mutationCount;
	/** nodeSequence[currentNodeSequence][nodeNr][siteNr] **/
	private int [][][] nodeSequence;
	private int [] currentNodeSequence;

	/** 
	 * track for every node how many mutations there are:
	 * branchMutationCount[nodeNr][stateTransition]
	 */
	protected int [][] branchMutationCount;
	/** 
	 * track for every node which fraction of the time is spent in a state
	 * to get the full length, multiply by branch length
	 * branchStateLength[nodeNr][state]
	 */
	protected double [][] branchStateLength;

	private int nodeCount;
	private int siteCount;
	private int stateCount, stateCountSquared;
	
	public int getStateCount() {
		return stateCount;
	}
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		data = dataInput.get();
		
		nodeCount = tree.getNodeCount();
		siteCount = data.getSiteCount();
		stateCount = data.getMaxStateCount();
		stateCountSquared = stateCount * stateCount;

		branchMutations = new List[nodeCount];
		for (int i = 0; i < nodeCount; i++) {
			branchMutations[i] = new ArrayList<>();
		}
		
		mutationCount = 0;	
		nodeSequence = new int[2][nodeCount][siteCount];
		currentNodeSequence = new int[nodeCount];
		
		for (BEASTInterface o : getOutputs()) {
			if (o instanceof EditList list) {
				editList = list;
			}
		}
		if (editList == null) {
			throw new RuntimeException("Expected an " + EditList.class.getName() + " to be present");
		}
	}

	protected void setNodeSequence(int nodeNr, int [] nodeSequence) {
		System.arraycopy(nodeSequence, 0, this.nodeSequence[currentNodeSequence[nodeNr]][nodeNr], 0, nodeSequence.length);
	}
	
	public int [] getNodeSequence(int nodeNr) {
		return this.nodeSequence[currentNodeSequence[nodeNr]][nodeNr];
	}
	
	public int getCharAt(int nodeNr, int siteNr) {
		return nodeSequence[currentNodeSequence[nodeNr]][nodeNr][siteNr];
	}

	/** operations on a MutationState: add, delete, replace **/
	public void addMutation(int siteNr, int nodeNr, double brancheFraction, int fromState, int toState) {
		MutationOnBranch mutation = new MutationOnBranch(nodeNr, brancheFraction, fromState, toState, siteNr);
		List<MutationOnBranch> list = branchMutations[nodeNr];
		int i = 0;
		while (i < list.size() && brancheFraction > list.get(i).brancheFraction) {
			i++;
		}
		list.add(i, mutation);
		mutationCount++;
	}
	
	public void moveBranchFraction(int siteNr, int nodeNr, MutationOnBranch mutation, double oldBranchFraction, double newBranchFraction) {
		startEditing(null);
		moveBranchFraction0(mutation, siteNr, nodeNr, newBranchFraction);
		editList.add(new Edit(EditType.moveBranchFraction, siteNr, nodeNr, mutation, oldBranchFraction, newBranchFraction));
	}
	
	protected void moveBranchFraction0(MutationOnBranch mutation, int siteNr, int nodeNr, double branchFraction) {
		
		mutation.brancheFraction = branchFraction;

		// ensure list remains in order of branch fractions
		Collections.sort(branchMutations[nodeNr]);

		// TODO: update branchStateLengths
		Log.warning("branchStateLengths are not updated -- do not use this method");
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
		hasStartedEditing = false;
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
		// editList.clear();
	}

	@Override
	public void restore() {
//		for (Edit e : editList) {
//			e.undo(this);
//		}
//		editList.clear();
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
	protected void reinitialise() {
		if (rootStateFreqs == null) {
			rootStateFreqs = new int[stateCount];
			branchMutationCount = new int[nodeCount][stateCountSquared];
			branchStateLength = new double[nodeCount][stateCount];
		} else {
			Arrays.fill(rootStateFreqs, 0);
			for (int i = 0; i < nodeCount; i++) {
				Arrays.fill(branchMutationCount[i], 0);
				Arrays.fill(branchStateLength[i], 0.0);
			}
		}

		int [] rootStates = collectStateLengths(tree.getRoot());
		for (int i : rootStates) {
			rootStateFreqs[i]++;
		}
	}

	private int [] collectStateLengths(Node node) {
		int [] states = null;
		if (node.isLeaf()) {
            states = getNodeSequence(node.getNr()).clone();
//            int taxonIndex = getTaxonIndex(node.getID(), data);
//            for (int i = 0; i < siteCount; i++) {
//            	int k = data.getPatternIndex(i);
//                int code = data.getPattern(taxonIndex, k);
//                int[] statesForCode = data.getDataType().getStatesForCode(code);
//                if (statesForCode.length == 1) {
//                    states[i] = statesForCode[0];
//                } else {
//                    states[i] = code; // Causes ambiguous states to be ignored.
//                }
//            }
		} else {
			int [] state0 = collectStateLengths(node.getLeft());
			states = collectStateLengths(node.getRight());
			if (debug) {
				// sanity check: make sure state and state0 are identical
				for (int i = 0; i < states.length; i++) {
					if (state0[i] != states[i] && state0[i] < stateCount && states[i] < stateCount) {
						throw new RuntimeException("Incompatible reconstruction of internal node states found");
					}
				}
				// sanity check: make sure state and nodeSequence are identical
				int nodeNr = node.getNr();
				for (int i = 0; i < states.length; i++) {
					if (nodeSequence[currentNodeSequence[nodeNr]][nodeNr][i] != states[i] && states[i] < stateCount) {
						states = collectStateLengths(node.getRight());
						throw new RuntimeException("Node sequences and reconstruction of internal node states are not compatible");
					}
				}
			}
			for (int i = 0; i < states.length; i++) {
				states[i] = Math.min(states[i], state0[i]);
			}
		}
		
		
		int [] stateCounts = new int[stateCount+1];
        for (int i : states) {
        	stateCounts[i]++;
        }
        
		// apply mutations
		int nodeNr = node.getNr();
		List<MutationOnBranch> list = branchMutations[nodeNr];
		Collections.sort(list);
		double prev = 0;
		for (MutationOnBranch m : list) {
			states[m.siteNr] = m.getToState();
			try {
				branchMutationCount[nodeNr][m.getFromState() * stateCount + m.getToState()]++;
			} catch (ArrayIndexOutOfBoundsException e) {
				int h = 3;
				h++;
			}
	        for (int i = 0; i < stateCount; i++) {
	        	branchStateLength[nodeNr][i] += (m.brancheFraction - prev) * stateCounts[i];
	        }
	        // going backward in time, so going toState => fromState
	        stateCounts[m.getFromState()]--;
	        stateCounts[m.getToState()]++;
	        prev = m.brancheFraction;
		}

        for (int i = 0; i < stateCount; i++) {
        	branchStateLength[nodeNr][i] += (1.0 - prev) * stateCounts[i];
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

	public void setBranchMutations(int nodeNr, List<MutationOnBranch> mutations) {
		// next commented out lines are incorrect
		// since branch lengths may have changed, so contribution to likelihood needs to be updated.
		// TODO: split recalculation due to mutations and those due to branch length changes more efficiently
//		if (branchMutations[nodeNr].size() == 0 && mutations.size() == 0) {
//			// no mutations, nothing has changed
//			return;
//		}
		startEditing(null);
		editList.add(new Edit(EditType.resample, nodeNr, branchMutations[nodeNr]));
		applyMutations(nodeNr, mutations);
	}

	public void setBranchMutationsAfterSlide(int nodeNr, List<MutationOnBranch> mutations) {
		startEditing(null);
		editList.add(new Edit(EditType.resampleAfterSlide, nodeNr, branchMutations[nodeNr]));
		applyMutations(nodeNr, mutations);
	}

	public void applyMutations(int nodeNr, List<MutationOnBranch> mutations) {
		// TODO: cache these calculations?
		
		branchMutations[nodeNr] = mutations;

		// update branchStateLength and branchMutationCount
		int [] states = getNodeSequence(nodeNr).clone();

		Arrays.fill(branchStateLength[nodeNr], 0.0);
		Arrays.fill(branchMutationCount[nodeNr], 0);

		// stateCounts[i] = number of sites in state i, which changes through the branch
		int [] stateCounts = new int[stateCount];
        for (int i : states) {
        	stateCounts[i]++;
        }
		
		// apply mutations
		Collections.sort(mutations);
		double prev = 0;
		for (MutationOnBranch m : mutations) {
			states[m.siteNr] = m.getToState();
        	branchMutationCount[nodeNr][m.getFromState() * stateCount + m.getToState()]++;
        	double delta = (m.brancheFraction - prev);
	        for (int i = 0; i < stateCount; i++) {
	        	branchStateLength[nodeNr][i] += delta * stateCounts[i];
	        }
	        // going backward in time, so going toState => fromState
	        stateCounts[m.getFromState()]--;
	        stateCounts[m.getToState()]++;
	        prev = m.brancheFraction;
		}

        for (int i = 0; i < stateCount; i++) {
        	branchStateLength[nodeNr][i] += (1.0 - prev) * stateCounts[i];
        }
        
        if (debug) {
    		int [] parentStates = getNodeSequence(tree.getNode(nodeNr).getParent().getNr());

            for (int i = 0; i < siteCount; i++) {
            	if (states[i] != parentStates[i]) {
            		throw new RuntimeException("states differ from parentstates after taking mutations on branch in account");
            	}
            }
        }
	}

	public int [] getNodeSequenceForUpdate(int nodeNr) {
		startEditing(null);
		editList.add(new Edit(EditType.setsequence, nodeNr, 0.0, 0.0));
		int k = currentNodeSequence[nodeNr];
		System.arraycopy(nodeSequence[k][nodeNr], 0, nodeSequence[1-k][nodeNr], 0, siteCount);
		return nodeSequence[k][nodeNr];
	}

	public void flipCurrentNodeSequence(int nodeNr) {
		currentNodeSequence[nodeNr] = 1 - currentNodeSequence[nodeNr];
	}

}
