package emat.likelihood;



import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;
import beast.base.inference.util.InputUtil;
import emat.substitutionmodel.EmatSubstitutionModel;

@Description("Likelihood for a explicit mutation annotated tree")
public class MutationStateTreeLikelihood extends GenericTreeLikelihood {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree", Validate.REQUIRED);
	final public Input<EditList> editListInput = new Input<>("editList", "list of edit actions on editable tree", Validate.REQUIRED);
	
	final public Input<EmatSubstitutionModel> substModelInput =
            new Input<>("substModel", "emat substitution model along branches in the beast.tree", null, Validate.REQUIRED);

	private TreeInterface tree;
	private MutationState state;
	// private GeneralSubstitutionModel substModel;
	private BranchRateModel clockModel;
	private int stateCount;
	
	private List<Edit> editList;
	private double deltaLogP;
	
	private double [][] branchLogP;
	private int [] currentBranchLogPInidicator;
	private boolean needsUpdate = true, neededUpdate = false;
	
	final static boolean debug = false;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		tree = treeInput.get();
		state = stateInput.get();
//		substModel = (GeneralSubstitutionModel)((SiteModel.Base)siteModelInput.get()).substModelInput.get();
//		substModel.setupRelativeRates();
//		substModel.setupRateMatrix();
		
		clockModel = branchRateModelInput.get();
		
		stateCount = dataInput.get().getMaxStateCount();
		editList = new ArrayList<>();
		
		deltaLogP = 0;

		int nodeCount = tree.getNodeCount();
		branchLogP = new double[2][nodeCount];
		currentBranchLogPInidicator = new int[nodeCount];
	}
	
	
	@Override
	public double calculateLogP() {
//System.err.println("\n\ncalculateLogP");		
		if (Double.isNaN(logP) || needsUpdate) {
			initLogP();
			needsUpdate = false;
			neededUpdate = true;
		} else {
			logP += deltaLogP;
		}
		if (debug) {
			double logP0 = logP;
			initLogP();
			if (Math.abs(logP - logP0) > 1e-6) {
				int h = 3;
				h++;
			}
		}
		return logP;
	}
	
	void initLogP() {
		logP = 0;
		
		// contribution of each branch in the tree
		for (int nodeNr = 0; nodeNr < tree.getNodeCount() - 1; nodeNr++) {
			double branchLogP = calculateLogPForBranch(nodeNr);
			if (needsUpdate) {
				flipBranchLogPInidicator(nodeNr);
			}
			if (debug && !needsUpdate && Math.abs(this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] - branchLogP) > 1e-6) {
				int h = 3;
				h++;
			}
			this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] = branchLogP;
			logP += branchLogP;
		}
		
		// contribution of root state
		double rootBranchLogP = 0;
		double [] freqs = substModelInput.get().getFrequencies();
		int [] rootStateFreqs = state.getRootStateFreqs();
		for (int i = 0; i < stateCount; i++) {
			rootBranchLogP += rootStateFreqs[i] * Math.log(freqs[i]);
		}
		int nodeNr = tree.getRoot().getNr();
		if (needsUpdate) {
			currentBranchLogPInidicator[nodeNr] = 1 - currentBranchLogPInidicator[nodeNr];
		}
		if (debug && !needsUpdate && Math.abs(this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] - rootBranchLogP) > 1e-6) {
			int h = 3;
			h++;
		}
		this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] = rootBranchLogP;
		logP += rootBranchLogP;
	}
	
	private void flipBranchLogPInidicator(int nodeNr) {
		currentBranchLogPInidicator[nodeNr] = 1 - currentBranchLogPInidicator[nodeNr];
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

	
	@Override
	public void store() {
		super.store();
		
		deltaLogP = 0;
		editList.clear();
		neededUpdate = false;
	}
	
	@Override
	public void restore() {
		deltaLogP = 0;
		for (Edit e : editList) {
			e.undo(this);
		}
		super.restore();
		if (neededUpdate) {
			for (int nodeNr = 0; nodeNr < tree.getNodeCount(); nodeNr++) {
				currentBranchLogPInidicator[nodeNr] = 1 - currentBranchLogPInidicator[nodeNr];
			}
			neededUpdate = false;
		}
		editList.clear();
	}
	
	@Override
	protected boolean requiresRecalculation() {
		deltaLogP = 0;
		editList.clear();
		editList.addAll(editListInput.get().list);
		for (Edit e : editList) {
			e.apply(this);
		}

		neededUpdate = false;
		if (InputUtil.isDirty(branchRateModelInput) || InputUtil.isDirty(siteModelInput)) {
			deltaLogP = 0;
			needsUpdate = true;
		}
		
		return super.requiresRecalculation();
	}

	/**
	 * For each site starting at state i ending in state j over branch length T
	 * L(path | R, i, j, T)     = [ Π_{k=0 to N_actual-1} ( R[s_k, s_{k+1}] * exp(R[s_k, s_k] * (t_{k+1} - t_k)) ) ] 
	 * 							* exp(R[s_{N_actual}, s_{N_actual}] * (T - t_{N_actual}))
	 * so 
	 * log L(path | R, i, j, T) = [ Σ_{k=0 to N_actual-1} ( log(R[s_k, s_{k+1}]) + R[s_k, s_k] * (t_{k+1} - t_k) ) ] 
	 * 				            + R[s_{N_actual}, s_{N_actual}] * (T - t_{N_actual})
	 * where     
	 * o N_actual: The number of actual state changes (not fictitious jumps).
     * o t_0 = 0, t_1, t_2, ..., t_{N_actual}: The times of the actual state changes.
     * o s_0 = i, s_1, s_2, ..., s_{N_actual} = j: The sequence of states visited, 
     *   where s_k is the state occupied in [t_k, t_{k+1}).
	 * o The final state s_{N_actual} is held until T, and s_{N_actual} must be j.
	 * 
	 * 
	 * log L(R | all paths) = Σ_{state pairs (u,v), u≠v} N_uv * log(R[u,v]) + Σ_{state u} T_u * R[u,u]
	 * o N_uv: The total number of observed direct transitions from state u to state v across all paths in your dataset.
	 * o T_u: The total time spent in state u across all paths in your dataset.
	 * o These N_uv and T_u are "sufficient statistics" for the CTMC.
	 */

	private double calculateLogPForBranch(int nodeNr) {
		double [][] rates = substModelInput.get().getRateMatrix();
		
		double length = tree.getNode(nodeNr).getLength();
		double clockRate = clockModel.getRateForBranch(tree.getNode(nodeNr));
		
		double branchLogP = 0;
		int [] branchMutationCount = state.branchMutationCount[nodeNr];
		double [] branchStateLength = state.branchStateLength[nodeNr];
		
//System.err.println(nodeNr + Arrays.toString(branchMutationCount) + Arrays.toString(branchStateLength));		
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				if (i == j) {
					branchLogP += rates[i][i] * clockRate * branchStateLength[i] * length;
				} else {
					branchLogP += Math.log(rates[i][j] * clockRate) * branchMutationCount[i*stateCount+j] ; 
				}
			}
		}
		
		return branchLogP;
	}
	
	public void moveBranchFraction(int nodeNr, double newBranchFractionx) {
		double branchLogP = calculateLogPForBranch(nodeNr);
		deltaLogP = -this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr];
		flipBranchLogPInidicator(nodeNr);
		this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] = branchLogP;
		deltaLogP += branchLogP;
	}



	public void moveNode(int nodeNr, double d) {
		Node node = tree.getNode(nodeNr);
		if (!node.isRoot()) {
			recalcBranchContribution(nodeNr);
		}

		final int left = node.getLeft().getNr();
		recalcBranchContribution(left);

		final int right = node.getRight().getNr();
		recalcBranchContribution(right);		
	}


	public void undoMoveNode(int nodeNr) {
		Node node = tree.getNode(nodeNr);
		if (!node.isRoot()) {	
			flipBranchLogPInidicator(nodeNr);
		}
		final int left = node.getLeft().getNr();
		flipBranchLogPInidicator(left);
		final int right = node.getRight().getNr();
		flipBranchLogPInidicator(right);
	}


	public void recalcBranchContribution(int nodeNr) {
		double branchLogP = calculateLogPForBranch(nodeNr);
		deltaLogP += -this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr];
		flipBranchLogPInidicator(nodeNr);
		this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] = branchLogP;
		deltaLogP += branchLogP;
	}

	public void undoMoveBranchFraction(int nodeNr) {
		flipBranchLogPInidicator(nodeNr);
	}

	public void undoBranchContribution(int nodeNr) {
		flipBranchLogPInidicator(nodeNr);
	}
	
//	public GeneralSubstitutionModel getSubstModel() {
//		return substModel;
//	}
}
