package emat.likelihood;



import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;

@Description("Likelihood for a explicit mutation annotated tree")
public class MutationStateTreeLikelihood extends GenericTreeLikelihood {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree", Validate.REQUIRED);
	final public Input<EditList> editListInput = new Input<>("editList", "list of edit actions on editable tree", Validate.REQUIRED);

	private TreeInterface tree;
	private MutationState state;
	private GeneralSubstitutionModel substModel;
	private int stateCount;
	
	private List<Edit> editList;
	private double deltaLogP;
	
	private double [][] branchLogP;
	private int [] currentBranchLogPInidicator;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		tree = treeInput.get();
		state = stateInput.get();
		substModel = (GeneralSubstitutionModel)((SiteModel.Base)siteModelInput.get()).substModelInput.get();
		substModel.setupRelativeRates();
		substModel.setupRateMatrix();
		
		stateCount = dataInput.get().getMaxStateCount();
		editList = new ArrayList<>();
		
		deltaLogP = 0;

		int nodeCount = tree.getNodeCount();
		branchLogP = new double[2][nodeCount];
		currentBranchLogPInidicator = new int[nodeCount];
	}
	
	
	@Override
	public double calculateLogP() {
		if (Double.isNaN(logP)) {
			initLogP();
		} else {
			logP += deltaLogP;
		}
		return logP;
	}
	
	void initLogP() {
		logP = 0;
		
		// contribution of each branch in the tree
		for (int nodeNr = 0; nodeNr < tree.getNodeCount() - 1; nodeNr++) {
			double branchLogP = calculateLogPForBranch(nodeNr);
			this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] = branchLogP;
			logP += branchLogP;
		}
		
		// contribution of root state
		double rootBranchLogP = 0;
		double [] freqs = substModel.getFrequencies();
		int [] rootStateFreqs = state.getRootStateFreqs();
		for (int i = 0; i < stateCount; i++) {
			rootBranchLogP += rootStateFreqs[i] * Math.log(freqs[i]);
		}
		int nodeNr = tree.getRoot().getNr();
		this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] = rootBranchLogP;
		logP += rootBranchLogP;
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
	}
	
	@Override
	public void restore() {
		deltaLogP = 0;
		for (Edit e : editList) {
			e.undo(this);
		}
		super.restore();
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

		return super.requiresRecalculation();
	}


	private double calculateLogPForBranch(int nodeNr) {
		double [][] rates = substModel.getRateMatrix();
		int siteCount = state.siteCount;
		double notMutatingRate = -1;

		List<MutationOnBranch> list = state.getMutationList(nodeNr);
		Node node = tree.getNode(nodeNr);
		double start = node.getHeight();
		double end = start;
		double totalLength = node.getParent().getHeight() - start;
		// process all mutations on this branch
		double branchLogP = 0;
		for (MutationOnBranch mutation : list) {
			end = mutation.brancheFraction * totalLength;
			double len = end - start;
			// contribution of branch not having a mutation
			double p = 1.0 - Math.exp(len * notMutatingRate);
			if (p > 0) {
				branchLogP += siteCount * Math.log(p);
			}
			// contribution of mutation
			int i = mutation.stateTransition % 4;
			int j = mutation.stateTransition / 4;
			branchLogP += Math.log(-rates[i][j] / rates[j][j]);
			start = end;
		}
		
		// contribution of branch not having a mutation
		double len = node.getParent().getHeight() - end;
		double p = 1.0 - Math.exp(len * notMutatingRate);
		if (p > 0) {
			branchLogP += siteCount * Math.log(p);
		}
		return branchLogP;
	}
	
	public void moveBranchFraction(int nodeNr, double newBranchFractionx) {
		double branchLogP = calculateLogPForBranch(nodeNr);
		deltaLogP = -this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr];
		currentBranchLogPInidicator[nodeNr] = 1-currentBranchLogPInidicator[nodeNr];
		this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] = branchLogP;
		deltaLogP += branchLogP;
	}


	public void undoMoveBranchFraction(int nodeNr) {
		currentBranchLogPInidicator[nodeNr] = 1-currentBranchLogPInidicator[nodeNr];
	}


	public void moveNode(int nodeNr, double d) {
		Node node = tree.getNode(nodeNr);
		if (!node.isRoot()) {
			double branchLogP = calculateLogPForBranch(nodeNr);
			deltaLogP = -this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr];
			currentBranchLogPInidicator[nodeNr] = 1-currentBranchLogPInidicator[nodeNr];
			this.branchLogP[currentBranchLogPInidicator[nodeNr]][nodeNr] = branchLogP;
			deltaLogP += branchLogP;
		}

		final int left = node.getLeft().getNr();
		double branchLogP = calculateLogPForBranch(left);

		deltaLogP += -this.branchLogP[currentBranchLogPInidicator[left]][left];
		currentBranchLogPInidicator[left] = 1-currentBranchLogPInidicator[left];
		this.branchLogP[currentBranchLogPInidicator[left]][left] = branchLogP;
		deltaLogP += branchLogP;

		final int right = node.getRight().getNr();
		branchLogP = calculateLogPForBranch(right);

		deltaLogP += -this.branchLogP[currentBranchLogPInidicator[right]][right];
		currentBranchLogPInidicator[right] = 1-currentBranchLogPInidicator[right];
		this.branchLogP[currentBranchLogPInidicator[right]][right] = branchLogP;
		deltaLogP += branchLogP;
		
	}


	public void undoMoveNode(int nodeNr) {
		Node node = tree.getNode(nodeNr);
		if (!node.isRoot()) {		
			currentBranchLogPInidicator[nodeNr] = 1-currentBranchLogPInidicator[nodeNr];
		}
		final int left = node.getLeft().getNr();
		currentBranchLogPInidicator[left] = 1-currentBranchLogPInidicator[left];
		final int right = node.getRight().getNr();
		currentBranchLogPInidicator[right] = 1-currentBranchLogPInidicator[right];
	}
}
