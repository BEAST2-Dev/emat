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

	private TreeInterface tree;
	private MutationState state;
	private GeneralSubstitutionModel substModel;
	private int stateCount;
	
	private List<Edit> editList;
	private double deltaLogP;
	
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
		double [][] rates = substModel.getRateMatrix();
		double [] freqs = substModel.getFrequencies();
		
		int [] rootStateFreqs = state.getRootStateFreqs();
		
		int siteCount = state.siteCount;
		double notMutatingRate = -1;
		
		// contribution of each branch in the tree
		for (int nodeNr = 0; nodeNr < tree.getNodeCount() - 1; nodeNr++) {
			// contribution of each part of the branch till a mutation is encountered
			Node node = tree.getNode(nodeNr);
			List<MutationOnBranch> list = state.getMutationList(nodeNr);
			double start = node.getHeight();
			double end = start;
			double totalLength = end - start;
			// process all mutations on this branch
			for (MutationOnBranch mutation : list) {
				end = mutation.brancheFraction * totalLength;
				double len = end - start;
				// contribution of branch not having a mutation
				logP += Math.log(1.0 - Math.exp(len * siteCount * notMutatingRate));
				// contribution of mutation
				int i = mutation.stateTransition % 4;
				int j = mutation.stateTransition / 4;
				logP += Math.log(-rates[i][j] / rates[j][j]);
				start = end;
			}
			
			// contribution of branch not having a mutation
			double len = node.getParent().getHeight() - end;
			logP += Math.log(1.0 - Math.exp(len * siteCount * notMutatingRate));
		}
		
		// contribution of root state
		for (int i = 0; i < stateCount; i++) {
			logP += rootStateFreqs[i] * Math.log(freqs[i]);
		}
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
	protected boolean requiresRecalculation() {
		MutationState state = stateInput.get();
		editList.clear();
		editList.addAll(state.editList);
		for (Edit e : editList) {
			e.apply(this);
		}

		return super.requiresRecalculation();
	}
}
