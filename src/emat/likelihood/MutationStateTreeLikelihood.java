package emat.likelihood;


import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;

@Description("Likelihood for a explicit mutation annotated tree")
public class MutationStateTreeLikelihood extends GenericTreeLikelihood {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree", Validate.REQUIRED);

	private TreeInterface tree;
	private MutationState state;
	private GeneralSubstitutionModel substModel;
	private int stateCount;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		tree = treeInput.get();
		state = stateInput.get();
		substModel = (GeneralSubstitutionModel)((SiteModel.Base)siteModelInput.get()).substModelInput.get();
		substModel.setupRelativeRates();
		substModel.setupRateMatrix();
		
		stateCount = dataInput.get().getMaxStateCount();
		
	}
	
	
	@Override
	public double calculateLogP() {
		logP = 0;
		double [][] rates = substModel.getRateMatrix();
		double [] freqs = substModel.getFrequencies();
		
		int [] mutationCounts = state.getMutationCounts();
		double [] lengths = state.getTotalStateLengths();
		int [] rootStateFreqs = state.getRootStateFreqs();
		
		int k = 0;
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				if (i == j) {
					logP += lengths[i] * rates[i][i];
				} else {
					logP += Math.log(-rates[i][j] / rates[j][j]) * mutationCounts[k++];
				}
			}
		}
		for (int i = 0; i < stateCount; i++) {
			logP += rootStateFreqs[i] * Math.log(freqs[i]);
		}
		
		return logP;
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

}
