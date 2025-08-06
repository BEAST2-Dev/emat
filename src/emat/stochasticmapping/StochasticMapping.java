package emat.stochasticmapping;

import java.util.ArrayList;
import java.util.List;

import emat.substitutionmodel.EmatSubstitutionModel;

public interface StochasticMapping {

    default public List<TimeStateInterval> generatePath(
    		EmatSubstitutionModel substModel,
            int [] startState,
            int [] endState,
            double totalTime) {
    	setSubstModel(substModel);
    	return generatePath(startState, endState, totalTime);
    }

    // public void setRatematrix(double [][] rateMatrixR);
    public void setSubstModel(EmatSubstitutionModel substModel);

    public List<TimeStateInterval> generatePath(
    		int site,
            int startState,
            int endState,
            double totalTime);
    
    default public List<TimeStateInterval> generatePath(
            int [] startState,
            int [] endState,
            double totalTime) {
    	List<TimeStateInterval> mutations = new ArrayList();
    	for (int i = 0; i < startState.length; i++) {
    		List<TimeStateInterval> stateMutations = generatePath(i, startState[i], endState[i], totalTime);
    		mutations.addAll(stateMutations);
    	}
    	return mutations;
    }
}
