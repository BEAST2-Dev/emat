package emat.stochasticmapping;

import java.util.ArrayList;
import java.util.List;

public interface StochasticMapping {

    default public List<TimeStateInterval> generatePath(
            double[][] rateMatrixR,
            int [] startState,
            int [] endState,
            double totalTime) {
    	setRatematrix(rateMatrixR);
    	return generatePath(startState, endState, totalTime);
    }

    public void setRatematrix(double [][] rateMatrixR);

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
