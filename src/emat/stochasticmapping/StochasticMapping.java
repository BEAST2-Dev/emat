package emat.stochasticmapping;

import java.util.List;

public interface StochasticMapping {

    default public List<TimeStateInterval> generatePath(
            double[][] rateMatrixR,
            int startState,
            int endState,
            double totalTime) {
    	setRatematrix(rateMatrixR);
    	return generatePath(startState, endState, totalTime);
    }

    public void setRatematrix(double [][] rateMatrixR);

    public List<TimeStateInterval> generatePath(
            int startState,
            int endState,
            double totalTime);
    
}
