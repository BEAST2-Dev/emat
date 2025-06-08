package emat.stochasticmapping;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.distribution.GeometricDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.util.Randomizer;

public class DelphiStochasticMapping implements StochasticMapping {
	
	

	@Override
	public void setRatematrix(double[][] rateMatrixR) {
		// TODO Auto-generated method stub

	}
	
	/** Jukes-Cantor model with mutation rate μ **/
	private double mu = 1.0;
	
	public void setMu(double mu) {
		this.mu = mu;
	}
	
	@Override
	public List<TimeStateInterval> generatePath(final int [] startState, final int [] endState, double totalTime) {
//		Algorithm 1 Efficient stochastic mapping for a branch of length τ using a Jukes-Cantor model with mutation rate μ ̃, with starting and ending states A at all L sites.
//		  p0 ← e−μ ̃τ
		double muT = mu * totalTime;
		double p0 = Math.exp(-muT);
//		p1 ←μ ̃τe−μ ̃τ
		double p1 = muT * Math.exp(-muT);
//		p∗ ←p0/(1−p1) 
		double pstar = p0 / (1 - p1);
//		l←0
		int l = 0;
		
		GeometricDistribution g = new GeometricDistribution(pstar);
		PoissonDistribution poisson = new PoissonDistribution(muT);
		// lowerPoisson = probability of sampling 0 or 1 from poisson
		double lowerPoisson = poisson.cumulativeProbability(1);
		
		int L = startState.length;
		
		List<TimeStateInterval> mutations = new ArrayList<>();
		
//		while l < L do
		while (l < L) {
//			∆ ∼ Geom(1 − p∗) 
			int delta = g.inverseCumulativeProbability(Randomizer.nextDouble());
//			l←l+∆
			l = l + delta;
//			if l<Lthen
			if (l < L) {
				List<TimeStateInterval> mutationsForSite = new ArrayList<>();
//				n ∼ KPois(μ ̃τ, 2)
//				▷ KPois(λ, k) is the restriction of Pois(λ) to at least k events.
				int n = poisson.inverseCumulativeProbability(lowerPoisson + (1-lowerPoisson) * Randomizer.nextDouble()); 
//				s0 ← A
				int state = startState[l];
				int [] states = new int[n+1];
				states[0] = state;
//				for i ← 1 to n do
				for (int i = 1; i <= n; i++) {
//					si ∼DUnif({A,C,G,T}\{si−1})
					int nextState = Randomizer.nextInt(3);
					if (nextState >= state) {
						nextState++;
					}
					states[i] = state;
					state = nextState;
//					end for
				}
//				if sn = A then
				if (state == endState[l]) {
//					t1 ≤ ... ≤ tn ∼ Unif(0,τ)
					double [] times = new double[n];
					for (int i = 1; i < n; i++) {
						times[i] = Randomizer.nextDouble() * totalTime;
					}
					times[n] = totalTime;
					Arrays.sort(times);
					for (int i = 0; i < n; i++) {
						mutationsForSite.add(new TimeStateInterval(states[i], states[i+1], times[i], times[i+1]));
					}					
//					Output mutations s0ls1 at t1, . . . , sn−1lsn at tn 
//					l←l+1
					l++;
//				else
//					(Rejection sampling starts anew at site l, not l + 1). 
//				end if
				}
//			end if 
			}
//		end while
		}
		return mutations;
	}

}
