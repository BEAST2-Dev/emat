package emat.stochasticmapping;

// Record to store path segments of a stochastic mapping
public record TimeStateInterval(int site, int state, double startTime, double endTime) {

	@Override
	public String toString() {
		return String.format("Site %d State %d from %.4f to %.4f (duration %.4f)", site, state, startTime, endTime,
				endTime - startTime);
	}

}
