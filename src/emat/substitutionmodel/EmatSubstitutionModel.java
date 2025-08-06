package emat.substitutionmodel;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.CalculationNode;

@Description("Wrapper for GeneralSubsitutionModel -- keeps track of data structures for stochastic mapping")
public class EmatSubstitutionModel extends CalculationNode {
	final public Input<GeneralSubstitutionModel> substModelInput =
            new Input<>("substModel", "substitution model along branches in the beast.tree", null, Validate.REQUIRED);

	
	private GeneralSubstitutionModel substModel;
	private boolean needsUpdate = true;

	protected int stateCount;

	public static int M_MAX_JUMPS = 20;
	
	protected double lambdaMax;
	protected double[][] qUnif;
	protected List<double[][]> qUnifPowers;

	public EmatSubstitutionModel() {		
	}
	
	public EmatSubstitutionModel(SubstitutionModel substitutionModel) {
		substModelInput.setValue(substitutionModel, this);
		initAndValidate();
	}

	@Override
	public void initAndValidate() {
		substModel = substModelInput.get();
		substModel.setupRelativeRates();
		substModel.setupRateMatrix();
		
		stateCount = substModel.getStateCount();
		
		setRatematrix(substModel.getRateMatrix());
	}

	public double[] getFrequencies() {
		return substModel.getFrequencies();
	}

	public double[][] getRateMatrix() {
		return substModel.getRateMatrix();
	}
	
	
	public void setRatematrix(double[][] rateMatrixR) {
		int numStates = rateMatrixR.length;

		// --- Step 0: Precomputation & Initialization ---
		lambdaMax = 0.0;
		for (int i = 0; i < numStates; i++) {
			if (-rateMatrixR[i][i] > lambdaMax) {
				lambdaMax = -rateMatrixR[i][i];
			}
		}

		qUnif = getQUnif(rateMatrixR, lambdaMax);

		// Precompute powers of Q_unif to avoid re-computation
		qUnifPowers = new ArrayList<>();
		qUnifPowers.add(identity(numStates)); // Q_unif^0

		for (int n = 0; n <= M_MAX_JUMPS; n++) {
			if (n > 0) {
				qUnifPowers.add(multiply(qUnifPowers.get(n - 1), qUnif));
			}
		}
		
		needsUpdate = false;
	}

	
    /** return matrix I+R/lambdaMax **/
    private double [][] getQUnif(double [][] rateMatrixR, double lambdaMax) {
//		return add(
//		        identity(numStates),
//		        multiplyByScalar(rateMatrixR, 1.0 / lambdaMax)
        int rows = rateMatrixR.length;
        int cols = rateMatrixR[0].length;
        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = rateMatrixR[i][j] /lambdaMax;
            }
            result[i][i] += 1.0;
        }
        return result;
    }
    
    
	public static double[][] multiply(double[][] a, double[][] b) {
        int rowsA = a.length;
        int colsA = a[0].length;
        int colsB = b[0].length;
        if (colsA != b.length) {
            throw new IllegalArgumentException("Matrix dimensions incompatible for multiplication.");
        }
        double[][] result = new double[rowsA][colsB];
        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < colsB; j++) {
                for (int k = 0; k < colsA; k++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return result;
    }

	public static double[][] identity(int size) {
        double[][] id = new double[size][size];
        for (int i = 0; i < size; i++) {
            id[i][i] = 1.0;
        }
        return id;
    }

	public double getLambdaMax() {
		if (needsUpdate) {
			setRatematrix(substModel.getRateMatrix());
		}
		return lambdaMax;
	}

	public List<double[][]> getQUnifPowers() {
		if (needsUpdate) {
			setRatematrix(substModel.getRateMatrix());
		}
		return qUnifPowers;
	}

	public int getStateCount() {
		return stateCount;
	}

	public double[][] getQUnifPowers(int remainingJumps) {
		if (needsUpdate) {
			setRatematrix(substModel.getRateMatrix());
		}
		return qUnifPowers.get(remainingJumps);
	}

	public double[][] getQUnif() {
		if (needsUpdate) {
			setRatematrix(substModel.getRateMatrix());
		}
		return qUnif;
	}	
	
	@Override
	protected void store() {
		super.store();
	}
	
	@Override
	protected void restore() {
		needsUpdate = true;
		super.restore();
	}
	
	@Override
	protected boolean requiresRecalculation() {
		needsUpdate = true;
		return true;
	}
}
