package emat.util;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.datatype.DataType;
import beast.base.core.Input.Validate;
import emat.likelihood.MutationState;

@Description("Logs sequences of internal nodes for a MutationState. "
		+ "This only makes sense when the tree topology is fixed.")
public class SequenceLogger extends BEASTObject implements Loggable {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree to be initialised", Validate.REQUIRED);

    private MutationState state;
    private DataType datatype;
    private int siteCount;
    
    @Override
    public void initAndValidate() {
		state = stateInput.get();
		datatype = state.dataInput.get().getDataType();
		siteCount = state.dataInput.get().getSiteCount();
	}

	@Override
	public void init(PrintStream out) {
		int n = state.dataInput.get().getTaxonCount();
		for (int i = n; i < 2*n-1; i++) {
			for (int j = 0; j < siteCount; j++) {
				out.append("sequence" + i + "_" + j + "\t");
			}
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		int n = state.dataInput.get().getTaxonCount();
		for (int i = n; i < 2*n-1; i++) {
			int [] states = state.getNodeSequence(i);
			for (int j = 0; j < siteCount; j++) {
				out.append(datatype.getCharacter(states[j])+"\t");
			}
		}
	}

	@Override
	public void close(PrintStream out) {
	}

}
