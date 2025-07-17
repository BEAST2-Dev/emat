package emat.likelihood;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;

@Description("Keep track of list of edit actions on a Mutation State")
public class EditList extends CalculationNode {
	final public Input<MutationState> mutationStateInput = new Input<>("mutationState", "mutationState on which this operation is performed", Validate.REQUIRED);
    
	private MutationState mutationState;
	
	public List<Edit> list = new ArrayList<>();
	
	@Override
	public void initAndValidate() {
		mutationState = mutationStateInput.get();
	}
	
	@Override
	protected void accept() {
		list.clear();
		super.accept();
	}
	
	@Override
	protected void store() {
		// list.clear();
		super.store();
	}
	
	@Override
	protected void restore() {
		for (Edit e : list) {
			e.undo(mutationState);
		}
		list.clear();
		super.restore();
	}

	public void add(Edit edit) {
		list.add(edit);
	}
}
