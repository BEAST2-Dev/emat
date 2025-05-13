package emat.likelihood;

import emat.likelihood.MutationState.EditType;

/** 
 * class for tracking edits to the MutationState 
 * used to restore a MutationState if necessary
 **/
class Edit {
	EditType type;
	int siteNr;
	int nodeNr;
	MutationOnBranch oldValue;
	MutationOnBranch newValue;
	

	
	public Edit(EditType type, int siteNr, int nodeNr, MutationOnBranch mutation) {
		this.type = type;
		this.siteNr = siteNr;
		this.nodeNr = nodeNr;
		this.oldValue = mutation;
		this.newValue = mutation;
	}



	public Edit(EditType type, int siteNr, int nodeNr, MutationOnBranch oldMutation,
			MutationOnBranch newMutation) {
		this.type = type;
		this.siteNr = siteNr;
		this.nodeNr = nodeNr;
		this.oldValue = oldMutation;
		this.newValue = newMutation;
	}



	void undo(MutationState state) {
		switch(type) {
		case add:
			state.deleteMutation0(siteNr, nodeNr, newValue);
			break;
		case delete:
			state.addMutation0(siteNr, nodeNr, oldValue.brancheFraction, oldValue.stateTransition);
			break;
		case replace:
			state.replaceMutation0(siteNr, nodeNr, oldValue, newValue);
			break;
		}
	}
}