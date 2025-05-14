package emat.likelihood;

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
	float oldBranchFraction;
	

	
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



	public Edit(EditType type, int siteNr, int nodeNr, MutationOnBranch mutation,
			float oldBranchFraction, float newBranchFraction) {
		this.type = type;
		this.siteNr = siteNr;
		this.nodeNr = nodeNr;
		this.oldValue = mutation;
		this.oldBranchFraction = oldBranchFraction;
		
	}



	void undo(MutationState state) {
		switch(type) {
		case addMutation:
			state.deleteMutation0(siteNr, nodeNr, newValue);
			break;
		case deleteMutation:
			state.addMutation0(siteNr, nodeNr, oldValue.brancheFraction, oldValue.stateTransition);
			break;
		case replaceMutation:
			state.replaceMutation0(siteNr, nodeNr, oldValue, newValue);
			break;
		case moveBranchFraction:
			state.moveBranchFraction0(oldValue, siteNr, nodeNr, oldBranchFraction);
		}
	}
}