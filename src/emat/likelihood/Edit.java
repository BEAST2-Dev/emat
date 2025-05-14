package emat.likelihood;

/** 
 * class for tracking edits to the MutationState 
 * used to restore a MutationState if necessary
 **/
class Edit {
	EditType type;
	int siteNr;
	int nodeNr;
	Object oldValue;
	Object newValue;
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
			state.deleteMutation0(siteNr, nodeNr, (MutationOnBranch)newValue);
			break;
		case deleteMutation:
			state.addMutation0(siteNr, nodeNr, ((MutationOnBranch)oldValue).brancheFraction, ((MutationOnBranch)oldValue).stateTransition);
			break;
		case replaceMutation:
			state.replaceMutation0(siteNr, nodeNr, (MutationOnBranch)oldValue, (MutationOnBranch)newValue);
			break;
		case moveBranchFraction:
			state.moveBranchFraction0((MutationOnBranch)oldValue, siteNr, nodeNr, oldBranchFraction);
		}
	}
}