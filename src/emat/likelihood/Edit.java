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

	public Edit(EditType type,  int nodeNr, double oldHeight, double newHeight) {
		this.type = type;
		this.nodeNr = nodeNr;
		this.oldValue = oldHeight;
		this.newValue = newHeight;
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
			break;
		case nodeHeightMove:
			break;
		case nni:
			break;
		}
	}


	void undo(EditableTree tree) {
		switch(type) {
		case addMutation:
		case deleteMutation:
		case replaceMutation:
		case moveBranchFraction:
			// no change to tree, so nothing to do 
			break;
		case nodeHeightMove:
			tree.undoHeight(nodeNr, (Double) oldValue);
			break;
		case nni:
			tree.undoNNI(nodeNr, (Double) oldValue);
			break;
		}
	}

	void apply(MutationStateTreeLikelihood state) {
		switch(type) {
		case addMutation:
			break;
		case deleteMutation:
			break;
		case replaceMutation:
			break;
		case moveBranchFraction:
			break;
		case nodeHeightMove:
			break;
		case nni:
			break;
		}
	}


}