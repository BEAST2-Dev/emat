package emat.likelihood;

import java.util.ArrayList;
import java.util.List;

/** 
 * class for tracking edits to the MutationState 
 * used to restore a MutationState if necessary
 **/
class Edit {
	EditType type;
	int siteNr;
	int nodeNr;
	int parentNr, targetNr, siblingNr;
	Object oldValue;
	Object newValue;
	MutationOnBranch mutation;
	

	public Edit(EditType type) {
		this.type = type;
	}
	
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
			double oldBranchFraction, double newBranchFraction) {
		this.type = type;
		this.siteNr = siteNr;
		this.nodeNr = nodeNr;
		this.mutation = mutation;
		this.oldValue = oldBranchFraction;
		this.newValue = newBranchFraction;
		
	}

	public Edit(EditType type,  int nodeNr, double oldHeight, double newHeight) {
		this.type = type;
		this.nodeNr = nodeNr;
		this.oldValue = oldHeight;
		this.newValue = newHeight;
	}


	public Edit(EditType type, int nodeNr, List<MutationOnBranch> list) {
		this.type = type;
		this.nodeNr = nodeNr;
		List<MutationOnBranch> copy = new ArrayList<>();
		copy.addAll(list);
		this.oldValue = copy;
	}

	
	Edit(int nodeNr, int siblingNr, int parentNr, int targetNr,
			double oldHeight) {
		this.type = EditType.spr;
		this.nodeNr = nodeNr;
		this.siblingNr = siblingNr;
		this.parentNr = parentNr;
		this.targetNr = targetNr;
	
		this.oldValue = oldHeight;
	}
	
	void undo(MutationState state) {
		switch(type) {
		case moveBranchFraction:
			state.moveBranchFraction0((MutationOnBranch)mutation, siteNr, nodeNr, (Double) oldValue);
			break;
		case nodeHeightMove:
			break;
		case nni:
			// TODO
			break;
		case spr:
			// should not get here : user SPREdit instead
			break;
		case resample:
			state.restoreMutations(nodeNr, (List<MutationOnBranch>) oldValue);
			break;
		}
	}


	void undo(EditableTree tree) {
		switch(type) {
		case moveBranchFraction:
			// no change to tree, so nothing to do 
			break;
		case nodeHeightMove:
			tree.undoHeight(nodeNr, (Double) oldValue);
			break;
		case nni:
			tree.undoNNI(nodeNr, (Double) oldValue);
			break;
		case spr:
			tree.undoSPR(nodeNr, targetNr, (Double) oldValue);
			break;
		case resample:
			// nothing to do
			break;

		}
	}

	void apply(MutationStateTreeLikelihood treelikelihood) {
		switch(type) {
		case moveBranchFraction:
			treelikelihood.moveBranchFraction(nodeNr, (double) (Double) newValue);
			break;
		case nodeHeightMove:
			treelikelihood.moveNode(nodeNr, (double) (Double) newValue);
			break;
		case nni:
			// TODO
			break;
		case spr:
			// should not get here : user SPREdit instead
			break;
		case resample:
			treelikelihood.recalcBranchContribution(nodeNr);
			break;
		}
	}



	public void undo(MutationStateTreeLikelihood treelikelihood) {
		switch(type) {
		case moveBranchFraction:
			treelikelihood.undoMoveBranchFraction(nodeNr);
			break;
		case nodeHeightMove:
			treelikelihood.undoMoveNode(nodeNr);
			break;
		case nni:
			// TODO
			break;
		case spr:
			// should not get here : user SPREdit instead
			break;
		case resample:
			treelikelihood.undoBranchContribution(nodeNr);
			break;
		}
	}


}