package emat.likelihood;

import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

@Description("Rooted time tree that allows editing operations for efficient restore")
public class EditableTree extends Tree {
	
	/** for tracking editing operations **/
	protected List<Edit> editList;

	public EditableTree() {
		nodeTypeInput.setValue(EditableNode.class.getName(), this);
	}
	
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		for (BEASTInterface o : getOutputs()) {
			if (o instanceof MutationState state) {
				for (BEASTInterface o2 : state.getOutputs()) {
					if (o2 instanceof EditList list) {
						editList = list.list;
					}
				}
			}
		}
		if (editList == null) {
			throw new RuntimeException("Could not find EditList in output");
		}
	}
	
	/** set height of a single node in the tree **/
	public void setHeight(int nodeNr, double height) {
		editList.add(new Edit(EditType.nodeHeightMove, nodeNr, m_nodes[nodeNr].getHeight(), height));
		m_nodes[nodeNr].setHeight(height);
	}

	public void undoHeight(int nodeNr, double height) {
		((EditableNode)m_nodes[nodeNr]).resetHeight(height);
	}

	/** apply nearest neighbour interchange move on clade above node
	 * and place it at specified height
	 *  **/
	public void doNNI(int nodeNr, double newHeight) {
		Node node = m_nodes[nodeNr];
		Node parent = node.getParent();
		if (parent == null) {
			throw new IllegalArgumentException("Cannot apply NNI to root");
		}
		Node sibling = parent.getLeft() == node ? parent.getRight() : parent.getLeft();
		Node grandParent = parent.getParent();
		if (grandParent == null) {
			throw new IllegalArgumentException("Cannot apply NNI to child of root");
		}
		Node other = grandParent.getLeft() == parent ? grandParent.getRight() : grandParent.getLeft();
		
		grandParent.removeChild(other);
		parent.removeChild(sibling);
		grandParent.addChild(sibling);
		parent.addChild(other);
		
		double oldHeight = parent.getHeight();
		parent.setHeight(newHeight);
		
		if (parent.getLength() < 0 || node.getLength() < 0) {
			throw new IllegalArgumentException("NNI gives negative branch lengths");
		}

		editList.add(new Edit(EditType.nni, nodeNr, oldHeight, newHeight));
	}
	
	public void undoNNI(int nodeNr, double height) {
		Node node = m_nodes[nodeNr];
		Node parent = node.getParent();
		if (parent == null) {
			throw new IllegalArgumentException("Cannot apply NNI to root");
		}
		Node sibling = parent.getLeft() == node ? parent.getRight() : parent.getLeft();
		Node grandParent = parent.getParent();
		if (grandParent == null) {
			throw new IllegalArgumentException("Cannot apply NNI to child of root");
		}
		Node other = grandParent.getLeft() == parent ? grandParent.getRight() : grandParent.getLeft();
		
		grandParent.getChildrenMutable().remove(other);
		parent.getChildrenMutable().remove(sibling);
		grandParent.getChildrenMutable().add(sibling);
		((EditableNode)sibling).resetParent(grandParent);
		parent.getChildrenMutable().add(other);
		((EditableNode)other).resetParent(parent);
		
		((EditableNode)parent).resetHeight(height);
		
		if (parent.getLength() < 0 || node.getLength() < 0) {
			throw new IllegalArgumentException("NNI gives negative branch lengths");
		}
	}
	
	
	public void doSPR(int subtreeNodeNr, int targetNodNr, double newHeight) {
		Node node = m_nodes[subtreeNodeNr];
		Node parent = node.getParent();
		if (parent == null) {
			throw new IllegalArgumentException("Cannot apply SPR to root");
		}
		Node sibling = parent.getLeft() == node ? parent.getRight() : parent.getLeft();
		Node grandParent = parent.getParent();
		if (grandParent == null) {
			throw new IllegalArgumentException("Cannot apply SPR to child of root");
		}
		
		Node targetNode = m_nodes[targetNodNr];
		Node targetParent = targetNode.getParent();
		if (targetParent == null) {
			throw new IllegalArgumentException("Cannot apply SPR to root");
		}
		
		grandParent.removeChild(parent);
		grandParent.addChild(sibling);
		targetParent.removeChild(targetNode);
		targetParent.addChild(parent);
		parent.removeChild(sibling);
		parent.addChild(targetNode);

		
		double oldHeight = parent.getHeight();
		parent.setHeight(newHeight);
		
		if (parent.getLength() < 0 || node.getLength() < 0) {
			throw new IllegalArgumentException("SPR gives negative branch lengths");
		}

		editList.add(new Edit(subtreeNodeNr, sibling.getNr(), parent.getNr(), targetNodNr, oldHeight));
	}
	
	public void undoSPR(int subtreeNodeNr, int targetNodNr, double oldHeight) {
		Node node = m_nodes[subtreeNodeNr];
		Node parent = node.getParent();
		if (parent == null) {
			throw new IllegalArgumentException("Cannot apply SPR to root");
		}
		Node sibling = parent.getLeft() == node ? parent.getRight() : parent.getLeft();
		Node grandParent = parent.getParent();
		if (grandParent == null) {
			throw new IllegalArgumentException("Cannot apply SPR to child of root");
		}
		
		Node targetNode = m_nodes[targetNodNr];
		Node targetParent = targetNode.getParent();
		if (targetParent == null) {
			throw new IllegalArgumentException("Cannot apply SPR to root");
		}
		
		grandParent.getChildrenMutable().remove(parent);
		parent.getChildrenMutable().remove(sibling);
		targetParent.getChildrenMutable().remove(targetNode);
		
		grandParent.getChildrenMutable().add(sibling);
		((EditableNode)sibling).resetParent(grandParent);

		targetParent.getChildrenMutable().add(parent);
		((EditableNode)parent).resetParent(targetParent);
		
		parent.getChildrenMutable().add(targetNode);
		((EditableNode)targetNode).resetParent(parent);
		
		parent.setHeight(oldHeight);
		
		if (parent.getLength() < 0 || node.getLength() < 0) {
			throw new IllegalArgumentException("SPR gives negative branch lengths");
		}
	}
	
	@Override
	protected void accept() {
		super.accept();
		hasStartedEditing = false;
	}
	
	@Override
	protected void store() {
	}
	
	@Override
	public void restore() {
		for (Edit e : editList) {
			e.undo(this);
		}
        hasStartedEditing = false;
	}
	
}
