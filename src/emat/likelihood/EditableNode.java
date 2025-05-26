package emat.likelihood;

import beast.base.evolution.tree.Node;

public class EditableNode extends Node {

	public void resetHeight(double height) {
		this.height = height;
	}
	
	public void resetParent(Node parent) {
		this.parent = parent;
	}

    public String toNewick(boolean onlyTopology) {
        final StringBuilder buf = new StringBuilder();
        if (!isLeaf()) {
            buf.append("(");
            boolean isFirst = true;
            for (Node child : getChildren()) {
                if (isFirst)
                    isFirst = false;
                else
                    buf.append(",");
                buf.append(child.toNewick(onlyTopology));
            }
            buf.append(")");

            if (getID() != null)
                buf.append(getID());
        } else {
            if (getID() != null)
                buf.append(getID());
            else
                buf.append(labelNr);
        }

        if (!onlyTopology) {
            buf.append(":").append(getLength());
        }
        return buf.toString();
    }

}
