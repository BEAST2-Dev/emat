package emat.operators;



import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import emat.likelihood.EditableNode;
import emat.likelihood.EditableTree;
import emat.substitutionmodel.EmatSubstitutionModel;


/**
 * Implements the Nearest Neighbor Interchange (NNI) operation. This particular
 * implementation assumes explicitly bifurcating trees. It works similar to the
 * Narrow Exchange but with manipulating the height of a node if necessary.
 * Adapted from BEAST1 dr.evomodel.operators.NNI by Sebastian Hoehna 
 */
@Description("Nearest Neighbor Interchange (NNI) operation")
public class NNIOperator extends SPR {


    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }
    
    @Override
    public double proposal() {
        final int nodeCount = tree.getNodeCount();
        final Node root = tree.getRoot();

        Node node;
        // 0. determine set of candidate nodes
        do {
          node = tree.getNode(Randomizer.nextInt(nodeCount));
        } while( root == node || node.getParent() == root );

        
        // get parent node
        final Node parent = node.getParent();
        // get parent of parent -> grant parent :)
        final Node grandParent = parent.getParent();
        // get other child of grant parent -> uncle
        Node uncle = getOtherChild(grandParent, parent);

        // change the height of my father to be randomly between my uncle's
        // heights and my grandfather's height
        // this is necessary for the hastings ratio to do also if the uncle is
        // younger anyway

        final double heightGrandfather = grandParent.getHeight();
        final double heightUncle = uncle.getHeight();
        final double heightI = node.getHeight();
        final double minHeightFather  = Math.max(heightI, heightUncle);
        final Node sibling = getOtherChild(parent, node);
        final double minHeightReverse = Math.max(heightI, sibling.getHeight());

        double ran;
        do {
            ran = FastRandomiser.nextDouble();
        } while( ran == 0.0 || ran == 1.0 );

        // now calculate the new height for father between the height of the
        // uncle and the grandparent
        final double newHeightFather = minHeightFather + (ran * (heightGrandfather - minHeightFather));

        // perform SPR move so parent of node becomes parent of uncle
        double logHR = subtreePruneRegraft((EditableNode) node, (EditableNode) uncle, newHeightFather, node.getParent().getHeight(), EmatSubstitutionModel.M_MAX_JUMPS);

        // hastings ratio = backward Prob / forward Prob
        logHR += Math.log((heightGrandfather - minHeightFather) / (heightGrandfather - minHeightReverse));

        return logHR;
    }

}
