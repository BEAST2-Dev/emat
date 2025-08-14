package emat.operators;



import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import emat.likelihood.Edit;
import emat.likelihood.EditableNode;
import emat.likelihood.MutationOnBranch;
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
        
        double [] logHR = {0};
        boolean targeted = false;
        if (FastRandomiser.nextDouble() < targetedInput.get()) {
        	node = MutationOperatorUtil.selectNodeByMutationCount(logHR, tree, state);
        	targeted = true;
        } else {
	        do {
	          node = tree.getNode(FastRandomiser.nextInt(nodeCount));
	        } while( root == node || node.getParent() == root );
        }

        
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
        logHR[0] += FastRandomiser.nextDouble() >  resampleProbabilityInput.get()
        		? subtreePruneRegraft((EditableNode) node, (EditableNode) uncle, newHeightFather, node.getParent().getHeight(), EmatSubstitutionModel.M_MAX_JUMPS)
        		: NNIAndResample((EditableNode) node, (EditableNode) uncle, newHeightFather, node.getParent().getHeight(), EmatSubstitutionModel.M_MAX_JUMPS);

        // hastings ratio = backward Prob / forward Prob
        logHR[0] += Math.log((heightGrandfather - minHeightFather) / (heightGrandfather - minHeightReverse));

        if (targeted) {
        	logHR[0] += MutationOperatorUtil.logHRUpdate(node, tree, state);
        }
        return logHR[0];
    }
    
    
    



	private static boolean debug = true;
    
    /**
	 * Perform an NNI move & resample grandparent sequence
	 * @param subtree = MRCA of the subtree to be removed
	 * @param targetBranch = MRCA above which the subtree will be grafted 
	 * @param newHeight = height at which the subtree will be grafted
	 */
	protected double NNIAndResample(EditableNode subtree, EditableNode targetBranch, double newHeight,  double oldHeight, final int M_MAX_JUMPS) {
		
		Node parent = subtree.getParent();
		int parentNr = parent.getNr();
		
		double logHR = 0;

		
		if (false || parent.getParent().isRoot()) {
			if (true) {
				return Double.NEGATIVE_INFINITY;
			}
			// NNI below root allows resampling all mutations
			// on branches below root as well as below parent
			Node root = parent.getParent();

	        logHR += -state.getMutationList(root.getLeft().getNr()).size() * Math.log(root.getLeft().getLength());
	        logHR += -state.getMutationList(root.getRight().getNr()).size() * Math.log(root.getRight().getLength());
	        logHR += -state.getMutationList(parent.getLeft().getNr()).size() * Math.log(parent.getLeft().getLength());
	        logHR += -state.getMutationList(parent.getRight().getNr()).size() * Math.log(parent.getRight().getLength());
			
			Edit e = tree.doSPR(subtree.getNr(), targetBranch.getNr(), newHeight);
			

			int [] rootSequence = state.getNodeSequenceForUpdate(root.getNr());
			int [] nodeSequence = state.getNodeSequenceForUpdate(parent.getNr());
			List<MutationOnBranch> branchMutationsLeftOfRoot = new ArrayList<>();
			List<MutationOnBranch> branchMutationsRightOfRoot = new ArrayList<>();
			List<MutationOnBranch> branchMutationsLeft = new ArrayList<>();
			List<MutationOnBranch> branchMutationsRight = new ArrayList<>();
			
			MutationOperatorUtil.resampleBelowRoot(parent, M_MAX_JUMPS,
					branchMutationsLeftOfRoot,
					branchMutationsRightOfRoot,
					branchMutationsLeft,
					branchMutationsRight,
					rootSequence,
					nodeSequence,
					state,
					substModel,
					clockModel);
			
			state.setBranchMutations(root.getLeft().getNr(), branchMutationsLeftOfRoot);
			state.setBranchMutations(root.getRight().getNr(), branchMutationsRightOfRoot);		
			state.setBranchMutations(parent.getLeft().getNr(), branchMutationsLeft);
			state.setBranchMutations(parent.getRight().getNr(), branchMutationsRight);

	        logHR += state.getMutationList(root.getLeft().getNr()).size() * Math.log(root.getLeft().getLength());
	        logHR += state.getMutationList(root.getRight().getNr()).size() * Math.log(root.getRight().getLength());
	        logHR += state.getMutationList(parent.getLeft().getNr()).size() * Math.log(parent.getLeft().getLength());
	        logHR += state.getMutationList(parent.getRight().getNr()).size() * Math.log(parent.getRight().getLength());
			return logHR;

		}
		
        logHR += -state.getMutationList(subtree.getNr()).size() * Math.log(subtree.getLength());
		
		// subtree's grand parent is not the root
		
		// split target branch mutations
		int targetNr = targetBranch.getNr();
		List<MutationOnBranch> targetMutations = state.getMutationList(targetNr);
		double targetLength = targetBranch.getLength();
		List<MutationOnBranch> newTargetMutations = new ArrayList<>();
		double f = targetLength / (newHeight - targetBranch.getHeight());
		double threshold = 1 / f;
		int [] states0 = state.getNodeSequence(targetNr);
		int [] states = state.getNodeSequenceForUpdate(parentNr);
		int [] parentNodeStates = states.clone();
		System.arraycopy(states0, 0, states, 0, states.length);
		int [] statesOrig = null;
		if (debug) {
			statesOrig = new int[states.length];
			System.arraycopy(states0, 0, statesOrig, 0, states.length);
		}
		for (MutationOnBranch m : targetMutations) {
			if (m.getBrancheFraction() < threshold) {
				states[m.siteNr()] = m.getToState();
				newTargetMutations.add(new MutationOnBranch(targetNr, f * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
			}
			if (debug) {
				statesOrig[m.siteNr()] = m.getToState();
			}
		}

		if (debug) {
			int [] gpstates = state.getNodeSequence(targetBranch.getParent().getNr());
			
			for (int i = 0; i < states.length; i++) {
				if (statesOrig[i] != gpstates[i]) {
					System.err.println("Something wrong with the state reconstruction at site "  + i);
				}
			}
		}		
		
		// resample mutations on branch above subtree (only if necessary)
		int nodeNr = subtree.getNr();
		int [] nodeStates = state.getNodeSequence(nodeNr);

		List<MutationOnBranch> nodeMutations = new ArrayList<>();
		double totalTime = (newHeight - subtree.getHeight()) * clockModel.getRateForBranch(subtree);
		double [] weightsN = MutationOperatorUtil.setUpWeights(totalTime, substModel.getLambdaMax(), M_MAX_JUMPS);
		double [] p = new double[M_MAX_JUMPS];

		boolean [] needsResampling = new boolean[states.length];
		for (int i = 0; i < states.length; i++) {
			if (parentNodeStates[i] != states[i]) {
				needsResampling[i] = true;
			}
		}
		// keep mutations on sites that do not differ
		List<MutationOnBranch> currentNodeMutations = state.getMutationList(nodeNr);
		//double volumeChange = (newHeight - subtree.getHeight())/(oldHeight - subtree.getHeight());
		for (MutationOnBranch m : currentNodeMutations) {
			if (!needsResampling[m.siteNr()]) {
				nodeMutations.add(m);
				//logHR += Math.log(volumeChange);
					
			}
		}
		// resample sites that do differ
		List<double[][]> qUnifPowers = substModel.getQUnifPowers();
		for (int i = 0; i < states.length; i++) {
			if (needsResampling[i]) {
				int nodeState = nodeStates[i]; 
				for (int r = 0; r < M_MAX_JUMPS; r++) {
					p[r] = weightsN[r] * qUnifPowers.get(r)[states[i]][nodeState];
				}			
				int N = FastRandomiser.randomChoicePDF(p);
				MutationOperatorUtil.generatePath(nodeNr, i, states[i], nodeState, N, qUnifPowers, nodeMutations);
			}
		}

		
		Edit e = tree.doSPR(subtree.getNr(), targetBranch.getNr(), newHeight);
		
		Node _node = parent.getParent();
		int _nodeNr = _node.getNr();
		List<MutationOnBranch> branchMutationsLeft = new ArrayList<>();
		List<MutationOnBranch> branchMutationsRight = new ArrayList<>();
		
		if (parent.getParent().isRoot()) {
			MutationOperatorUtil.resampleRoot(parent.getParent(), M_MAX_JUMPS, 
					branchMutationsLeft, 
					branchMutationsRight, 
					nodeStates, 
					state, substModel, clockModel);
		} else {
			int [] nodeSequence = state.getNodeSequenceForUpdate(_nodeNr);
			List<MutationOnBranch> branchMutations = new ArrayList<>();
			MutationOperatorUtil.resample(_node, M_MAX_JUMPS, 
				branchMutations, 
				branchMutationsLeft, 
				branchMutationsRight,
				nodeSequence,
				state, substModel, clockModel);
			state.setBranchMutations(_nodeNr, branchMutations);

			state.setBranchMutations(nodeNr, nodeMutations);
			state.setBranchMutations(targetNr, newTargetMutations);
		}

		state.setBranchMutations(_node.getLeft().getNr(), branchMutationsLeft);
		state.setBranchMutations(_node.getRight().getNr(), branchMutationsRight);		
		
        logHR += state.getMutationList(subtree.getNr()).size() * Math.log(subtree.getLength());
        
        //return 0;
        return logHR;
	}

}
