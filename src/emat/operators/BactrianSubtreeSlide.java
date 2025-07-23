package emat.operators;

/*
* File SubtreeSlide.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
/*
 * SubtreeSlideOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */


import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import emat.likelihood.Edit;
import emat.likelihood.EditableNode;
import emat.likelihood.EditableTree;
import emat.likelihood.MutationOnBranch;


/**
 * Implements the subtree slide move.
 */
@Description("As sub-tree slide but with Bactrian kernel to determine size of step. " +
		"Moves the height of an internal node along the branch. " +
        "If it moves up, it can exceed the root and become a new root. " +
        "If it moves down, it may need to make a choice which branch to " +
        "slide down into.")
public class BactrianSubtreeSlide extends SPR {

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by tree-height/log2(n-taxa).", -1.0);

    private double size;
    private double limit;
    private KernelDistribution kernelDistribution;

    @Override
    public void initAndValidate() {
        size = sizeInput.get();
        limit = limitInput.get();
        kernelDistribution = kernelDistributionInput.get();
        
        super.initAndValidate();
    }

    /**
     * Do a probabilistic subtree slide move.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	
		// TODO: only do this when substModel changes
		substModel.setupRateMatrix();
		setRatematrix(substModel.getRateMatrix());

		
        final EditableTree tree = (EditableTree) state.treeInput.get(); 

        double logHR = 0;

        Node i;
        // 1. choose a random node avoiding root
        i = tree.getNode(Randomizer.nextInt(tree.getNodeCount()-1));

        final Node parent = i.getParent();
        final Node sibling = getOtherChild(parent, i);
        final Node grandParent = parent.getParent();

        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldHeight = parent.getHeight();
        final double newHeight = oldHeight + delta;

        // 3. if the move is up
        if (delta > 0) {

            // 3.1 if the topology will change
            if (grandParent != null && grandParent.getHeight() < newHeight) {
                // find new parent
                Node newParent = grandParent;
                Node newChild = parent;
                while (newParent.getHeight() < newHeight) {
                    newChild = newParent;
                    newParent = newParent.getParent();
                    if (newParent == null) break;
                }
                // the moved node 'p' would become a child of 'newParent'
                //

                // 3.1.1 if creating a new root
                if (newChild.isRoot()) {
                	// TODO: implement;
                	return Double.NEGATIVE_INFINITY;
//                    replace(p, CiP, newChild);
//                    replace(PiP, p, CiP);
//
//                    p.setParent(null);
//                    tree.setRoot(p);
                }
                // 3.1.2 no new root
                else {
                	logHR += subtreePruneRegraft((EditableNode)i, (EditableNode)newChild, newHeight, oldHeight);
                }

                // p.setHeight(newHeight);
            	tree.setHeight(parent.getNr(), newHeight);

                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(newChild, oldHeight, null);
                //System.out.println("possible sources = " + possibleSources);

                logHR += -Math.log(possibleSources);

            } else {
                // just change the node height
            	logHR += slide((EditableNode) parent, newHeight);
            }
        }
        // 4 if we are sliding the subtree down.
        else {

            // 4.0 is it a valid move?
            if (i.getHeight() > newHeight) {
                return Double.NEGATIVE_INFINITY;
            }

            // 4.1 will the move change the topology
            if (sibling.getHeight() > newHeight) {

                final List<Node> newChildren = new ArrayList<>();
                final int possibleDestinations = intersectingEdges(sibling, newHeight, newChildren);

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                final int childIndex = Randomizer.nextInt(newChildren.size());
                final Node newChild = newChildren.get(childIndex);
                final Node newParent = newChild.getParent();

                // 4.1.1 if p was root
                if (parent.isRoot()) {
                    // new root is CiP
                	// TODO: implement;
                	return Double.NEGATIVE_INFINITY;
//                    replace(p, CiP, newChild);
//                    replace(newParent, newChild, p);
//
//                    CiP.setParent(null);
//                    tree.setRoot(CiP);

                } else {
                	logHR += subtreePruneRegraft((EditableNode)i, (EditableNode)newChild, newHeight, oldHeight);
                }

            	tree.setHeight(parent.getNr(), newHeight);

                logHR += Math.log(possibleDestinations);
            } else {
                logHR += slide((EditableNode)parent, newHeight);
            }
        }
        return logHR;
    }

    /** 
     * slide parent of subtree to newHeight -- without chaning topology, 
     * but with resampling mutations on branch above subtree (only if necessary) 
     * **/
	protected double slide(EditableNode subtree, double newHeight) {
		
		Node parent = subtree.getParent();
		Node sibling = getOtherChild(parent, subtree);
		int siblingNr = sibling.getNr();
		int parentNr = parent.getNr();
		
		double oldHeight = parent.getHeight();
		double newLength = newHeight - subtree.getHeight();
		double oldLength = oldHeight - subtree.getHeight();
		
		double logHR = 0;

		// split sibling and parent mutations so that they remain in the same place
		List<MutationOnBranch> siblingMutations = state.getMutationList(siblingNr);
		double siblingLength = sibling.getLength();
		List<MutationOnBranch> parentMutations = state.getMutationList(parent.getNr());
		double parentLength = parent.getLength();
		
		int [] states0 = state.getNodeSequence(siblingNr);
		int [] states = state.getNodeSequenceForUpdate(parentNr);
		int [] parentNodeStates = states.clone();
		System.arraycopy(states0, 0, states, 0, states.length);

		List<MutationOnBranch> newSiblingMutations = new ArrayList<>();
		List<MutationOnBranch> newParentMutations = new ArrayList<>();
		if (oldHeight > newHeight) {
			double f = (oldHeight - sibling.getHeight())/(siblingLength + parentLength);
			double f2 = (parent.getHeight() - sibling.getHeight()) / (parent.getHeight() - newHeight); 
			// all sibling mutations remain on sibling mutation branch
			for (MutationOnBranch m : siblingMutations) {
				newSiblingMutations.add(new MutationOnBranch(siblingNr, f * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
			}
			// parent mutations may remain be moved to sibling mutation branch
			double threshold = 1 / f;
			for (MutationOnBranch m : parentMutations) {
				if (m.getBrancheFraction() < threshold) {
					newSiblingMutations.add(new MutationOnBranch(siblingNr, f * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
					states[m.siteNr()] = m.getToState();
				} else {
					newParentMutations.add(new MutationOnBranch(parentNr, f2 * (m.brancheFraction() - threshold), m.getFromState(), m.getToState(), m.siteNr()));
				}
			}

		} else {
			 // oldHeight < newHeight
			double f = (oldHeight - sibling.getHeight())/(siblingLength + parentLength);
			double f2 = (parent.getHeight() - sibling.getHeight()) / (parent.getHeight() - newHeight); 
			double threshold = 1 / f;
			// sibling mutations may remain be moved to parent mutation branch
			for (MutationOnBranch m : siblingMutations) {
				if (m.getBrancheFraction() < threshold) {
					newSiblingMutations.add(new MutationOnBranch(siblingNr, f * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
					states[m.siteNr()] = m.getToState();
				} else {
					newParentMutations.add(new MutationOnBranch(parentNr, f2 * (m.brancheFraction() - threshold), m.getFromState(), m.getToState(), m.siteNr()));
				}
			}

			// all parent  mutations remain on parent mutation branch
			for (MutationOnBranch m : parentMutations) {
				newParentMutations.add(new MutationOnBranch(siblingNr, f * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
			}
			
		}
		double f = siblingLength/(siblingLength + parentLength);
		for (MutationOnBranch m : siblingMutations) {
			newSiblingMutations.add(new MutationOnBranch(siblingNr, f * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
		}
		for (MutationOnBranch m : parentMutations) {
			newSiblingMutations.add(new MutationOnBranch(siblingNr, f + (1-f) * m.brancheFraction(), m.getFromState(), m.getToState(), m.siteNr()));
		}


		// resample mutations on branch above subtree (only if necessary)
		int nodeNr = subtree.getNr();
		int [] nodeStates = state.getNodeSequence(nodeNr);

		double totalTime = (newHeight - subtree.getHeight()) * clockModel.getRateForBranch(subtree);
		double [] weightsN = setUpWeights(totalTime);
		double [] p = new double[M_MAX_JUMPS];

		boolean [] needsResampling = new boolean[states.length];
		for (int i = 0; i < states.length; i++) {
			if (parentNodeStates[i] != states[i]) {
				needsResampling[i] = true;
			}
		}
		// keep mutations on sites that do not differ
		List<MutationOnBranch> currentNodeMutations = state.getMutationList(nodeNr);
		List<MutationOnBranch> newNodeMutations = new ArrayList<>();
		double volumeChange = newLength / oldLength;
		for (MutationOnBranch m : currentNodeMutations) {
			if (!needsResampling[m.siteNr()]) {
				newNodeMutations.add(m);
				logHR += Math.log(volumeChange);
					
			}
		}
		// resample sites that do differ
		for (int i = 0; i < states.length; i++) {
			if (needsResampling[i]) {
				int nodeState = nodeStates[i]; 
				for (int r = 0; r < M_MAX_JUMPS; r++) {
					p[r] = weightsN[r] * qUnifPowers.get(r)[states[i]][nodeState];
				}			
				int N = FastRandomiser.randomChoicePDF(p);
				generatePath(nodeNr, i, states[i], nodeState, N, newNodeMutations);
			}
		}

		
		state.setBranchMutations(siblingNr, newSiblingMutations);
		state.setBranchMutations(nodeNr, newNodeMutations);
		state.setBranchMutations(parentNr, newParentMutations);

		EditableTree tree = (EditableTree) state.treeInput.get();
    	tree.setHeight(subtree.getNr(), newHeight);

        return logHR;
	}

    
    

	
	

	protected double getDelta() {
    	return kernelDistribution.getRandomDelta(0, Double.NaN, size);
    }

    protected int intersectingEdges(Node node, double height, List<Node> directChildren) {
        final Node parent = node.getParent();

        if (parent.getHeight() < height) return 0;

        if (node.getHeight() < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }

        if (node.isLeaf()) {
            return 0;
        } else {
            final int count = intersectingEdges(node.getLeft(), height, directChildren) +
                    intersectingEdges(node.getRight(), height, directChildren);
            return count;
        }
    }

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(size);
            final double f = Math.exp(delta);
            if( limit > 0 ) {
                final TreeInterface tree = state.treeInput.get();
                final double h = tree.getRoot().getHeight();
                final double k = Math.log(tree.getLeafNodeCount()) / Math.log(2);
                final double lim = (h / k) * limit;
                if( f <= lim ) {
                    size = f;
                }
            } else {
               size = f;
            }
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return size;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        size = value;
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newDelta = size * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newDelta);
        } else return "";
    }

}
