package emat.util;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Loggable;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beastclassic.evolution.likelihood.AncestralStateTreeLikelihood;
import emat.likelihood.MutationOnBranch;
import emat.operators.SimpleMutationOnNodeResampler;

@Description("As MutationTreeWithMetaDataLogger, but re-simulating mutations on every log action -- "
		+ "good for comparing with Felsenstein's algorithm.")
public class MutationTreeWithMetaDataLogger2 extends AncestralStateTreeLikelihood implements Function, Loggable {
	
	
	private GeneralSubstitutionModel substModel;
	private BranchRateModel clockModel;
	private int siteCount;

	// for generating mutation mappings
    private double [][] pMatrix;
    private double [][] qMatrix;
    private double gamma;

	@Override
	public void init(PrintStream out) {
		((Tree)treeInput.get()).init(out);
		
		dataType =  dataInput.get().getDataType();
		SiteModelInterface siteModel = siteModelInput.get();
		substModel = (GeneralSubstitutionModel) ((SiteModel) siteModel).getSubstitutionModel();
		clockModel = branchRateModelInput.get();
	
		substModel.setupRateMatrix();
		qMatrix = substModel.getRateMatrix();
        this.gamma = SimpleMutationOnNodeResampler.calculateUniformisationRate(qMatrix);
        this.pMatrix = SimpleMutationOnNodeResampler.createDtmsTransitionMatrix(qMatrix, gamma);
		
		siteCount = dataInput.get().getSiteCount();
	}
	
	@Override
	public void log(long sample, PrintStream out) {
		calculateLogP();
		redrawAncestralStates();
        out.print("tree STATE_" + sample + " = ");
        TreeInterface tree = treeInput.get();
        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot()));
        out.print(";");
	}

	
    String toNewick(Node node) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft()));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight()));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }

        
        if (!node.isRoot()) {
    	    int [] patternstates = getStatesForNode(treeInput.get(), node);
    	    int [] parentPatternstates = getStatesForNode(treeInput.get(), node.getParent());
    	    double branchLength = node.getLength() * clockModel.getRateForBranch(node);
    	    
    	    List<MutationOnBranch> mutations = new ArrayList<>();
    	    for (int i = 0; i < siteCount; i++) {
    	    	int startState = patternstates[dataInput.get().getPatternIndex(i)];
    	    	int endState = parentPatternstates[dataInput.get().getPatternIndex(i)];
    	    	SimpleMutationOnNodeResampler.generatePath(node.getNr(), i, endState, startState, branchLength, mutations, pMatrix, gamma);	
    	    }
    	    
    		buf.append("[&mutationcount=" + mutations.size());

    		Set<Integer> sitesWithMutations = new HashSet<>();
    		for (MutationOnBranch m : mutations) {
    			sitesWithMutations.add(m.siteNr());
    		}
    		int multiSiteMutations = mutations.size() - sitesWithMutations.size();
    		buf.append(",multiSiteMutations=" + multiSiteMutations + "]");
    	    
        }

	    buf.append(':');
        buf.append(node.getLength());
        return buf.toString();
    }
    
	@Override
	public void close(PrintStream out) {
		((Tree)treeInput.get()).close(out);
	}	

}

