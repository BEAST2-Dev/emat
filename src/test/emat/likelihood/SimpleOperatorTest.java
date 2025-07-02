package test.emat.likelihood;


import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GTR;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import emat.likelihood.EditList;
import emat.likelihood.EditableTree;
import emat.likelihood.MutationState;
import emat.likelihood.MutationStateTreeLikelihood;
import emat.likelihood.ParsimonyMutationStateInitialiser;
import emat.operators.BactrianNodeOperator;
import emat.operators.MutationMover;
import emat.operators.MutationOnNodeResampler;
import emat.operators.SPR;

public class SimpleOperatorTest {

	@Test
	public void testOperators() {
        Sequence a = new Sequence("A", "A A C G T TT");
        Sequence b = new Sequence("B", "A C C G T CC");
        Sequence c = new Sequence("C", "A A C G T TT");

        Alignment data = new Alignment();
        data.initByName("sequence", a, "sequence", b, "sequence", c, "dataType", "nucleotide");

		TreeParser newick = new TreeParser();
		newick.initByName("taxa", data,
                "newick", "((A:1,B:1):1,C:2)",
                "IsLabelledNewick", true);

		EditableTree tree = new EditableTree();
        tree.assignFrom(newick);
        
        MutationState mutationState = new MutationState();
        EditList editList = new EditList();
        editList.mutationStateInput.setValue(mutationState, editList);
        mutationState.initByName("tree", tree, "data", data);

        // sets EditList in tree
        tree.initAndValidate();
        
        ParsimonyMutationStateInitialiser init = new ParsimonyMutationStateInitialiser();
        init.initByName("mutationState", mutationState, "tree", tree);
        init.initStateNodes();
        
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "0.3 0.2 0.2 0.3");

        GTR gtr = new GTR();
        gtr.initByName("frequencies", freqs);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", gtr);
        
        StrictClockModel clockModel = new StrictClockModel();
        RealParameter clockRate = new RealParameter("1.0");
        clockModel.initByName("clock.rate", clockRate);

        MutationStateTreeLikelihood likelihood = new MutationStateTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "mutationState", mutationState, 
        		"editList", editList,
        		"branchRateModel", clockModel);
        
        double logP = likelihood.calculateLogP();
        
        State state = new State();
        state.initByName("stateNode", mutationState, "stateNode", tree);
        state.initialise();
        state.setPosterior(likelihood);
        
        state.store(0);

        Randomizer.setSeed(127);
//        {   
//	        MutationMover operator = new MutationMover();
//	        operator.initByName("mutationState", mutationState, "weight", 1.0);
//	        operator.proposal();
//            state.storeCalculationNodes();
//	        state.checkCalculationNodesDirtiness();
//	        double logP2 = likelihood.calculateLogP();
//	        assertNotEquals(logP, logP2);
//	        
//	        state.restore();
//	        state.restoreCalculationNodes();
//	        double logP3 = likelihood.calculateLogP();
//	        
//	        assertEquals(logP, logP3, 1e-10);
//        }
        
        {
	        BactrianNodeOperator operator = new BactrianNodeOperator();
	        operator.initByName("tree", tree, "weight", 1.0);
	        operator.proposal();
            state.storeCalculationNodes();
	        state.checkCalculationNodesDirtiness();
	        double logP2 = likelihood.calculateLogP();
	        assertNotEquals(logP, logP2);
	        
	        state.restore();
	        state.restoreCalculationNodes();
	        double logP3 = likelihood.calculateLogP();
	        
	        assertEquals(logP, logP3, 1e-10);
        }

//        {
//	        NNIOperator operator = new NNIOperator();
//	        operator.initByName("tree", tree, "weight", 1.0);
//	        operator.proposal();
//            state.storeCalculationNodes();
//	        state.checkCalculationNodesDirtiness();
//	        double logP2 = likelihood.calculateLogP();
//	        assertNotEquals(logP, logP2);
//	        
//	        state.restore();
//	        state.restoreCalculationNodes();
//	        double logP3 = likelihood.calculateLogP();
//	        
//	        assertEquals(logP, logP3, 1e-10);
//        }

//        {
//        	double r = -0;
//        	
//            Randomizer.setSeed(127);
//	        SPR operator = new SPR();
//	        operator.initByName("tree", tree, "weight", 1.0, "mutationState", mutationState, "likelihood", likelihood);
//	        operator.proposal();
//            state.storeCalculationNodes();
//	        state.checkCalculationNodesDirtiness();
//	        double logP2 = likelihood.calculateLogP();
//	        assertNotEquals(logP, logP2);
//	        
//	        state.restore();
//	        state.restoreCalculationNodes();
//	        double logP3 = likelihood.calculateLogP();
//	        
//	        assertEquals(logP, logP3, 1e-10);
//        }
//
//        {
//            Randomizer.setSeed(127);
//            MutationOnBranchResampler operator = new MutationOnBranchResampler();
//	        operator.initByName("weight", 1.0, "mutationState", mutationState, "likelihood", likelihood);
//	        operator.proposal();
//            state.storeCalculationNodes();
//	        state.checkCalculationNodesDirtiness();
//	        double logP2 = likelihood.calculateLogP();
//	        assertNotEquals(logP, logP2);
//	        
//	        state.restore();
//	        state.restoreCalculationNodes();
//	        double logP3 = likelihood.calculateLogP();
//	        
//	        assertEquals(logP, logP3, 1e-10);
//        }

        {
            Randomizer.setSeed(127);
            MutationOnNodeResampler operator = new MutationOnNodeResampler();
	        operator.initByName("weight", 1.0, "mutationState", mutationState, "likelihood", likelihood);
	        operator.proposal();
            state.storeCalculationNodes();
	        state.checkCalculationNodesDirtiness();
	        double logP2 = likelihood.calculateLogP();
	        assertNotEquals(logP, logP2);
	        
	        state.restore();
	        state.restoreCalculationNodes();
	        double logP3 = likelihood.calculateLogP();
	        
	        assertEquals(logP, logP3, 1e-10);
        }
	}
	
}
