package test.emat.likelihood;

import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GTR;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.State;
import beast.base.util.Randomizer;
import emat.likelihood.EditList;
import emat.likelihood.MutationState;
import emat.likelihood.MutationStateTreeLikelihood;
import emat.likelihood.ParsimonyMutationStateInitialiser;
import emat.operators.MutationMover;

public class MutationMoverTest {

	@Test
	public void testMutationMover() {
        Sequence a = new Sequence("A", "A A C G T TT");
        Sequence b = new Sequence("B", "A C C G T CC");
        Sequence c = new Sequence("C", "A A C G T TT");

        Alignment data = new Alignment();
        data.initByName("sequence", a, "sequence", b, "sequence", c, "dataType", "nucleotide");

		TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((A:1,B:1):1,C:2)",
                "IsLabelledNewick", true);
        
        MutationState mutationState = new MutationState();
        EditList editList = new EditList();
        editList.mutationStateInput.setValue(mutationState, editList);
        mutationState.initByName("tree", tree, "data", data);
        
        ParsimonyMutationStateInitialiser init = new ParsimonyMutationStateInitialiser();
        init.initByName("mutationState", mutationState, "tree", tree);
        init.initStateNodes();
        
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "0.3 0.2 0.2 0.3");

        GTR gtr = new GTR();
        gtr.initByName("frequencies", freqs);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", gtr);

        MutationStateTreeLikelihood likelihood = new MutationStateTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "mutationState", mutationState, "editList", editList);
        
        double logP = likelihood.calculateLogP();

        
        State state = new State();
        state.initByName("stateNode", mutationState);
        state.initialise();
        state.setPosterior(likelihood);
        
        state.store(0);

        Randomizer.setSeed(127);
        
        MutationMover operator = new MutationMover();
        operator.initByName("mutationState", mutationState, "weight", 1.0);
        operator.proposal();
        state.checkCalculationNodesDirtiness();
        double logP2 = likelihood.calculateLogP();
        assertNotEquals(logP, logP2);
        
        state.restore();
        state.restoreCalculationNodes();
        double logP3 = likelihood.calculateLogP();
        
        assertEquals(logP, logP3, 1e-10);
	}
	
}
