package test.emat.likelihood;

import java.util.Arrays;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GTR;
import beast.base.evolution.tree.TreeParser;
import emat.likelihood.AncestralStateMutationStateInitialiser;
import emat.likelihood.EditList;
import emat.likelihood.MutationState;
import emat.likelihood.MutationStateTreeLikelihood;
import emat.likelihood.ParsimonyMutationStateInitialiser;

public class SimpleLikelihoodTest {

	
	@Test
	public void testLikelihood() {
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
        
//        ParsimonyMutationStateInitialiser init = new ParsimonyMutationStateInitialiser();
//        init.initByName("mutationState", mutationState, "tree", tree);
//        init.initStateNodes();
        
        
        
        System.out.println("root states: " + Arrays.toString(mutationState.getRootStateFreqs()));;
//        System.out.println("state lengths: " + Arrays.toString(mutationState.getTotalStateLengths()));;
//        System.out.println("mutation counts: " + Arrays.toString(mutationState.getMutationCounts()));;
        
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "0.25 0.25 0.25 0.25");

        GTR gtr = new GTR();
        gtr.initByName("frequencies", freqs);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", gtr);

        AncestralStateMutationStateInitialiser init = new AncestralStateMutationStateInitialiser();
        init.initByName("mutationState", mutationState, "tree", tree, "data", data, "siteModel", siteModel);
        init.initStateNodes();

        MutationStateTreeLikelihood likelihood = new MutationStateTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "mutationState", mutationState, "editList", editList);
        
        double logP = likelihood.calculateLogP();
        System.out.println(logP);
   	}
}
