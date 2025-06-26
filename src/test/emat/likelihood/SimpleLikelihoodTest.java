package test.emat.likelihood;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GTR;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beastfx.app.inputeditor.XMLImporter;
import emat.likelihood.AncestralStateMutationStateInitialiser;
import emat.likelihood.EditList;
import emat.likelihood.MutationState;
import emat.likelihood.MutationStateTreeLikelihood;
import emat.likelihood.ParsimonyMutationStateInitialiser;
import emat.util.MutationTreeWithMetaDataLogger;

public class SimpleLikelihoodTest {
	
	@Before
    public void setUp() {
		Randomizer.setSeed(1273);
    }


	
	@Test
	public void testLikelihoodInitByParsimony() {
		testLikelihood(0);
	}

	@Test
	public void testLikelihoodInitByAncestralReconstruction() {
		testLikelihood(1);
	}
	
	public void testLikelihood(int initMode) {
		
        Sequence a = new Sequence("A", "A A C G T TT -");
        Sequence b = new Sequence("B", "A C C G T CC -");
        Sequence c = new Sequence("C", "A A C G T TT A");

        a = new Sequence("A", "T-");
        b = new Sequence("B", "C-");
        c = new Sequence("C", "TA");

        
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
                
        
        System.out.println("root states: " + Arrays.toString(mutationState.getRootStateFreqs()));;
//        System.out.println("state lengths: " + Arrays.toString(mutationState.getTotalStateLengths()));;
//        System.out.println("mutation counts: " + Arrays.toString(mutationState.getMutationCounts()));;
        
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "0.25 0.25 0.25 0.25");

        GTR gtr = new GTR();
        gtr.initByName("frequencies", freqs);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", gtr);
        
        RealParameter clockRate = new RealParameter("1.0");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", clockRate);

        if (initMode > 0) {
        	AncestralStateMutationStateInitialiser init = new AncestralStateMutationStateInitialiser();
        	init.initByName("mutationState", mutationState, "tree", tree, "data", data, "siteModel", siteModel, "branchRateModel", clockModel);
	        init.initStateNodes();
        } else {
	        ParsimonyMutationStateInitialiser init = new ParsimonyMutationStateInitialiser();
	        init.initByName("mutationState", mutationState, "tree", tree);
	        init.initStateNodes();
        }
        

        MutationStateTreeLikelihood likelihood = new MutationStateTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "mutationState", mutationState, "editList", editList);
        
        double logP = likelihood.calculateLogP();
        System.out.println(logP);
   	}
	
		
	final static String dnaNewickTree = "((((Carp:0.11106063102912841,Loach:0.11106063102912841):0.06374314048186747,Frog:0.17480377151099588):0.04959439049312739,Chicken:0.22439816200412327):0.03095577892707152,((((Cow:0.10545881932623141,Whale:0.10545881932623141):0.01983872729398753,Seal:0.12529754662021894):0.037258475892112036,(Mouse:0.08708082642962245,Rat:0.08708082642962245):0.07547519608270853):0.05755867568285161,Human:0.2201146981951826):0.0352392427360122):0.0";
	final static String dnaFile = "examples/dna.xml";
	
	@Test
	public void testDNALikelihoodInitByParsimony() throws FileNotFoundException {
		testYangLikelihood(0, dnaNewickTree, dnaFile);
	}

	@Test
	public void testDNALikelihoodInitByAncestralReconstruction() throws FileNotFoundException {
		testYangLikelihood(1, dnaNewickTree, dnaFile);
	}
	
	final static String yangNewickTree = "((((((Balaenoptera_physalus:0.1317614013511199,Physeter_macrocephalus:0.29729950142864925):0.1978230084935262,Hippopotamus_amphibius:0.404548910412583):0.07915994982945396,Bos_tarus:0.3203021342425879):0.09295186380834797,((Canis_familiaris:0.368143707788285,Ursus_americanus:0.2818638512656598):0.09623390798084952,Felis_catus:0.3443774925057492):0.11166286579227092):0.03084971467984743,(Equus_caballus:0.25928830507424405,Rhinoceros_unicornis:0.2450943168319102):0.14481535193085104):0.08688257364452157,((Callithrix_jacchus:0.7708442409231275,(((Gorilla_gorilla:0.11200121192857693,(Homo_sapiens:0.10021261041786533,Pan_troglodytes:0.0866150994855941):0.042822465692259115):0.05285138771394582,Pongo_pygmaeus:0.18281384178550608):0.14775954779802758,Macaca_mulatta:0.4541320335152548):0.4431781922929532):0.3385390967282038,(((((Cheirogaleus_medius:0.18941226567755987,(((((((Microcebus_berthae:0.01938597077880111,Microcebus_myoxinus:0.024858344372513885):0.0028383800744640197,Microcebus_rufus_1:0.015618025459872342):0.014859461945712749,(Microcebus_rufus_2:0.04016010534297809,Microcebus_sambiranensis:0.0430912943072993):0.005229480435073985):0.008019535757913365,Microcebus_tavaratra:0.05121426102695925):0.029691694156483095,Microcebus_ravelobensis:0.0735982604998342):0.018398350139204322,(Microcebus_griseorufus:0.06958823337701125,Microcebus_murinus:0.07342351488381116):0.029713473243328514):0.11579429366721261,Mirza_coquereli:0.1787874439372752):0.0949029505126372):0.09606565799470823,Lepilemur_edwardsi:0.3762328786839756):0.028804849353506135,(((Eulemur_mongoz:0.18786994118752498,(Hapalemur_griseus:0.12470197573006947,Lemur_catta:0.088864706112877):0.049366376882922114):0.03341352279016996,Varecia_variegata:0.274397812811469):0.07633796453539599,Propithecus_tattersalli:0.22935575994044355):0.032873353084618806):0.17476939704063443,(Galago_crassicaudatus:0.2868664380903533,Loris_tardigradus:0.3195335325695139):0.23956891421891102):0.05345144243283673,Daubentonia_madagascariensis:0.4861670375495506):0.11584851784089789):0.052778202133946905):0.0";
	final static String yangFile = "examples/yang.xml";

	
	@Test
	public void testYangLikelihoodInitByParsimony() throws FileNotFoundException {
		testYangLikelihood(0, yangNewickTree, yangFile);
	}

	@Test
	public void testYangLikelihoodInitByAncestralReconstruction() throws FileNotFoundException {
		testYangLikelihood(1, yangNewickTree, yangFile);
	}
	
	public void testYangLikelihood(int initMode, String newick, String file) throws FileNotFoundException {
		XMLImporter importer = new XMLImporter();
		List<BEASTInterface> list = importer.loadFile(new File(file));
        Alignment data = (Alignment) list.get(0);

		TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", newick,
                "IsLabelledNewick", true);
        
        MutationState mutationState = new MutationState();
        EditList editList = new EditList();
        editList.mutationStateInput.setValue(mutationState, editList);
        mutationState.initByName("tree", tree, "data", data);
        
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "0.3632	0.3628	0.0629	0.2111");

        GTR gtr = new GTR();
        gtr.initByName("frequencies", freqs, "rateAC", "0.0394", "rateAG", "1.1736", "rateAT", "0.0897", "rateCG", "0.0817", "rateCT", "0.1151");

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", gtr);

        RealParameter clockRate = new RealParameter("1.0");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", clockRate);

        if (initMode > 0) {
            AncestralStateMutationStateInitialiser init = new AncestralStateMutationStateInitialiser();
        	init.initByName("mutationState", mutationState, "tree", tree, "data", data, "siteModel", siteModel, "branchRateModel", clockModel);
	        init.initStateNodes();
        } else {
            ParsimonyMutationStateInitialiser init = new ParsimonyMutationStateInitialiser();
            init.initByName("mutationState", mutationState, "tree", tree);
            init.initStateNodes();
        }

        MutationTreeWithMetaDataLogger logger = new MutationTreeWithMetaDataLogger();
        logger.initByName("mutationState", mutationState, "tree", tree);

        PrintStream out = new PrintStream("/tmp/yang.tree");
        logger.init(out);
        logger.log(0, out);
        logger.close(out);
   	}

	
	public static void main(String[] args) throws FileNotFoundException {
		new SimpleLikelihoodTest().testDNALikelihoodInitByParsimony();
	}
	
}
