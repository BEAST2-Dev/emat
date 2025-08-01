package emat.util;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import emat.likelihood.MutationOnBranch;
import emat.likelihood.MutationState;

@Description("Logs mutation information on trees + the usual metadata")
public class MutationTreeWithMetaDataLogger extends BEASTObject implements Loggable {
	final public Input<MutationState> stateInput = new Input<>("mutationState", "mutation state for the tree to be initialised", Validate.REQUIRED);

    final public Input<Tree> treeInput = new Input<>("tree", "tree to be logged", Validate.REQUIRED);
    // TODO: make this input a list of valuables
    final public Input<List<Function>> parameterInput = new Input<>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());
    final public Input<BranchRateModel.Base> clockModelInput = new Input<>("branchratemodel", "rate to be logged with branches of the tree");
    final public Input<Boolean> substitutionsInput = new Input<>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    final public Input<Integer> decimalPlacesInput = new Input<>("dp", "the number of decimal places to use writing branch lengths, rates and real-valued metadata, use -1 for full precision (default = full precision)", -1);
    final public Input<Boolean> sortTreeInput = new Input<>("sort", "whether to sort the tree before logging.", true);
    final public Input<Boolean> reportMutationsInput = new Input<>("reportMutations", "report individual mutations on branches", false);


    private boolean substitutions = false;

    private DecimalFormat df;
    private boolean sortTree;

    private MutationState state;
	
    @Override
    public void initAndValidate() {
		state = stateInput.get();

		int dp = decimalPlacesInput.get();
        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }

    	// without substitution model, reporting substitutions == reporting branch lengths 
        if (clockModelInput.get() != null) {
        	substitutions = substitutionsInput.get();
        }

        // default is to sort the tree
        sortTree = sortTreeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        treeInput.get().init(out);
    }

    @Override
    public void log(long sample, PrintStream out) {
        // make sure we get the current version of the inputs
        Tree tree = (Tree) treeInput.get().getCurrent();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode) {
        		metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        	}
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + sample + " = ");

        if (sortTree) {
            tree.getRoot().sort();
        }

        out.print(toNewick(tree.getRoot(), metadata, branchRateModel));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
    }

    /**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }

    String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList, branchRateModel));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList, branchRateModel));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }
		StringBuffer buf2 = new StringBuffer();
		
		buf2.append("[&");
		List<MutationOnBranch> mutations = state.getMutationList(node.getNr());
		buf2.append("mutationcount=" + mutations.size());

		Set<Integer> sitesWithMutations = new HashSet<>();
		for (MutationOnBranch m : mutations) {
			sitesWithMutations.add(m.siteNr());
		}
		int multiSiteMutations = mutations.size() - sitesWithMutations.size();
		buf2.append(",multiSiteMutations=" + multiSiteMutations);
		
		if (mutations.size() > 0 && reportMutationsInput.get()) {
			buf2.append(",mutations=\"");
			for (int i = 0; i < mutations.size()-1; i++) {
				buf2.append(mutations.get(i).toString() + ",");
			}
			buf2.append(mutations.get(mutations.size()-1).toString() + "\"");
		}
		
		if (metadataList.size() > 0) {
			for (Function metadata : metadataList) {
				if (metadata instanceof Parameter<?>) {
					Parameter<?> p = (Parameter<?>) metadata;
					int dim = p.getMinorDimension1();
					if (p.getMinorDimension2() > node.getNr()) {
						buf2.append(",");
						buf2.append(((BEASTObject) metadata).getID());
						buf2.append('=');
						if (dim > 1) {
							buf2.append('{');
							for (int i = 0; i < dim; i++) {
								if (metadata instanceof RealParameter) {
									RealParameter rp = (RealParameter) metadata;
									appendDouble(buf2, rp.getMatrixValue(node.getNr(), i));
								} else {
									buf2.append(p.getMatrixValue(node.getNr(), i));
								}
								if (i < dim - 1) {
									buf2.append(',');
								}
							}
							buf2.append('}');
						} else {
							if (metadata instanceof RealParameter) {
								RealParameter rp = (RealParameter) metadata;
								appendDouble(buf2, rp.getArrayValue(node.getNr()));
							} else {
								buf2.append(metadata.getArrayValue(node.getNr()));
							}
						}
					} else {
					
					}
				} else {
					if (metadata.getDimension() > node.getNr()) {
						buf2.append(",");
						buf2.append(((BEASTObject) metadata).getID());
						buf2.append('=');
						buf2.append(metadata.getArrayValue(node.getNr()));
					}
				}
			}
			if (buf2.length() > 2 && branchRateModel != null) {
				buf2.append(",");
			}
		}
		if (branchRateModel != null) {
			buf2.append("rate=");
			appendDouble(buf2, branchRateModel.getRateForBranch(node));
		}
		buf2.append(']');
		if (buf2.length() > 3) {
			buf.append(buf2.toString());
		}
        buf.append(":");
        if (substitutions) {
            appendDouble(buf, node.getLength() * branchRateModel.getRateForBranch(node));
        } else {
            appendDouble(buf, node.getLength());
        }
        return buf.toString();
    }


    @Override
    public void close(PrintStream out) {
        treeInput.get().close(out);
    }

}

 