import gctree
import gctree.phylip_parse as pp
import click
from collections import Counter
import ete3
import historydag as hdag
import pickle


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def _cli():
    """
    A collection of tools for calculating parsimony scores of newick trees, and
    using them to create a history DAG
    """
    pass


def load_fasta(fastapath):
    """Load a fasta file as a dictionary, with sequence ids as keys and sequences as values."""
    fasta_map = {}
    with open(fastapath, "r") as fh:
        seqid = None
        for line in fh:
            if line[0] == ">":
                seqid = line[1:].strip()
                if seqid in fasta_map:
                    raise ValueError(
                        "Duplicate records with matching identifier in fasta file"
                    )
                else:
                    fasta_map[seqid] = ""
            else:
                if seqid is None and line.strip():
                    raise ValueError(
                        "First non-blank line in fasta does not contain identifier"
                    )
                else:
                    fasta_map[seqid] += line.strip()
    return fasta_map


def build_tree(
    newickstring,
    fasta_map,
    newickformat=1,
    reference_id=None,
    reference_sequence=None,
    ignore_internal_sequences=False,
):
    """Load an ete tree from a newick string, and add 'sequence' attributes from fasta.

    If internal node sequences aren't specified in the newick string and fasta data,
    internal node sequences will be fully ambiguous (contain repeated N's).

    Arguments:
        newickstring: a newick string
        fasta_map: a dictionary with sequence id keys matching node names in `newickstring`, and sequence values.
        newickformat: the ete format identifier for the passed newick string. See ete docs
        reference_id: key in `fasta_map` corresponding to the root sequence of the tree, if root sequence is fixed.
        reference_sequence: fixed root sequence of the tree, if root sequence is fixed
        ignore_internal_sequences: if True, sequences at non-leaf nodes specified in fasta_map and newickstring will be ignored, and all internal sequences will be fully ambiguous."""
    tree = ete3.Tree(newickstring, format=newickformat)
    # all fasta entries should be same length
    seq_len = len(next(iter(fasta_map.values())))
    ambig_seq = "N" * seq_len
    for node in tree.traverse():
        if node.is_root() and reference_sequence is not None:
            node.add_feature("sequence", reference_sequence)
        elif node.is_root() and reference_id is not None:
            node.add_feature("sequence", fasta_map[reference_id])
        elif (not node.is_leaf()) and ignore_internal_sequences:
            node.add_feature("sequence", ambig_seq)
        elif node.name in fasta_map:
            node.add_feature("sequence", fasta_map[node.name])
        else:
            node.add_feature("sequence", ambig_seq)
    return tree


def build_trees_from_files(newickfiles, fastafile, **kwargs):
    """Same as `build_tree`, but takes a list of filenames containing newick strings, and a filename for a fasta file, and returns a generator on trees"""
    fasta_map = load_fasta(fastafile)
    trees = []
    for newick in newickfiles:
        with open(newick, "r") as fh:
            newick = fh.read()
        yield build_tree(newick, fasta_map, **kwargs)


def parsimony_score(tree):
    """returns the parsimony score of a (disambiguated) ete tree.
    Tree must have 'sequence' attributes on all nodes."""
    return sum(
        gctree.utils.hamming_distance(node.up.sequence, node.sequence)
        for node in tree.iter_descendants()
    )


def parsimony_scores_from_files(*args, **kwargs):
    """returns the parsimony scores of trees specified by newick files and a fasta file.
    Arguments match `build_trees_from_files`."""
    trees = build_trees_from_files(*args, **kwargs)
    trees = [pp.disambiguate(tree) for tree in trees]
    return [parsimony_score(tree) for tree in trees]


def build_dag_from_trees(trees):
    """Build a history DAG from trees containing a `sequence` attribute on all nodes.
    unifurcations in the provided trees will be deleted."""
    trees = [tree.copy() for tree in trees]
    for tree in trees:
        for node in tree.iter_descendants():
            if len(node.children) == 1:
                node.delete(prevent_nondicotomic=False)
        if len(tree.children) == 1:
            newchild = tree.add_child()
            newchild.add_feature("sequence", tree.sequence)
    return hdag.history_dag_from_etes(
        trees,
        ["sequence"],
    )


@_cli.command("parsimony_scores")
@click.argument("treefiles", nargs=-1)
@click.option(
    "-f",
    "--fasta-file",
    required=True,
    help="Filename of a fasta file containing sequences appearing on nodes of newick tree",
)
@click.option(
    "-r",
    "--root-id",
    default=None,
    help="The fasta identifier of the fixed root of provided trees. May be omitted if there is no fixed root sequence.",
)
@click.option(
    "-F",
    "--newick-format",
    default=1,
    help="Newick format of the provided newick file. See http://etetoolkit.org/docs/latest/reference/reference_tree.html#ete3.TreeNode",
)
@click.option(
    "-i",
    "--include-internal-sequences",
    is_flag=True,
    help="include non-leaf node labels, and associated sequences in the fasta file.",
)
@click.option(
    "-d",
    "--save-to-dag",
    default=None,
    help="Combine loaded and disambiguated trees into a history DAG, and save pickled DAG to provided path.",
)
def _cli_parsimony_score_from_files(
    treefiles,
    fasta_file,
    root_id,
    newick_format,
    include_internal_sequences,
    save_to_dag,
):
    """Print the parsimony score of one or more newick files"""
    fasta_map = load_fasta(fasta_file)
    parsimony_counter = Counter()
    trees = []
    for tree, treepath in zip(
        build_trees_from_files(
            treefiles,
            fasta_file,
            reference_id=root_id,
            ignore_internal_sequences=(not include_internal_sequences),
        ),
        treefiles,
    ):
        print(treepath)
        print(parsimony_score(pp.disambiguate(tree)))
        if save_to_dag is not None:
            trees.append(tree)
    if save_to_dag is not None:
        with open(save_to_dag, "wb") as fh:
            fh.write(pickle.dumps(build_dag_from_trees(trees)))


def summarize_dag(dag):
    """print summary information about the provided history DAG."""
    print("DAG contains")
    print("trees: ", dag.count_trees())
    print("nodes: ", len(list(dag.preorder())))
    print("edges: ", sum(len(list(node.children())) for node in dag.preorder()))
    print("parsimony scores: ", dag.weight_count())


@_cli.command("summarize-dag")
@click.argument("dagpath")
def _cli_summarize_dag(dagpath):
    """print summary information about the provided history DAG."""
    with open(dagpath, "rb") as fh:
        dag = pickle.load(fh)
    summarize_dag(dag)


if __name__ == "__main__":
    _cli()
