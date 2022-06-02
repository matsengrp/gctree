import gctree
import random
import gctree.phylip_parse as pp
import click
from collections import Counter
import ete3
import historydag as hdag
import pickle
import numpy as np

_yey = [
    [0, 1, 1, 1, 1],
    [1, 0, 1, 1, 1],
    [1, 1, 0, 1, 1],
    [1, 1, 1, 0, 1],
    [1, 1, 1, 1, 0],
]

delta_vectors = {
    code: np.array(
        [
            0 if base in gctree.utils.ambiguous_dna_values[code] else 1
            for base in gctree.utils.bases
        ]
    )
    for code in gctree.utils.ambiguous_dna_values
}
code_vectors = {
    code: np.array(
        [
            0 if base in gctree.utils.ambiguous_dna_values[code] else float("inf")
            for base in gctree.utils.bases
        ]
    )
    for code in gctree.utils.ambiguous_dna_values
}
cost_adjust = {
    base: np.array([int(not i == j) for j in range(5)])
    for i, base in enumerate(gctree.utils.bases)
}


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def _cli():
    """
    A collection of tools for calculating parsimony scores of newick trees, and
    using them to create a history DAG
    """
    pass


def sankoff_upward(tree, gap_as_char=False):
    """Compute Sankoff cost vectors at nodes in a postorder traversal,
    and return best possible parsimony score of the tree.

    Args:
        gap_as_char: if True, the gap character ``-`` will be treated as a fifth character. Otherwise,
            it will be treated the same as an ``N``."""
    if gap_as_char:

        def translate_base(char):
            return char

    else:

        def translate_base(char):
            if char == "-":
                return "N"
            else:
                return char

    seq_len = len(tree.sequence)
    adj_arr = np.array([_yey] * seq_len)
    # First pass of Sankoff: compute cost vectors
    for node in tree.traverse(strategy="postorder"):
        node.add_feature(
            "cv",
            np.array(
                [code_vectors[translate_base(base)].copy() for base in node.sequence]
            ),
        )
        if not node.is_leaf():
            child_costs = []
            for child in node.children:
                stacked_child_cv = np.stack((child.cv,) * 5, axis=1)
                total_cost = adj_arr + stacked_child_cv
                child_costs.append(np.min(total_cost, axis=2))
            child_cost = np.sum(child_costs, axis=0)
            node.cv = (
                np.array([code_vectors[translate_base(base)] for base in node.sequence])
                + child_cost
            )
    return np.sum(np.min(tree.cv, axis=1))


def disambiguate(tree, random_state=None, remove_cvs=False, adj_dist=False):
    """Randomly resolve ambiguous bases using a two-pass Sankoff Algorithm on
    subtrees of consecutive ambiguity codes."""
    if random_state is None:
        random.seed(tree.write(format=1))
    else:
        random.setstate(random_state)

    seq_len = len(tree.sequence)
    adj_arr = np.array([_yey] * seq_len)
    sankoff_upward(tree)
    # Second pass of Sankoff: choose bases
    preorder = list(tree.traverse(strategy="preorder"))
    for node in preorder:
        new_seq = []
        base_indices = []
        for idx in range(seq_len):
            min_cost = min(node.cv[idx])
            base_index = random.choice(
                [i for i, val in enumerate(node.cv[idx]) if val == min_cost]
            )
            base_indices.append(base_index)

        # Adjust child cost vectors
        adj_vec = adj_arr[np.arange(seq_len), base_indices]
        for child in node.children:
            child.cv += adj_vec
        new_seq = [gctree.utils.bases[base_index] for base_index in base_indices]
        node.sequence = "".join(new_seq)

    if remove_cvs:
        for node in tree.traverse():
            try:
                node.del_feature("cv")
            except (AttributeError, KeyError):
                pass
    if adj_dist:
        tree.dist = 0
        for node in tree.iter_descendants():
            node.dist = gctree.utils.hamming_distance(node.up.sequence, node.sequence)
    return tree


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
                    fasta_map[seqid] += line.strip().upper()
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
    ambig_seq = "?" * seq_len
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


def parsimony_scores_from_topologies(newicks, fasta_map, gap_as_char=False, **kwargs):
    """returns a generator on parsimony scores of trees specified by newick strings and fasta.
    additional keyword arguments are passed to `build_tree`."""
    # eliminate characters for which there's no diversity:
    informative_sites = [
        idx for idx, chars in enumerate(zip(*fasta_map.values())) if len(set(chars)) > 1
    ]
    newfasta = {
        key: "".join(oldseq[idx] for idx in informative_sites)
        for key, oldseq in fasta_map.items()
    }
    yield from (
        sankoff_upward(build_tree(newick, newfasta, **kwargs), gap_as_char=gap_as_char)
        for newick in newicks
    )


def parsimony_scores_from_files(treefiles, fastafile, **kwargs):
    """returns the parsimony scores of trees specified by newick files and a fasta file.
    Arguments match `build_trees_from_files`. Additional keyword arguments are passed to
    ``parsimony_scores_from_topologies``."""

    def load_newick(npath):
        with open(npath, "r") as fh:
            return fh.read()

    return parsimony_scores_from_topologies(
        (load_newick(path) for path in treefiles), load_fasta(fastafile), **kwargs
    )


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
@click.option(
    "-g",
    "--gap-as-char",
    is_flag=True,
    help="Treat gap character `-` as a fifth character. Otherwise treated as ambiguous `N`.",
)
def _cli_parsimony_score_from_files(
    treefiles,
    fasta_file,
    root_id,
    newick_format,
    include_internal_sequences,
    save_to_dag,
    gap_as_char,
):
    """Print the parsimony score of one or more newick files"""
    pscore_gen = parsimony_scores_from_files(
        treefiles,
        fasta_file,
        reference_id=root_id,
        gap_as_char=gap_as_char,
        ignore_internal_sequences=(not include_internal_sequences),
    )

    for score, treepath in zip(pscore_gen, treefiles):
        print(treepath)
        print(score)


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
