import gctree
import gctree.phylip_parse as pp
import click
from collections import Counter
import ete3

@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def _cli():
    """
    A collection of tools for calculating parsimony scores of newick trees, and
    using them to create a history DAG
    """
    pass

def load_fasta(fastapath):
    fasta_map = {}
    with open(fastapath, 'r') as fh:
        seqid = None
        for line in fh:
            if line[0] == '>':
                seqid = line[1:].strip()
                if seqid in fasta_map:
                    raise ValueError("Duplicate records with matching identifier in fasta file")
                else:
                    fasta_map[seqid] = ""
            else:
                if seqid is None and line.strip():
                    raise ValueError("First non-blank line in fasta does not contain identifier")
                else:
                    fasta_map[seqid] += line.strip()
    return fasta_map

def build_tree(newickstring, fasta_map, newickformat=1, reference_id=None, reference_sequence=None, ignore_internal_sequences=False):
    tree = ete3.Tree(newickstring, format=newickformat)
    # all fasta entries should be same length
    seq_len = len(next(iter(fasta_map.values())))
    ambig_seq = 'N' * seq_len
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

def build_tree_from_files(newickfile, fastafile, **kwargs):
    with open(newickfile, 'r') as fh:
        newick = fh.read()
    return build_tree(newick, load_fasta(fastafile), **kwargs)

def parsimony_score(tree):
    return sum(gctree.utils.hamming_distance(node.up.sequence, node.sequence)
               for node in tree.iter_descendants())

def parsimony_score_from_files(*args, **kwargs):
    tree = build_tree_from_files(*args, **kwargs)
    tree = pp.disambiguate(tree)
    return parsimony_score(tree)

@_cli.command("parsimony_scores")
@click.argument("treefiles", nargs=-1)
@click.option("-f", "--fasta-file", required=True, help="Filename of a fasta file containing sequences appearing on nodes of newick tree")
@click.option("-r", "--root-id", default=None, help="The fasta identifier of the fixed root of provided trees. May be omitted if there is no fixed root sequence.")
@click.option("-F", "--newick-format", default=1, help="Newick format of the provided newick file. See http://etetoolkit.org/docs/latest/reference/reference_tree.html#ete3.TreeNode")
@click.option("-i", "--include-internal-sequences", is_flag=True, help=" non-leaf node labels, and associated sequences in the fasta file.")
def _cli_parsimony_score_from_files(treefiles, fasta_file, root_id, newick_format, include_internal_sequences):
    """Print the parsimony score of one or more newick files"""
    fasta_map = load_fasta(fasta_file)
    parsimony_counter = Counter()
    for treepath in treefiles:
        with open(treepath, 'r') as fh:
            tree = build_tree(fh.read(), fasta_map, newickformat=newick_format, reference_id=root_id, ignore_internal_sequences=(not include_internal_sequences))
        print(treepath)
        print(parsimony_score(pp.disambiguate(tree)))

if __name__ == "__main__":
    _cli()
