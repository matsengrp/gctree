#! /usr/bin/env python
# -*- coding: utf-8 -*-


from ete3 import Tree, TreeNode, NodeStyle, TreeStyle, faces, TextFace, add_face_to_node, CircleFace, AttrFace, CircleFace, PieChartFace


def try_cast(k):
    try:
        return int(k)
    except:
        pass
    try:
        return float(k)
    except ValueError:
        return k


def is_number(k):
    try:
        int(k)
        return True
    except:
        pass
    try:
        float(k)
        return True
    except ValueError:
        return False


def to_string(s):
    try:
        return str(s)
    except:
        return False


def to_float(s):
    try:
        return float(s)
    except:
        return False


def to_int(s):
    try:
        return int(s)
    except:
        return False


def to_int_or_float(s):
    try:
        return int(s)
    except:
        pass
    try:
        return float(s)
    except:
        return False


def ifhasthenget(obj, attr):
    try:
        return getattr(obj, attr)
    except:
        return False


def collapse_low_freq(tree, fk, cl_except):
    for node in tree.traverse():
        if node.is_leaf() is False and hasattr(node, 'cl') is False:
            for child in node.get_children():
                # Skip those node to keep:
                if has_true_attrlist(child, cl_except):
                    continue
                elif not keep_lineage(node, child, fk, cl_except):
                    for linchild in child.traverse():
                        linchild.add_features(cl=True)
                elif child.is_leaf() and fk[0] > child.frequency and fk[1] > 0:
                    child.add_features(cl=True)

    return tree


def has_true_attrlist(node, attrlist):
    if [1 for a in attrlist if ifhasthenget(node, a)]:
        return True
    else:
        return False


def keep_lineage(first_node, node, fk, cl_except):
    '''
    Search for a lineage decending from a node that is longer
    than a given threshold, max_linlen. If node exception is meeet
    in cl_except e.g. a feature that should never be collapsed
    the searchis terminated, likewise it is terminated when the
    threshold of a lineage is found, that is:
    linlen > max_linlen True

                         child
                      -----*
                     |
    first_node     node
        *------------*
                     |       child
                      ---------*
                    <-------------->
                      linlen includes 2 nodes
    '''
    if has_true_attrlist(node, cl_except):
        return True
    if node.is_leaf() and first_node is not node:
        f_max = 0
        k = 1
        while True:
            if node.frequency > f_max:
                f_max = node.frequency
            node = node.up
            if node is first_node:
                return [f_max, k]
            else:
                k += 1
    elif node.is_leaf() and first_node is node:
        return [node.frequency, 1]
    for child in node.get_children():
        fk_rec = keep_lineage(first_node, child, fk, cl_except)
        # Terminate the search because the max_linlen is surpased:
        if fk_rec is False:
            continue
        if fk_rec is True or fk_rec[0] > fk[0] or fk_rec[1] > fk[1]:
            return True
    return False


def read_plot_config(fnam):
    import json
    '''
    ### File format is JSON ###
    {
        "frequency": "LNobs", 
        "node_color": [
            [
                "hit_affinities", 
                "firebrick"
            ], 
            [
                "node_color", 
                "common_links"
            ]
        ], 
        "node_label": "hit_affinities", 
        "node_shape": [
            [
                "hit_affinities", 
                "square"
            ]
        ]
    }
    '''
    allowed_keys = ['node_color', 'node_shape', 'frequency', 'node_label', 'node_size', 'no_default', 'prune_max_branch_length', 'collapse_low_freq', 'collapse_syn']
    with open(fnam) as json_data:
        d = json.load(json_data)

    not_found = [k for k in d if k not in allowed_keys]
    if not_found:
        raise RuntimeError('Key not found in plot config: [{0}]\nPossible keys: [{1}]'
                           .format(", ".join(str(i) for i in not_found), ", ".join(str(i) for i in allowed_keys)))

    return d


def read_tree_stats(fnam):
    '''
    Read a comma separated table into a dict of dicts.
    First line will be interpretted as the header with column names,
    the following lines will be interpretted as entries with the first
    element being the node name to map back to.
    '''
    with open(fnam) as fh:
        lines = fh.readlines()
    lines[:] = [line.strip() for line in lines]
    features = lines[0].split(',')
    lines.pop(0)
    tstat = dict()
    for line in lines:
        cols = line.split(',')
        assert(cols[0] not in tstat)
        tstat[cols[0]] = {features[j+1]: try_cast(c) for j, c in enumerate(cols[1:])}
    return tstat


def hamming_distance(seq1, seq2):
    '''Hamming distance between two sequences of equal length'''
    return sum(x != y for x, y in zip(seq1, seq2))


def tree_render_minimum(tree):
    '''Set the minimum settings on a ete3 tree for rendering a GCtree.'''
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.rotation = 90
    tree.ladderize()

    return tree, ts


def tree_render_default(tree, frame=None):
    '''Set the default settings on a ete3 tree for rendering a GCtree.'''
    import sys
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    import scipy
    # Set the minimum requirements:
    tree, ts = tree_render_minimum(tree)
    # Then add the defaults:
    for node in tree.traverse():
        nstyle = NodeStyle()
        # Set node size according the allele frequency:
        small_node = 5
        if node.frequency == 0:
            nstyle['size'] = small_node
            nstyle['fgcolor'] = 'grey'
        else:
            nstyle['size'] = small_node*scipy.sqrt(node.frequency) + small_node
            nstyle['fgcolor'] = 'black'

        # Mark clades with stop codon and distinguish sense/missense mutations:
        if node.up is not None:
            if set(node.sequence.upper()) == set('ACGT'):
                if frame is not None:
                    aa = Seq(node.sequence[(frame-1):(frame-1+(3*(((len(node.sequence)-(frame-1))//3))))],
                             generic_dna).translate()
                    aa_parent = Seq(node.up.sequence[(frame-1):(frame-1+(3*(((len(node.sequence)-(frame-1))//3))))],
                                    generic_dna).translate()
                    nonsyn = hamming_distance(aa, aa_parent)
                    if '*' in aa:
                        nstyle['bgcolor'] = 'red'
                    if nonsyn > 0:
                        nstyle['hz_line_color'] = 'black'
                        nstyle['hz_line_width'] = nonsyn
                    else:
                        nstyle['hz_line_type'] = 1
        node.set_style(nstyle)

    return tree, ts


def tree_render_user(tree, frame=None, tree_features=None, namecolor=None):
    '''
    Base function to add rendering attributes to an ete3 tree.
    The function can plot tree is default mode and/or take user input.
    '''
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    import scipy

    if tree_features is None:  # No user defined tree features
        return tree_render_default(tree, frame=frame)
    elif 'no_default' in tree_features and tree_features['no_default']:
        tree, ts = tree_render_minimum(tree, frame=frame)
        no_default = True
    else:
        tree, ts = tree_render_default(tree, frame=frame)
        no_default = False

    for node in tree.traverse():
        nstyle = NodeStyle()

        small_node = 5
        if not no_default and node.frequency == 0 and not hasattr(node, 'cl'):
            nstyle['size'] = small_node
            nstyle['fgcolor'] = 'grey'
        elif not no_default and hasattr(node, 'cl'):
            nstyle['size'] = 2
            nstyle['fgcolor'] = 'mediumblue'
        else:
            # Assign node size:
            if 'node_size' in tree_features and ifhasthenget(node, tree_features['node_size']):
                nstyle['size'] = small_node*scipy.sqrt(getattr(node, tree_features['node_size'])) + small_node
                nstyle['fgcolor'] = 'black'  # Keep default color unless otherwise is specified below
            # Assign user defined colors:
            if 'node_color' in tree_features:
                for t in tree_features['node_color']:
                    node_attr = ifhasthenget(node, t[0])
                    if node_attr:
                        nstyle['fgcolor'] = t[1]
                        if nstyle['size'] < small_node:
                            nstyle['size'] = small_node
            if namecolor is not None and node.name in namecolor:
                # print(node.name)
                nstyle['fgcolor'] = 'DarkGreen'

        if 'node_shape' in tree_features:
            for t in tree_features['node_shape']:
                node_attr = ifhasthenget(node, t[0])
                if node_attr:
                    # Resize the area of the square to match a circle:
                    if t[1] == 'square':
                        nstyle['size'] = scipy.sqrt(scipy.pi) * nstyle['size'] / 2
                    nstyle['shape'] = t[1]

        # Mark clades with stop codon and distinguish sense/missense mutations:
        if node.up is not None:
            if set(node.sequence.upper()) == set('ACGT'):
                if frame is not None:
                    aa = Seq(node.sequence[(frame-1):(frame-1+(3*(((len(node.sequence)-(frame-1))//3))))],
                             generic_dna).translate()
                    aa_parent = Seq(node.up.sequence[(frame-1):(frame-1+(3*(((len(node.sequence)-(frame-1))//3))))],
                                    generic_dna).translate()
                    nonsyn = hamming_distance(aa, aa_parent)
                    if '*' in aa:
                        nstyle['bgcolor'] = 'red'
                    if nonsyn > 0:
                        nstyle['hz_line_color'] = 'black'
                        nstyle['hz_line_width'] = nonsyn
                    else:
                        nstyle['hz_line_type'] = 1

        node.set_style(nstyle)

    # Add TextFace:
    def my_layout(node):
        N = None
        # If user specified labels are specified add these:
        if 'node_label' in tree_features:
            for t in tree_features['node_label']:
                node_attr = ifhasthenget(node, t)
                if node_attr and to_string(node_attr):
                    if t == 'name' and hasattr(node, 'cl'):
                        continue
                    textface = to_string(node_attr)
                    N = TextFace(textface, fsize=12, fgcolor='black')
                    if 'e-' in textface and to_float(textface):  # Stupid exception because ete3 rendering doesn't understand scientific notation
                        N = TextFace('%2.1e    ' % to_float(textface), fsize=12, fgcolor='black')

        # Add default values:
        elif not no_default and node.frequency > 1 and not hasattr(node, 'cl'):
            N = TextFace(node.frequency, fsize=12, fgcolor='black')  # Default: use node frequency as TextFace

        if N is not None:
            N.rotation = -90
            faces.add_face_to_node(N, node, 0, position='branch-top')
    ts.layout_fn = my_layout

    return tree, ts


def prune_long_branches(tree, cutlength):
    '''
    Detach the whole decending clade if a branch length
    is encountered above a given cutlength.
    '''
    for node in tree.traverse():
        if node.dist > cutlength:
            node.detach()
    return tree


def update_after_collapse(node_up, node, tree_features=None):
    '''
    Update the feature according to its type and specifications
    in the tree_features dict
    '''
    feature_set = node_up.features.union(node.features)
    feature_set.remove('dist')
    feature_set.remove('name')
    for feature in feature_set:
        if hasattr(node, feature) and not hasattr(node_up, feature):
            setattr(node_up, feature, getattr(node, feature))

        elif hasattr(node, feature) and hasattr(node_up, feature) and getattr(node_up, feature) == '':
            setattr(node_up, feature, getattr(node, feature))
        elif hasattr(node, feature) and hasattr(node_up, feature) and getattr(node, feature) == '':
            pass
        elif hasattr(node, feature) and hasattr(node_up, feature):
            if is_number(getattr(node, feature)) and feature in tree_features['collapse_syn']:
                if tree_features['collapse_syn'][feature] == 'sum':
                    setattr(node_up, feature, (getattr(node, feature) + getattr(node_up, feature)))
                elif tree_features['collapse_syn'][feature] == 'mean':
                    setattr(node_up, feature, (getattr(node, feature) + getattr(node_up, feature))/2)
                else:
                    raise Exception('Do not know this option:', tree_features['collapse_syn'][feature])
            elif is_number(getattr(node, feature)):
                setattr(node_up, feature, (getattr(node, feature) + getattr(node_up, feature)))
            else:
                pass
    return node_up


def collapse_syn(tree, frame=None, tree_features=None):
    '''
    Recalculate the branch lengths as hamming distance between amino acid sequences.
    Then collapse all branches with zero length and update the features by either
    taking the sum of the mean, if the feature is a number, if the feature is a string
    then just keep the parent feature.
    '''
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    tree.dist = 0  # no branch above root
    for node in tree.iter_descendants():
        aa = Seq(node.sequence[(frame-1):(frame-1+(3*(((len(node.sequence)-(frame-1))//3))))],
                 generic_dna).translate()
        aa_parent = Seq(node.up.sequence[(frame-1):(frame-1+(3*(((len(node.sequence)-(frame-1))//3))))],
                        generic_dna).translate()
        node.dist = hamming_distance(aa, aa_parent)

    for node in tree.get_descendants():
        if node.dist == 0:
            if node.name and not node.up.name:
                node.up.name = node.name
            node.up = update_after_collapse(node.up, node, tree_features=tree_features)
            node.delete(prevent_nondicotomic=False)

    return tree


def make_tree(tree, outfile, tree_features_file=None, statsfile=None, frame=None, namecolor=None):
    try:
        import cPickle as pickle
    except:
        import pickle
    from gctree import CollapsedForest, CollapsedTree
    '''
    This function wraps the rendering of an ete3 tree with custom user options
    and the adding to the default GCtree rendering.
    '''
    if tree_features_file is not None and statsfile is not None:
        default_plot = False
    elif tree_features_file is None and statsfile is None:
        default_plot = True
    else:
        raise RuntimeError('Both tree_features_file and statsfile must be provided.')
    # No features, no fun, just render the default tree:
    if default_plot is True:
        tree, ts = tree_render_user(tree, frame=frame)
        tree.render(outfile, tree_style=ts)
        return

    # Read in a JSON file that specifies what to plot and where:
    tree_features = read_plot_config(tree_features_file)
    # Read in the new features and add them to the tree:
    tstat = read_tree_stats(statsfile)

    for node in tree.traverse():
        if node.name not in tstat:
            continue
        node.add_features(**tstat[node.name])

    # Prune of clades with very long branch lengths:
    if 'prune_max_branch_length' in tree_features and to_float(tree_features['prune_max_branch_length']):
        tree = prune_long_branches(tree, to_float(tree_features['prune_max_branch_length']))

    # Collapse synonymous mutations upwards into single nodes:
    if frame is not None and 'collapse_syn' in tree_features:
        tree = collapse_syn(tree, frame=frame, tree_features=tree_features)
    elif 'collapse_syn' in tree_features and frame is None:
        raise Exception('Cannot collapse synonymous mutations when frame parameter is not specified.')

    # Collapse tree nodes with low frequencies:
    if 'collapse_low_freq' in tree_features:
        try:
            f = int(tree_features['collapse_low_freq']['f'])
            k = int(tree_features['collapse_low_freq']['k'])
            if 'except' in tree_features['collapse_low_freq']:
                cl_except = list(tree_features['collapse_low_freq']['except'])
            else:
                cl_except = list()
        except Exception as e:
            raise RuntimeError('Could not cast f and k to integers, or maybe they are not defined at all?:', e)
        tree = collapse_low_freq(tree, [f, k], cl_except)

    tree, ts = tree_render_user(tree, frame=frame, tree_features=tree_features, namecolor=namecolor)
    tree.render(outfile, tree_style=ts)
    collapsed_tree = CollapsedTree(tree=tree)
    forest = CollapsedForest(forest=[collapsed_tree])
    with open(outfile[:-4] + '.p', 'wb') as f:
        pickle.dump(forest, f)


def main():
    import sys, os
    sys.path.append(os.path.abspath('/'.join(os.path.realpath(__file__).split('/')[:-2]) + "/bin"))
    import argparse
    import pickle
    from gctree import CollapsedTree

    parser = argparse.ArgumentParser(description='Custom plotting of a GCtree output tree.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--forrest_file", required=True, type=str, dest="forrest_file", help="Output pickled forrest file from GCtree.")
    parser.add_argument("--tree_numb", required=True, type=int, dest="tree_numb", metavar="INT", help="Index specifying which tree to pick from the forrest.")
    parser.add_argument("--outfile", required=True, type=str, dest="outfile", metavar="FILENAME", help="Output tree filename.")
    parser.add_argument("--frame", type=int, dest="frame", help="Which frame is the node sequences in? If unspecified frame is not enforced and stop codons are allowed.")
    parser.add_argument("--config", type=str, dest="config", metavar="JSON", help="JSON formatted file that specifies which features to plot and how.")
    parser.add_argument("--statsfile", type=str, dest="statsfile", help="Comman separated list of features linked to the node name.")
    args = parser.parse_args()

    with open(args.forrest_file) as forest_fh:
        forest = pickle.load(forest_fh)
    collapsed_tree = CollapsedTree(tree=forest.forest[args.tree_numb].tree, frame=args.frame)
    make_tree(collapsed_tree.tree, args.outfile, tree_features_file=args.config, statsfile=args.statsfile, frame=args.frame)


if __name__ == '__main__':
    main()
