#! /usr/bin/env python
# -*- coding: utf-8 -*-


def try_cast(k):
    try:
        return int(k)
    except:
        pass
    try:
        return float(k)
    except ValueError:
        return k


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
    from ete3 import TreeStyle
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
    from ete3 import NodeStyle
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


def tree_render_user(tree, frame=None, tree_features=None, **kwargs):
    '''
    Base function to add rendering attributes to an ete3 tree.
    The function can plot tree is default mode and/or take user input.
    '''
    from ete3 import NodeStyle, faces, TextFace
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

        node.set_style(nstyle)

    # Add TextFace:
    def my_layout(node):
        N = None
        # Add default values:
        if not no_default and node.frequency > 1 and not hasattr(node, 'cl'):
            N = TextFace(node.frequency, fsize=12, fgcolor='black')  # Default: use node frequency as TextFace

        # If user specified labels are specified add these:
        if 'node_label' in tree_features:
            for t in tree_features['node_label']:
                node_attr = ifhasthenget(node, t)
                if node_attr and to_string(node_attr):
                    textface = to_string(node_attr)
                    N = TextFace(textface, fsize=12, fgcolor='black')
                    if 'e-' in textface and to_float(textface):  # Stupid exception because ete3 rendering doesn't understand scientific notation
                        N = TextFace('%2.1e    ' % to_float(textface), fsize=12, fgcolor='black')

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


def collapse_syn(tree):
    pass


def make_tree(tree, outfile, tree_features_file=None, statsfile=None, frame=None):
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
    if 'prune_branch_length' in tree_features and to_float(tree_features['prune_branch_length']):
        tree = prune_long_branches(tree, to_float(tree_features['prune_max_branch_length']))

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

    #### Collapse synonymous mutations upwards into single nodes:
    if 'collapse_syn' in tree_features:
        tree = collapse_syn(tree)
        # Implemnt this soon

    tree, ts = tree_render_user(tree, frame=frame, tree_features=tree_features)
    tree.render(outfile, tree_style=ts)



