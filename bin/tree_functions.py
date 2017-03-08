#! /usr/bin/env python
# -*- coding: utf-8 -*-


def lineage_length_bigger_than(first_node, node, max_linlen):
    '''
    Search for a lineage decending from a node that is longer
    than a given threshold, max_linlen. Terminate is such a
    lineage is found, that is:
    linlen > max_linlen True

                         child
                      -----*
                     |
    first_node     node
        *------------*
                     |       child
                      ---------*
                    <------------->
                      linlen includes 2 nodes
    '''
    if node.is_leaf() and first_node is not node:
        linlen = 1
        while True:
            node = node.up
            if node is first_node:
                return linlen
            else:
                linlen += 1
    elif node.is_leaf() and first_node is node:
        return 1
    for child in node.get_children():
        linlen = lineage_length_bigger_than(first_node, child, max_linlen)
        # Terminate the search because the max_linlen is surpased:
        if linlen is True or linlen > max_linlen:
            return True

    return False


def longest_lineage(first_node, node):
    if node.is_leaf() and node.up is not None and first_node is not node:
        linlen = 2  # Include both nodes e.g. *--*--* is a distance of 3
        while True:
            node = node.up
            if node is first_node:
                return linlen
            linlen += 1
    elif (first_node is node or node.up is None) and node.is_leaf():
        return 1
    max_linlen = 0
    for child in node.get_children():
        linlen = longest_lineage(first_node, child)
        if linlen > max_linlen:
            max_linlen = linlen
    return max_linlen


def sum_of_lineage_lengths(first_node, node):
    if node.is_leaf() and node.up is not None and first_node is not node:
        linlen = 2  # Include both nodes e.g. *--*--* is a distance of 3
        while True:
            node = node.up
            if node is first_node:
                return linlen
            linlen += 1
    elif (first_node is node or node.up is None) and node.is_leaf():
        return 1
    sum_linlen = 0
    for child in node.get_children():
        sum_linlen += sum_of_lineage_lengths(first_node, child)
    return sum_linlen
