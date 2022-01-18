**********************
Implementation Details
**********************

A note about abundances:
########################

In the original gctree package (v3.2.1 and before), an ambiguous tree from
dnapars, with edges collapsed and nodes named bijectively, was disambiguated in
one of the possible ways, and abundances were added **using node names**. This
meant that only leaves had nonzero abundance annotations. Internal nodes whose
disambiguated sequence matched a leaf sequence had abundance 0, since the names
of those nodes were not leaf names.

During collapsing, in CollapsedTree init, abundances were **added** up the
tree, which was correct, since only leaves had nonzero abundance.

Now, ambiguous dnapars trees will be put in a history DAG prior to
disambiguation. Disambiguation and collapsing occur inside the DAG, and
abundances will be added to each node's `attr` attribute dictionary according
to that node's **sequence**. Since the DAG is collapsed, only leaves and
leaf-adjacent nodes may have a leaf sequence label and hence a nonzero
abundance. In general, there's no guarantee that other isolated internal nodes
don't have a leaf sequence, but in the parsimony setting this isn't possible.
Just for safety though, nonzero abundances will only be added to leaves and
leaf adjacent nodes, by explicit filter.

With this new system, abundances can't be added up the tree during
CollapsedTree collapse. Instead they could be left alone (since leaf adjacent
nodes with a leaf sequence already have the correct abundance). In practice,
they will be transferred up the tree, to make the code slightly more resilient
to unexpected tree structures, or user-supplied trees in which only leaves have
correct abundance annotations. That is, if an edge is collapsed, it will be
assumed that if either node has nonzero abundance, it has the correct abundance
for that sequence, and the remaining node will carry the maximum abundance of
the original two nodes.
