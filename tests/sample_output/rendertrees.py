import pickle

with open("littleforest.p", "rb") as fh:
    forest = pickle.load(fh)

forest.forest = tuple(
    sorted(
        forest.forest,
        key=lambda x: -x.ll(0.4961832069459249, 0.3658705215787903, build_cache=False)[
            0
        ],
    )
)

for index, ctree in enumerate(forest.forest):
    ctree.render(outfile=f"gctree.out.{index}.tree.svg", idlabel=True)
