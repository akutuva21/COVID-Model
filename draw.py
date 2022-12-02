from pyboolnet import file_exchange
from pyboolnet import state_transition_graphs as stgs
from pyboolnet.interaction_graphs import primes2igraph, igraph2image, add_style_interactionsigns

primes = file_exchange.bnet2primes("BN.bnet")
# print(file_exchange.primes2bnet(primes))
igraph = primes2igraph(primes)
igraph.graph["splines"] = "curved"
igraph.graph["label"] = "BN Model!"
for x in igraph.nodes:
    if "GF" in x:
        x["shape"] = "square"
        x["fillcolor"] = "lightblue"
add_style_interactionsigns(igraph)
for engine in ["dot", "neato", "fdp", "sfdp", "circo", "twopi"]:
    try:
        igraph2image(igraph, f'bnet_{engine}.png', layout_engine = engine)
    except:
        pass

with open(r"BN.bnet", 'r') as fp:
    lines = fp.readlines()
all_components = [line.split(',')[0] for line in lines]

exclude = ["Virus", "ACE2"]
initial = "".join([str(0) if x not in exclude else str(1)
                  for x in all_components])
print(initial)

update = "synchronous"
stg = stgs.primes2stg(primes, update, initial_states=initial)
for engine in ["dot", "twopi"]:
    try:
        stgs.stg2image(stg, f'stg_{engine}.pdf', layout_engine=engine)
    except:
        pass

'''total = 5
activities = [initial] + ["".join([str(np.random.choice([0, 1]))
                                   for x in np.arange(len(all_components))])
                          for y in np.arange(total - 1)]
# IG.activities2animation(igraph, activities, "animation.gif",
# delay = 60, loop = 3)

steady_states, cyclic_attractors = attractors.compute_attractors_tarjan(stg)
print(steady_states, cyclic_attractors)'''