from pyboolnet import interaction_graphs as IG
from pyboolnet import file_exchange
from pyboolnet import state_transition_graphs as stg
import numpy as np

primes = file_exchange.bnet2primes("BN.bnet")
igraph = IG.primes2igraph(primes)
igraph.graph["splines"] = "curved"
igraph.graph["label"] = "BN Model!"
IG.add_style_interactionsigns(igraph)
IG.igraph2image(igraph, "bnet.png")

state = stg.random_state(primes)
local_igraph = IG.local_igraph_of_state(primes, state)
IG.add_style_interactionsigns(local_igraph)
IG.igraph2image(local_igraph, "local_igraph.pdf")

with open(r"BN.bnet", 'r') as fp:
    lines = len(fp.readlines())
activities = ["".join([str(np.random.choice([0,1])) for x in range(0, lines)]) for y in np.arange(0, 100)]
IG.activities2animation(igraph, activities, "animation.gif")