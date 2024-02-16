import networkx as nx
import matplotlib.pyplot as plt
import re, inspect

# Define the lambda functions
operations = {
        "ACE2": lambda: x.__setitem__("ACE2", OR(NOT(x["Virus"]), x["FOXO3A"], NOT(x["ADAM_17"]))), # added adam-17
        "ADAM_17": lambda: x.__setitem__("ADAM_17", OR(x["ANG_2_T1R"], x["HIF_1a"])), # added hif-1a
        "AKT": lambda: x.__setitem__("AKT", OR(x["mTORC2"], x["PI3K"], x["FOXO3A"])))
    }

G = nx.DiGraph()

for node, operation in operations.items():
    G.add_node(node)
    
    func_str = inspect.getsource(operation)
    
    pattern = r'\["(.*?)"\]'
    pointed_nodes = re.findall(pattern, func_str)
    print([match for match in pointed_nodes])
    
    for pointed_node in pointed_nodes:
        G.add_edge(node, pointed_node)

nx.draw(G, with_labels=True)
plt.show()