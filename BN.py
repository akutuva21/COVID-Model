import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns
import copy

def BN(components, n_iter=25, orders=500, n_states=3):
    '''Generates a boolean network based on the provided components

    Parameters
    ----------
    components : list
    n_iter : int, optional
        Number of iterations to run the network for (default is 25)
    orders : int, optional
        Number of samples to take for each iteration (default is 250)
        Averaged at the end of the function
    n_states : int
        Number of states for each component (default is 3)
    
    Raises
    ----------
        NotImplementedError
            If the component is not implemented
            If an 'if' statement was not executed (check names)
        ValueError
            If the number of orders or iterations is less than 1
            If the number of orders or iterations is not an integer
    '''

    # Generalized AND function for non-boolean states
    def AND(*args):
        # check that all args are positive
        for arg in args:
            if arg < 0:
                assert False, "AND function requires all arguments are positive"
        return min(args)

    # Generalized OR function for non-boolean states
    def OR(*args):
        # check that all args are positive
        for arg in args:
            if arg < 0:
                assert False, "AND function requires all arguments are positive"
        return sum(args) / len(args)
    
    # Generalized NOT function for non-boolean states
    def NOT(value, n=n_states-1):
        return n - value
    
    def g(x, index, update_func, n_states=n_states):
        current_state = x[index]
        update_func()
        new_state = x[index]

        if new_state > current_state:
            current_state += 1
        elif new_state < current_state:
            current_state -= 1

        if current_state < 0:
            current_state = 0
        elif current_state > (n_states - 1):
            current_state = n_states - 1

        return current_state

    operations = {
        "ACE2": lambda: x.__setitem__("ACE2", x["FOXO3A"] - OR(x["Virus"], x["ADAM_17"])), # added adam-17
        "ADAM_17": lambda: x.__setitem__("ADAM_17", OR(x["ANG_2_T1R"], x["HIF_1a"])), # added hif-1a
        "AKT": lambda: x.__setitem__("AKT", OR(x["mTORC2"], x["PI3K"], x["FOXO3A"])),
        # "AMPK": lambda: x.__setitem__("AMPK", ),
        "ANG_2": lambda: x.__setitem__("ANG_2", NOT(x["ACE2"])),
        "ANG_2_T1R": lambda: x.__setitem__("ANG_2_T1R", x["ANG_2"]),
        "Apoptosis": lambda: x.__setitem__("Apoptosis", OR(x["CASP8"], x["CASP9"])),
        # "Autophagy": lambda: x.__setitem__("Autophagy", NOT(x["mTORC1"])),
        "BCL_2": lambda: x.__setitem__("BCL_2", OR(x["NFKB"], x["CREB_1"], x["HIF_1a"]) - x["p53"]),
        # "Bax_Bak": lambda: x.__setitem__("Bax_Bak", NOT(x["BCL_2"])),
        "CASP1": lambda: x.__setitem__("CASP1", x["NLRP3"]),
        "CASP8": lambda: x.__setitem__("CASP8", OR(x["FADD"], x["p53"]) - OR(x["C_FLIP"], x["FOXO3A"])),
        # OR(AND(AND(x["FADD"], NOT(x["C_FLIP"])), NOT(x["FOXO3A"])), x["p53"])),
        # OR(x["FADD"], x["p53"]) - OR(x["C_FLIP"], x["FOXO3A"]),
        "CASP9": lambda: x.__setitem__("CASP9", x["tBid"]), #OR(x["Bax_Bak"], x["tBid"]))),
        "C_FLIP": lambda: x.__setitem__("C_FLIP", x["NFKB"] - x["FOXO3A"]),
        "CREB_1": lambda: x.__setitem__("CREB_1", x["AKT"]),
        "FADD": lambda: x.__setitem__("FADD", x["TNFR"]),
        "FOXO3A": lambda: x.__setitem__("FOXO3A", OR(x["STAT3"], x["MAPK_p38"], x["Nutr_Depr"]) - OR(x["IKKB a/b"], x["AKT"])),
        "HIF_1a": lambda: x.__setitem__("HIF_1a", AND(OR(x["NFKB"], x["mTORC1"]), OR(x["ROS"], x["Hypoxia"]))),
        "Hypoxia": lambda: x.__setitem__("Hypoxia", x["Hypoxia"]),
        "IFN a/b": lambda: x.__setitem__("IFN a/b", x["IRF3"]),
        "IFNR": lambda: x.__setitem__("IFNR", x["IFN a/b"]),
        "IL1": lambda: x.__setitem__("IL1", OR(x["CASP1"], x["CASP8"], x["NFKB"])),
        "IL1R": lambda: x.__setitem__("IL1R", x["IL1"]),
        "IL6": lambda: x.__setitem__("IL6", OR(x["NFKB"], x["MAPK_p38"])),
        "IL6R": lambda: x.__setitem__("IL6R", x["IL6"]),
        "IRF3": lambda: x.__setitem__("IRF3", x["RIG1"]),
        "IKKB a/b": lambda: x.__setitem__("IKKB a/b", OR(x["TLR4"], x["IL1R"], x["TNFR"], x["AKT"])),
        "ISG": lambda: x.__setitem__("ISG", x["STAT1"] - x["Viral_Repl"]),
        "MAPK_p38": lambda: x.__setitem__("MAPK_p38", OR(x["ANG_2_T1R"], x["TLR4"], x["ROS"])),
        # "MLKL": lambda: x.__setitem__("MLKL", x["RIPK1&3"]),
        "mTORC1": lambda: x.__setitem__("mTORC1", x["AKT"] - x["p53"]), #OR(x["IKKB a/b"])),
        "mTORC2": lambda: x.__setitem__("mTORC2", x["PI3K"] - x["mTORC1"]),        
        # "Necroptosis": lambda: x.__setitem__("Necroptosis", x["MLKL"]),
        "NFKB": lambda: x.__setitem__("NFKB", OR(x["IKKB a/b"], x["ROS"]) - x["FOXO3A"]),
        "NLRP3": lambda: x.__setitem__("NLRP3", AND(x["NFKB"], x["RIG1"])),
        "Nutr_Depr": lambda: x.__setitem__("Nutr_Depr", x["Nutr_Depr"]),
        # "Pyroptosis": lambda: x.__setitem__("Pyroptosis", x["CASP1"]),
        "p53": lambda: x.__setitem__("p53", OR(x["Hypoxia"], x["Nutr_Depr"]) - x["Virus"]),
        "PI3K": lambda: x.__setitem__("PI3K", OR(x["TLR4"], x["ROS"], x["IL6R"]) - x['PTEN']),
        "PTEN": lambda: x.__setitem__("PTEN", x["p53"]),
        "RIG1": lambda: x.__setitem__("RIG1", x["Viral_Repl"]),
        # "RIPK1&3": lambda: x.__setitem__("RIPK1&3", OR(x["RIG1"], x["TLR4"], x["FADD"])),
        "ROS": lambda: x.__setitem__("ROS", OR(x["ANG_2_T1R"], x["MAPK_p38"]) - x["FOXO3A"]),
        "SIL6R": lambda: x.__setitem__("SIL6R", OR(x["ADAM_17"], x["IL6"])),
        "STAT1": lambda: x.__setitem__("STAT1", x["IFNR"]), #x["IFNR"]
        "STAT3": lambda: x.__setitem__("STAT3", x["IL6R"]),
        "tBid": lambda: x.__setitem__("tBid", OR(x["CASP8"], x["ROS"])),
        "TLR4": lambda: x.__setitem__("TLR4", x["Virus"]),
        "TNF": lambda: x.__setitem__("TNF", OR(x["ADAM_17"], x["NFKB"], x["MAPK_p38"])),
        "TNFR": lambda: x.__setitem__("TNFR", x["TNF"]),
        # "TSC2": lambda: x.__setitem__("TSC2", (-x["AKT"]) + x["p53"]), # AND(NOT(x["IKKB a/b"]))),
        "Viral_Repl": lambda: x.__setitem__("Viral_Repl", (OR(x["Virus"], x["mTORC1"]) - x["ISG"])),
        "Virus": lambda: x.__setitem__("Virus", x["Virus"])
    }

    temp_mat = np.zeros((n_iter, len(operations), orders))
    modified_operations = {}

    for order in np.arange(orders):
        x = {component: 0 for component in components}
        x["Virus"] = n_states - 1
        x["IFN a/b"] = n_states - 1
        x["ACE2"] = n_states - 1
        x["TSC2"] = n_states - 1
        x["Nutr_Depr"] = 0 # n_states - 1
        x["Hypoxia"] = 0

        for i in np.arange(n_iter):
            indices = np.random.permutation(len(components))
            for idx in indices:
                key = components[idx]
                if key in operations:
                    component, operation_lambda = key, operations[key]
                    modified_operations[key] = lambda x, component=component, operation_lambda=operation_lambda: x.__setitem__(component, g(x, component, operation_lambda, n_states))
                    modified_operations[key](x)
                else:
                    raise NotImplementedError(f"Operation for '{key}' not defined in 'operations'")
            indices = np.sort(indices)
            temp_mat[i, :, order] = [x[components[a]] for a in indices]

    # print components in alphabetical order in the form of [a, b, c, ...]
    # print("[" + ", ".join([f"\"{x}\"" for x in sorted(components, key=lambda x: x.lower())]) + "]")

    return temp_mat

def main():
    np.random.seed(0)
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    plt.rcParams['figure.dpi'] = 300
    n_states = 3

    n_iter = 15
    orders = 500
    # consider simplification in the future, add negative feedback loops

    components = ["ACE2", "ADAM_17", "AKT", "ANG_2", "ANG_2_T1R", 
                  "Apoptosis", "BCL_2", "C_FLIP", "CASP1", "CASP8", 
                  "CASP9", "CREB_1", "FADD", "FOXO3A", "HIF_1a", 
                  "Hypoxia", "IFN a/b", "IFNR", "IKKB a/b", "IL1", 
                  "IL1R", "IL6", "IL6R", "IRF3", "ISG", "MAPK_p38", 
                  "mTORC1", "mTORC2", "NFKB", "NLRP3", "Nutr_Depr", 
                  "p53", "PI3K", "PTEN", "RIG1", "ROS", "SIL6R", 
                  "STAT1", "STAT3", "tBid", "TLR4", "TNF", "TNFR", 
                  "Viral_Repl", "Virus"]
    comp_edit = copy.deepcopy(components)

    if components is None:
        raise NotImplementedError("No components provided")
    if (n_iter < 1) or (orders < 1):
        raise ValueError("Number of iterations and orders must be positive")
    if (n_iter % 1 != 0) or (orders % 1 != 0):
        raise ValueError("Number of iterations and orders must be integers")

    mat = np.average(BN(comp_edit, n_iter, orders, n_states), axis=2)

    yticklabels = [str(x) for x in 1 + np.arange(n_iter)]

    cmap = cm.get_cmap('icefire_r')

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'truncated_%s' % cmap.name,
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    edited_comp = [x.replace("_", " ").replace("a/b", "α/β")
                   for x in components]
    mat_trunc = mat[0:n_iter]
    yticklabels_trunc = yticklabels[0:n_iter]
    ax = sns.heatmap(mat_trunc, cmap=truncate_colormap(cmap, 0.25, 0.75, n=200),
                     linewidths=.05, xticklabels=edited_comp,
                     yticklabels=yticklabels_trunc, vmin=0, vmax=n_states-1, alpha=1)
    ax.tick_params(axis='y', which='major', labelsize=10)

    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks(np.arange(0, n_states, 1))
    colorbar.set_ticklabels([str(int(i)) for i in np.arange(0, n_states, 1)])

    ax.set_xlabel('Model Component', fontsize=12)
    ax.set_ylabel('Iteration Number', fontsize=12)
    ax.set_title(
        f'Model Component Activation in COVID-19 with {orders} Samples',
        fontsize=14)
    
    # zip(components, mat_trunc))

    plt.tight_layout()
    fig.savefig('plot.png')


if __name__ == '__main__':
    main()
    print("Complete")