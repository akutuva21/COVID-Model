import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns
import copy

def BN(components, n_iter=25, orders=500):
    '''Generates a boolean network based on the provided components

    Parameters
    ----------
    components : list
    n_iter : int, optional
        Number of iterations to run the network for (default is 25)
    orders : int, optional
        Number of samples to take for each iteration (default is 250)
        Averaged at the end of the function
    
    Raises
    ----------
        NotImplementedError
            If the component is not implemented
            If an 'if' statement was not executed (check names)
        ValueError
            If the number of orders or iterations is less than 1
            If the number of orders or iterations is not an integer
    '''

    if components is None:
        raise NotImplementedError("No components provided")
    if (n_iter < 1) or (orders < 1):
        raise ValueError("Number of iterations and orders must be positive")
    if (n_iter % 1 != 0) or (orders % 1 != 0):
        raise ValueError("Number of iterations and orders must be integers")
    
    temp_mat = np.zeros((n_iter, len(components), orders))
    operations = {
        "ACE2": lambda: x.__setitem__("ACE2", not x["Virus"] or x["FOXO3A"]),
        "ADAM_17": lambda: x.__setitem__("ADAM_17", x["ANG_2_T1R"]),
        "AKT": lambda: x.__setitem__("AKT", (x["mTORC2"] or x["PI3K"]) and not x["FOXO3A"]),
        "ANG_2": lambda: x.__setitem__("ANG_2", not x["ACE2"]),
        "ANG_2_T1R": lambda: x.__setitem__("ANG_2_T1R", x["ANG_2"]),
        "Apoptosis": lambda: x.__setitem__("Apoptosis", x["CASP8"] or x["CASP9"]),
        "Autophagy": lambda: x.__setitem__("Autophagy", not x["mTORC1"]),
        "BCL_2": lambda: x.__setitem__("BCL_2", (x["NFKB"] or x["CREB_1"] or x["HIF_1a"])),
        "Bax_Bak": lambda: x.__setitem__("Bax_Bak", not x["BCL_2"]),
        "CASP1": lambda: x.__setitem__("CASP1", x["NLRP3"]),
        "CASP8": lambda: x.__setitem__("CASP8", (x["FADD"] and not x["C_FLIP"]) and not x["FOXO3A"]),
        "CASP9": lambda: x.__setitem__("CASP9", x["Bax_Bak"] or x["tBid"]),
        "C_FLIP": lambda: x.__setitem__("C_FLIP", x["NFKB"] and not x["FOXO3A"]),
        "CREB_1": lambda: x.__setitem__("CREB_1", x["AKT"]),
        "FADD": lambda: x.__setitem__("FADD", x["TNFR"]),
        "FOXO3A": lambda: x.__setitem__("FOXO3A", (x["STAT3"] or x["MAPK_p38"]) and not (x["IKKB a/b"] or x["AKT"])),
        "HIF_1a": lambda: x.__setitem__("HIF_1a", (x["NFKB"] or x["mTORC1"]) and x["ROS"]),
        "IFN a/b": lambda: x.__setitem__("IFN a/b", (x["IRF3"]) and not (x["Viral_Repl"])),
        "IFNR": lambda: x.__setitem__("IFNR", x["STAT1"]),
        "IL1": lambda: x.__setitem__("IL1", x["CASP1"] or x["CASP8"] or x["NFKB"]),
        "IL1R": lambda: x.__setitem__("IL1R", x["IL1"]),
        "IL6": lambda: x.__setitem__("IL6", x["NFKB"] or x["MAPK_p38"]),
        "IL6R": lambda: x.__setitem__("IL6R", x["IL6"]),
        "IRF3": lambda: x.__setitem__("IRF3", x["RIG1"]),
        "IKKB a/b": lambda: x.__setitem__("IKKB a/b", (x["TLR4"] or x["IL1R"] or x["TNFR"] or x["AKT"])),
        "ISG": lambda: x.__setitem__("ISG", x["STAT1"]),
        "MAPK_p38": lambda: x.__setitem__("MAPK_p38", (x["ANG_2_T1R"] or x["TLR4"] or x["ROS"])),
        "MLKL": lambda: x.__setitem__("MLKL", x["RIPK1&3"]),
        "mTORC1": lambda: x.__setitem__("mTORC1", not x["TSC2"]),
        "mTORC2": lambda: x.__setitem__("mTORC2", x["PI3K"]),
        "Necroptosis": lambda: x.__setitem__("Necroptosis", x["MLKL"]),
        "NFKB": lambda: x.__setitem__("NFKB", (x["IKKB a/b"] or x["ROS"]) and not x["FOXO3A"]),
        "NLRP3": lambda: x.__setitem__("NLRP3", x["NFKB"] and x["RIG1"]),
        "Pyroptosis": lambda: x.__setitem__("Pyroptosis", x["CASP1"]),
        "PI3K": lambda: x.__setitem__("PI3K", x["TLR4"] or x["ROS"] or x["IL6R"]),
        "RIG1": lambda: x.__setitem__("RIG1", x["Viral_Repl"]),
        "RIPK1&3": lambda: x.__setitem__("RIPK1&3", (x["RIG1"] or x["TLR4"] or x["FADD"])),
        "ROS": lambda: x.__setitem__("ROS", (x["ANG_2_T1R"] or x["MAPK_p38"]) and not x["FOXO3A"]),
        "SIL6R": lambda: x.__setitem__("SIL6R", x["ADAM_17"] or x["IL6"]),
        "STAT1": lambda: x.__setitem__("STAT1", x["IFNR"]),
        "STAT3": lambda: x.__setitem__("STAT3", x["IL6R"]),
        "tBid": lambda: x.__setitem__("tBid", x["CASP8"] or x["ROS"]),
        "TLR4": lambda: x.__setitem__("TLR4", x["Virus"]),
        "TNF": lambda: x.__setitem__("TNF", x["ADAM_17"] or x["NFKB"] or x["MAPK_p38"]),
        "TNFR": lambda: x.__setitem__("TNFR", x["TNF"]),
        "TSC2": lambda: x.__setitem__("TSC2", not x["AKT"]),
        "Viral_Repl": lambda: x.__setitem__("Viral_Repl", (x["Virus"] and x["mTORC2"]) and not x["ISG"]),
        "Virus": lambda: x.__setitem__("Virus", x["Virus"])
    }

    x = {component: False for component in components}
    x["Virus"] = True
    x["ACE2"] = True
    
    for order in np.arange(orders):
        for i in np.arange(n_iter):
            indices = np.random.permutation(len(components))
            for idx in indices:
                key = components[idx]
                if key in operations:
                    operations[key]()
                else:
                    raise NotImplementedError(f"Operation for '{key}' not defined in 'operations'")
            indices = np.sort(indices)
            temp_mat[i, :, order] = [x[components[a]] for a in indices]

    return temp_mat

def main():
    np.random.seed(0)
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    plt.rcParams['figure.dpi'] = 300

    components = ["Virus", "Viral_Repl", "ACE2", "ANG_2", 
                  "ANG_2_T1R", "ADAM_17", "TLR4", "RIG1", "NFKB", "SIL6R",
                  "IKKB a/b", "TNF", "IRF3", "STAT1", "STAT3", "IL6", "IL6R",
                  "ISG", "C_FLIP", "IFN a/b", "NLRP3", "CASP1", "FOXO3A",
                  "IFNR", "BCL_2", "tBid", "Bax_Bak", "CASP9", "ROS", "TNFR",
                  "FADD", "Pyroptosis", "IL1", "IL1R", "MLKL", "Necroptosis",
                  "RIPK1&3", "CASP8", "Apoptosis", "PI3K", "AKT",
                  "TSC2", "mTORC1", "mTORC2", "CREB_1",
                  "Autophagy", "MAPK_p38", "HIF_1a"]
    components.sort()
    comp_edit = copy.deepcopy(components)
    n_iter = 25
    orders = 250
    # mat = np.zeros((n_iter + 1, len(components), orders))
    # consider simplification in the future, add negative feedback loops

    mat = np.average(BN(comp_edit, n_iter, orders), axis=2)

    yticklabels = [str(x) for x in 1 + np.arange(n_iter)]
    # # yticklabels.append("Average")

    cmap = cm.get_cmap('icefire_r')

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'truncated_%s' % cmap.name,
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    edited_comp = [x.replace("_", " ").replace("a/b", "α/β")
                   for x in components]
    mat_trunc = mat[0:10]
    yticklabels_trunc = yticklabels[0:10]
    ax = sns.heatmap(mat_trunc, cmap=truncate_colormap(cmap, 0.25, 0.75, n=200),
                     linewidths=.05, xticklabels=edited_comp,
                     yticklabels=yticklabels_trunc, vmin=0, vmax=1, alpha=1)
    ax.tick_params(axis='y', which='major', labelsize=10)

    # colorbar = ax.collections[0].colorbar
    # xmin, xmax, delta = np.min(mat), np.max(mat), 0.1
    # colorbar.set_ticks(np.arange(xmin, xmax + delta, delta))

    ax.set_xlabel('Model Component', fontsize=12)
    ax.set_ylabel('Iteration Number', fontsize=12)
    ax.set_title(
        f'Model Component Activation in COVID-19 with {orders} Samples',
        fontsize=14)

    plt.tight_layout()
    fig.savefig('plot.png')


if __name__ == '__main__':
    main()
    print("Complete")