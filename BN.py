import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns
import copy

# Define constants
N_STATES = 3
N_ITER = 25
ORDERS = 500

def OR(*args):
    return np.mean(args)
def AND(*args):
    return np.min(args)
def NOT(value):
    return N_STATES - value

class Pneumocyte:
    def __init__(self, x):
        self.x = x
    
    def ACE2(self):
        return self.x["FOXO3A"] - OR(self.x["Virus"], self.x["ADAM_17"])
    def ADAM_17(self):
        return OR(self.x["ANG_2_T1R"], self.x["HIF_1a"])
    def AKT(self):
        return OR(self.x["mTORC2"], self.x["PI3K"], self.x["FOXO3A"])
    def ANG_2(self):
        return NOT(self.x["ACE2"])
    def ANG_2_T1R(self):
        return self.x["ANG_2"]
    def Apoptosis(self):
        return OR(self.x["CASP8"], self.x["CASP9"])
    def BCL_2(self):
        return OR(self.x["NFKB"], self.x["CREB_1"], self.x["HIF_1a"]) - self.x["p53"]
    def CASP1(self):
        return self.x["NLRP3"]
    def CASP8(self):
        return OR(self.x["FADD"], self.x["p53"]) - OR(self.x["C_FLIP"], self.x["FOXO3A"])
    def CASP9(self):
        return OR(self.x["C_FLIP"], self.x["FOXO3A"]) - self.x["BCL_2"]
    def C_FLIP(self):
        return self.x["NFKB"] - self.x["FOXO3A"]
    def CREB_1(self):
        return self.x["AKT"]
    def FADD(self):
        return self.x["TNFR"]
    def FOXO3A(self):
        return OR(self.x["STAT3"], self.x["MAPK_p38"], self.x["Nutr_Depr"]) - OR(self.x["IKKB_a_b"], self.x["AKT"])
    def HIF_1a(self):
        return AND(OR(self.x["NFKB"], self.x["mTORC1"]), OR(self.x["ROS"], self.x["Hypoxia"]))
    def Hypoxia(self):
        return self.x["Hypoxia"]
    def IFN_a_b(self):
        return self.x["IRF3"]
    def IFNR(self):
        return self.x["IFN_a_b"]
    def IL1(self):
        return OR(self.x["CASP1"], self.x["CASP8"], self.x["NFKB"])
    def IL1R(self):
        return self.x["IL1"]
    def IL6(self):
        return OR(self.x["NFKB"], self.x["MAPK_p38"])
    def IL6R(self):
        return self.x["IL6"]
    def IRF3(self):
        return self.x["RIG1"]
    def IKKB_a_b(self):
        return OR(self.x["TLR4"], self.x["IL1R"], self.x["TNFR"], self.x["AKT"])
    def ISG(self):
        return self.x["STAT1"] - self.x["Viral_Repl"]
    def MAPK_p38(self):
        return OR(self.x["ANG_2_T1R"], self.x["TLR4"], self.x["ROS"])
    def mTORC1(self):
        return self.x["AKT"] - self.x["p53"]
    def mTORC2(self):
        return self.x["PI3K"] - self.x["mTORC1"]
    def NFKB(self):
        return OR(self.x["IKKB_a_b"], self.x["ROS"]) - self.x["FOXO3A"]
    def NLRP3(self):
        return AND(self.x["NFKB"], self.x["RIG1"])
    def Nutr_Depr(self):
        return self.x["Nutr_Depr"]
    def p53(self):
        return OR(self.x["Hypoxia"], self.x["Nutr_Depr"]) - self.x["Virus"]
    def PI3K(self):
        return OR(self.x["TLR4"], self.x["ROS"], self.x["IL6R"]) - self.x['PTEN']
    def PTEN(self):
        return self.x["p53"]
    def RIG1(self):
        return self.x["Viral_Repl"]
    def ROS(self):
        return OR(self.x["ANG_2_T1R"], self.x["MAPK_p38"]) - self.x["FOXO3A"]
    def SIL6R(self):
        return OR(self.x["ADAM_17"], self.x["IL6"])
    def STAT1(self):
        return self.x["IFNR"]
    def STAT3(self):
        return self.x["IL6R"]
    def tBid(self):
        return OR(self.x["CASP8"], self.x["ROS"])
    def TLR4(self):
        return self.x["Virus"]
    def TNF(self):
        return OR(self.x["ADAM_17"], self.x["NFKB"], self.x["MAPK_p38"])
    def TNFR(self):
        return self.x["TNF"]
    def Viral_Repl(self):
        return (OR(self.x["Virus"], self.x["mTORC1"]) - self.x["ISG"])
    def Virus(self):
        return self.x["Virus"]
    
    def get_oper(self, component):
       return lambda: self.x.__setitem__(component, getattr(self, component)())

    def update_state(self, x, index, update_func, n_states=N_STATES):
        self.x = x
        current_state = self.x[index]
        update_func()
        new_state = self.x[index]

        if new_state > current_state:
            current_state += 1
        elif new_state < current_state:
            current_state -= 1

        if current_state < 0:
            current_state = 0
        elif current_state > (n_states - 1):
            current_state = n_states - 1

        return current_state

def BN(components, n_iter=N_ITER, orders=ORDERS, n_states=N_STATES):
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

    cell = Pneumocyte({component: 0 for component in components})
    operations = {component: cell.get_oper(component) for component in components}
    
    temp_mat = np.zeros((n_iter, len(operations), orders))
    modified_operations = {}

    for order in np.arange(orders):
        cell.x = {component: 0 for component in components}
        cell.x["Virus"] = n_states - 1
        cell.x["IFN_a_b"] = n_states - 1
        cell.x["ACE2"] = n_states - 1
        cell.x["TSC2"] = n_states - 1
        cell.x["Nutr_Depr"] = 0 # n_states - 1
        cell.x["Hypoxia"] = 0

        for i in np.arange(n_iter):
            indices = np.random.permutation(len(operations))
            for idx in indices:
                key = components[idx]
                if key in operations:
                    component, operation_lambda = key, operations[key]
                    cell.x[component] = cell.update_state(cell.x, component, operation_lambda, n_states)
                else:
                    raise NotImplementedError(f"Operation for '{key}' not defined in 'operations'")
            indices = np.sort(indices)
            temp_mat[i, :, order] = [cell.x[components[a]] for a in indices]

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
                  "Hypoxia", "IFN_a_b", "IFNR", "IKKB_a_b", "IL1", 
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

    mat_trunc = mat[0:n_iter]
    yticklabels_trunc = yticklabels[0:n_iter]
    ax = sns.heatmap(mat_trunc, cmap=truncate_colormap(cmap, 0.25, 0.75, n=200),
                     linewidths=.05, xticklabels=components,
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
    
    plt.tight_layout()
    fig.savefig('plot.png')


if __name__ == '__main__':
    main()
    print("Complete")