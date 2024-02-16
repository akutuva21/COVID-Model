import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns

# Define constants
N_STATES = 3
N_ITER = 15
ORDERS = 500

def AVG(*args):
    return np.max(args)
    # return np.mean(args)
def OR(*args):
    return np.max(args)
def AND(*args):
    return np.min(args)
def NOT(value):
    return N_STATES - value

class Pneumocyte:
    def __init__(self, x=None):
        components = ["ACE2", "ADAM_17", "AKT", "ANG_2", "ANG_2_T1R", 
                      "Apoptosis", "BCL_2", "CASP1", "CASP8", "CASP9", "C_FLIP", 
                      "CREB_1", "FADD", "FOXO3A", "HIF_1a", "Hypoxia", "IFN_a_b", 
                      "IFNR", "IL1", "IL1R", "IL6", "IL6R", "IRF3", "IKKB_a_b", 
                      "ISG", "MAPK_p38", "mTORC1", "mTORC2", "NFKB", "NLRP3", 
                      "Nutr_Depr", "p53", "PI3K", "PTEN", "RIG1", "ROS", "SIL6R", 
                      "STAT1", "STAT3", "tBid", "TLR4", "TNF", "TNFR", "Viral_Repl", 
                      "Virus"]
        if x is None:
            self.x = {component: 0 for component in components}
        else:
            self.x = x

    def ACE2(self):
        return self.x["FOXO3A"] - AVG(self.x["Virus"], self.x["ADAM_17"])
    def ADAM_17(self):
        return OR(self.x["ANG_2_T1R"], self.x["HIF_1a"])
    def AKT(self):
        return AVG(self.x["mTORC2"], self.x["PI3K"], self.x["FOXO3A"])
    def ANG_2(self):
        return 0 # not expressed in-vitro
    def ANG_2_T1R(self):
        return self.x["ANG_2"]
    def Apoptosis(self):
        return AVG(self.x["CASP8"], self.x["CASP9"])
    def BCL_2(self):
        return AVG(self.x["NFKB"], self.x["CREB_1"], self.x["HIF_1a"]) - self.x["p53"]
    def CASP1(self):
        return self.x["NLRP3"]
    def CASP8(self):
        return AVG(self.x["FADD"], self.x["p53"]) - AVG(self.x["C_FLIP"], self.x["FOXO3A"])
    def CASP9(self):
        return AVG(self.x["C_FLIP"], self.x["FOXO3A"]) - self.x["BCL_2"]
    def C_FLIP(self):
        return self.x["NFKB"] - self.x["FOXO3A"]
    def CREB_1(self):
        return self.x["AKT"]
    def FADD(self):
        return self.x["TNFR"]
    def FOXO3A(self):
        return AVG(self.x["STAT3"], self.x["MAPK_p38"], self.x["Nutr_Depr"]) - AVG(self.x["IKKB_a_b"], self.x["AKT"])
    def HIF_1a(self):
        return AND(AVG(self.x["NFKB"], self.x["mTORC1"]), 
                   AVG(self.x["ROS"], self.x["Hypoxia"]))
    def Hypoxia(self):
        return self.x["Hypoxia"]
    def IFN_a_b(self):
        return self.x["IRF3"]
    def IFNR(self):
        return self.x["IFN_a_b"]
    def IL1(self):
        return AVG(self.x["CASP1"], self.x["CASP8"], self.x["NFKB"])
    def IL1R(self):
        return self.x["IL1"]
    def IL6(self):
        return AVG(self.x["NFKB"], self.x["MAPK_p38"])
    def IL6R(self):
        return self.x["IL6"]
    def IRF3(self):
        return self.x["RIG1"]
    def IKKB_a_b(self):
        return AVG(self.x["TLR4"], self.x["IL1R"], 
                  self.x["TNFR"], self.x["AKT"])
    def ISG(self):
        return self.x["STAT1"] - self.x["Viral_Repl"]
    def MAPK_p38(self):
        return OR(self.x["ANG_2_T1R"], self.x["TLR4"], self.x["ROS"])
    def mTORC1(self):
        return self.x["AKT"] - self.x["p53"]
    def mTORC2(self):
        return self.x["PI3K"] - self.x["mTORC1"]
    def NFKB(self):
        return AVG(self.x["IKKB_a_b"], self.x["ROS"]) - self.x["FOXO3A"]
    def NLRP3(self):
        return AND(self.x["NFKB"], self.x["RIG1"])
    def Nutr_Depr(self):
        return self.x["Nutr_Depr"]
    def p53(self):
        return AVG(self.x["Hypoxia"], self.x["Nutr_Depr"]) - self.x["Virus"]
    def PI3K(self):
        return AVG(self.x["TLR4"], self.x["ROS"], 
                  self.x["IL6R"]) - self.x['PTEN']
    def PTEN(self):
        return self.x["p53"]
    def RIG1(self):
        return self.x["Viral_Repl"]
    def ROS(self):
        return OR(self.x["ANG_2_T1R"], 
                  self.x["MAPK_p38"]) - self.x["FOXO3A"]
    def SIL6R(self):
        return AVG(self.x["ADAM_17"], self.x["IL6"])
    def STAT1(self):
        return self.x["IFNR"]
    def STAT3(self):
        return self.x["IL6R"]
    def tBid(self):
        return AVG(self.x["CASP8"], self.x["ROS"])
    def TLR4(self):
        return self.x["Virus"]
    def TNF(self):
        return AVG(self.x["ADAM_17"], self.x["NFKB"], self.x["MAPK_p38"])
    def TNFR(self):
        return self.x["TNF"]
    def Viral_Repl(self):
        return AVG(self.x["Virus"], self.x["mTORC1"]) - self.x["ISG"]
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
    '''
    components: list of components in the model
    n_iter: number of iterations to run the model
    orders: number of orders to run the model
    n_states: number of states for each component
    '''

    cell = Pneumocyte()
    operations = {component: cell.get_oper(component) for component in components}
    indices = np.arange(len(operations))
    temp_mat = np.zeros((n_iter, len(operations), orders))

    for order in np.arange(orders):
        cell.x = {component: 0 for component in components}
        cell.x["Virus"] = n_states - 1
        cell.x["IFN_a_b"] = n_states - 1
        cell.x["ACE2"] = n_states - 1
        cell.x["TSC2"] = n_states - 1
        cell.x["Nutr_Depr"] = 0 # n_states - 1
        cell.x["Hypoxia"] = 0

        for i in np.arange(n_iter):
            np.random.shuffle(indices)
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
    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    plt.rcParams['figure.dpi'] = 300

    components = list(Pneumocyte().x.keys())
    mat = np.average(BN(components, N_ITER, ORDERS, N_STATES), axis=2)

    cmap = cm.get_cmap('icefire_r')

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'truncated_%s' % cmap.name,
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    yticklabels = [str(x) for x in 1 + np.arange(N_ITER)]
    components = [x.replace("_a_b", " α/β").replace("_", " ") for x in components]

    ax = sns.heatmap(mat, cmap=truncate_colormap(cmap, 0.25, 0.75, n=200),
                     linewidths=.05, xticklabels=components,
                     yticklabels=yticklabels, vmin=0, vmax=N_STATES-1, alpha=1)
    ax.tick_params(axis='y', which='major', labelsize=10)

    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks(np.arange(0, N_STATES, 1))
    colorbar.set_ticklabels([str(int(i)) for i in np.arange(0, N_STATES, 1)])

    ax.set_xlabel('Model Component', fontsize=12)
    ax.set_ylabel('Iteration Number', fontsize=12)
    ax.set_title(
        f'Model Component Activation in COVID-19 with {ORDERS} Samples',
        fontsize=14)
    
    plt.tight_layout()
    fig.savefig('plot.png')

if __name__ == '__main__':
    main()
    print("Complete")