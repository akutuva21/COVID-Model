import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns
import time

# Define constants
N_STATES = 2
MAX_STATES = N_STATES - 1
N_ITER = 45
ORDERS = 500

def AVG(*args):
    return np.mean(args)
def OR(*args):
    return np.mean(args)
def AND(*args):
    return np.min(args)
def NOT(value):
    return MAX_STATES - value

class Pneumocyte:
    MODEL = 0
    SCENARIO  = 0
    def __init__(self, x=None):

        components = ["ACE2", "ADAM_17", "AKT", "AMPK", "ANG_2", "ANG_2_T1R",
                      "Apoptosis", "BCL_2", "CASP1", "CASP8", "CASP9", "C_FLIP",
                      "CREB_1", "FADD", "FOXO3A", "HIF_1a", "Hypoxia", "IFN_a_b",
                      "IFNR", "IL1", "IL1R", "IL6", "IL6R", "IRF3", "IKKB_a_b",
                      "ISG", "MAPK_p38", "mTORC1", "mTORC2", "NFKB", "NLRP3",
                      "Nutr_Depr", "p53", "PI3K", "PTEN", "RIG1", "ROS", "SIL6R",
                      "STAT1", "STAT3", "TLR4", "TNF", "TNFR", "TSC2", "Viral_Repl", 
                      "Virus"]
        if x is None:
            self.x = {component: 0 for component in components}
            self.x["ACE2"] = N_STATES - 1
            self.x["Virus"] = N_STATES - 1
        else:
            self.x = x

    def ACE2(self):
        return self.x["FOXO3A"] - AVG(self.x["Virus"], self.x["ADAM_17"])
    def ADAM_17(self):
        return AVG(self.x["ANG_2_T1R"], self.x["HIF_1a"])
    def AKT(self):
        return AVG(self.x["mTORC2"], self.x["FOXO3A"], self.x["PI3K"]) - self.x["PTEN"]
    def AMPK(self):
        return self.x["Nutr_Depr"] - self.x["AKT"]
    def ANG_2(self):
        return NOT(self.x["ACE2"]) # not expressed in-vitro
    def ANG_2_T1R(self):
        return self.x["ANG_2"]
    def Apoptosis(self):
        return AVG(self.x["CASP8"], self.x["CASP9"])
    def BCL_2(self):
        return AVG(self.x["NFKB"], self.x["CREB_1"]) - self.x["p53"]
    def CASP1(self):
        return self.x["NLRP3"]
    def CASP8(self):
        return AVG(self.x["FADD"], self.x["p53"]) - AVG(self.x["C_FLIP"], self.x["FOXO3A"])
    def CASP9(self):
        return self.x["FOXO3A"] - self.x["BCL_2"]
    def C_FLIP(self):
        return self.x["NFKB"] - self.x["FOXO3A"]
    def CREB_1(self):
        return self.x["AKT"]
    def FADD(self):
        return self.x["TNFR"]
    def FOXO3A(self):
        return AVG(self.x["STAT3"], self.x["MAPK_p38"], self.x["AMPK"]) - AVG(self.x["IKKB_a_b"], self.x["AKT"])
    def HIF_1a(self):
        return AVG(self.x["ROS"], self.x["Hypoxia"])
    def Hypoxia(self):
        return self.x["Hypoxia"]
    def IFN_a_b(self):
        if Pneumocyte.MODEL > 0:
            if Pneumocyte.SCENARIO > 0:
                return AVG(self.x["IRF3"], self.x["mTORC1"]) - self.x["Viral_Repl"]
            else:
                return AVG(self.x["IRF3"], self.x["mTORC1"])
        else:
            if Pneumocyte.SCENARIO > 0:
                return self.x["IRF3"] - self.x["Viral_Repl"]
            else:
                return self.x["IRF3"]
    def IFNR(self):
        return self.x["IFN_a_b"]
    def IL1(self):
        if Pneumocyte.MODEL == MAX_STATES:
            return AVG(self.x["CASP1"], self.x["CASP8"], self.x["NFKB"], self.x["mTORC1"]) 
        else:
            return AVG(self.x["CASP1"], self.x["CASP8"], self.x["NFKB"])
    def IL1R(self):
        return self.x["IL1"]
    def IL6(self):
        if Pneumocyte.MODEL == MAX_STATES:
            return AVG(self.x["NFKB"], self.x["MAPK_p38"], self.x["mTORC1"]) 
        else:
            return AVG(self.x["NFKB"], self.x["MAPK_p38"])
    def IL6R(self):
        return self.x["IL6"]
    def IRF3(self):
        return self.x["RIG1"]
    def IKKB_a_b(self):
        return AVG(self.x["TLR4"], self.x["IL1R"],
                   self.x["TNFR"], self.x["AKT"])
    def ISG(self):
        if Pneumocyte.SCENARIO == MAX_STATES:
            return self.x["STAT1"] - self.x["Viral_Repl"] 
        else:
            return self.x["STAT1"]
    def MAPK_p38(self):
        return AVG(self.x["ANG_2_T1R"], self.x["TLR4"], self.x["ROS"])
    def mTORC1(self):
        return AVG(NOT(self.x["TSC2"]), NOT(self.x["FOXO3A"]))
    def mTORC2(self):
        return AVG(self.x["PI3K"], self.x["AMPK"])
    def NFKB(self):
        return AVG(self.x["IKKB_a_b"], self.x["ROS"]) - self.x["FOXO3A"]
    def NLRP3(self):
        return AVG(self.x["NFKB"], self.x["RIG1"])
    def Nutr_Depr(self):
        return self.x["Nutr_Depr"]
    def p53(self):
        return AVG(self.x["Hypoxia"], self.x["AMPK"]) - self.x["Viral_Repl"]
    def PI3K(self):
        return AVG(self.x["TLR4"], self.x["ROS"],
                   self.x["IL6R"])
    def PTEN(self):
        return self.x["p53"]
    def RIG1(self):
        return self.x["Viral_Repl"]
    def ROS(self):
        return AVG(self.x["ANG_2_T1R"],
                  self.x["MAPK_p38"]) - self.x["FOXO3A"]
    def SIL6R(self):
        return AVG(self.x["ADAM_17"], self.x["IL6R"])
    def STAT1(self):
        return self.x["IFNR"]
    def STAT3(self):
        return self.x["IL6R"]
    def TLR4(self):
        return self.x["Virus"]
    def TNF(self):
        if Pneumocyte.MODEL == MAX_STATES:
            return AVG(self.x["ADAM_17"], self.x["NFKB"], self.x["MAPK_p38"], self.x["mTORC1"]) 
        else:
            return AVG(self.x["ADAM_17"], self.x["NFKB"], self.x["MAPK_p38"])
    def TNFR(self):
        return self.x["TNF"]
    def TSC2(self):
        return AVG(self.x["p53"], self.x["AMPK"], self.x["HIF_1a"]) - self.x["AKT"]
    def Viral_Repl(self):
        return AVG(self.x["Virus"], self.x["mTORC1"]) - self.x["ISG"]
    def Virus(self):
        return self.x["Virus"]

    def set_oper(self, component):
        return lambda: self.x.__setitem__(component, getattr(self, component)())

    def update_state(self, x, index, update_func):
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
        elif current_state > MAX_STATES:
            current_state = MAX_STATES

        return current_state

def QN(components, n_iter=N_ITER, orders=ORDERS, n_states=N_STATES, model=0, scenario=0, nutr_depr=0, hypoxia=0):
    '''
    components: list of components in the model
    n_iter: number of iterations to run the model
    orders: number of orders to run the model
    n_states: number of states for each component
    '''

    cell = Pneumocyte()
    operations = {component: cell.set_oper(component) for component in components}
    
    for component in components:
        if component not in operations:
            raise NotImplementedError(f"Operation for '{component}' not defined in 'operations'")

    temp_mat = np.zeros((n_iter, len(components), orders))

    for order in range(orders):
        cell.x = {component: 0 for component in components}
        cell.x["Virus"] = MAX_STATES
        cell.x["ACE2"] = MAX_STATES
        cell.x["IFN_a_b"] = 0 # MAX_STATES
        cell.x["Nutr_Depr"] = nutr_depr # MAX_STATES
        cell.x["Hypoxia"] = hypoxia # MAX_STATES
        Pneumocyte.MODEL = model
        Pneumocyte.SCENARIO = scenario

        indices = np.arange(len(components))

        for i in range(n_iter):
            np.random.shuffle(indices)
            for idx in indices:
                key = components[idx]
                component, oper_func = key, operations[key]
                cell.x[component] = cell.update_state(cell.x, component, oper_func)
            
            temp_mat[i, :, order] = [cell.x[comp] for comp in components]

    return temp_mat

def main():
    nutr_cond = ["base", "starvation"]
    hyp_cond = ["normoxia", "hypoxia"]
    use_nutr_cond = False
    use_hyp_cond = False

    # base_name = ["mtorc"]
    SIZE = N_ITER
    for nutr_depr in range(len(nutr_cond)) if use_nutr_cond else [0]:
        nname = nutr_cond[nutr_depr] if use_nutr_cond else nutr_cond[0]
        
        for hypoxia in range(len(hyp_cond)) if use_hyp_cond else [0]:
            hname = hyp_cond[hypoxia] if use_hyp_cond else hyp_cond[0]

            for model in [0, MAX_STATES]:
                for scenario in [0, MAX_STATES]:
                    print("# " + nname + " + " + hname + " model " + str(model) + " scenario " + str(scenario))
                    name = nname + "." + hname + "." + str(model) + "." + str(scenario) + "."
                    matt = np.zeros((SIZE,3,5))
                    comp = ["APO", "IFN", "IL6", "TNF", "Vir.rep"]
                    comp_idx = [6, 17, 21, 42, 45]
                    rows = [15, 30, 45]

                    for ii in range(SIZE):
                        # np.random.seed(0)
                        # fig, ax = plt.subplots(1, 1, figsize=(9, 7))
                        # plt.rcParams['figure.dpi'] = 300

                        components = list(Pneumocyte().x.keys())
                        mat = np.average(QN(components, N_ITER, ORDERS, N_STATES, model=model, scenario=scenario, nutr_depr=nutr_depr), axis=2)
                        for jj, row in enumerate([r - 1 for r in rows]):
                            matt[ii, jj, :] = mat[row, comp_idx]
                        # cmap = cm.get_cmap('icefire_r')

                        # def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
                        #     new_cmap = colors.LinearSegmentedColormap.from_list(
                        #         'truncated_%s' % cmap.name,
                        #         cmap(np.linspace(minval, maxval, n)))
                        #     return new_cmap

                        # yticklabels = [str(x) for x in 1 + np.arange(N_ITER)]
                        # components = [x.replace("_a_b", " α/β").replace("_", " ") for x in components]

                        # ax = sns.heatmap(mat, cmap=truncate_colormap(cmap, 0.25, 0.75, n=200),
                        #                 linewidths=.05, xticklabels=components,
                        #                 yticklabels=yticklabels, vmin=0, vmax=N_STATES-1, alpha=1)
                        # ax.tick_params(axis='y', which='major', labelsize=10)

                        # colorbar = ax.collections[0].colorbar
                        # colorbar.set_ticks(np.arange(0, N_STATES, 1))
                        # colorbar.set_ticklabels([str(int(i)) for i in np.arange(0, N_STATES, 1)])

                        # ax.set_xlabel('Model Component', fontsize=12)
                        # ax.set_ylabel('Iteration Number', fontsize=12)
                        # ax.set_title(
                        #    f'Model Component Activation in COVID-19 with {ORDERS} Samples',
                        #    fontsize=14)

                        # plt.tight_layout()
                        # fig.savefig('QN_plot.png')

                    # print each value of components and the corresponding value in the last row of the matrix
                    # for i in range(len(components)):
                    #    print(components[i], mat[-1, i])
                    for i in range(len(comp)):
                        for k, kk in enumerate(rows):  # Loop through the three rows (14, 29, 44)
                            print(f"{name}{comp[i]}_{kk} = c(", end="")
                            for j in range(SIZE):
                                if j == 0:
                                    print(matt[j, k, i], end="")
                                else:
                                    print(f", {matt[j, k, i]}", end="")
                            print(")")
                    print("\n")

    #### why 3?
    for j in range(3):
        for x in comp:
            print("boxplot(cbind(base." + str(0) + "." + str(j) + "." + x +
                  ", base." + str(1) + "." + str(j) + "." + x +
                  ", base." + str(2) + "."  + str(j) + "."  + x +
                  ", starvation." + str(0) + "."  + str(j) + "." + x +
                  ", starvation."  + str(1) + "." + str(j) + "." + x +
                  ", starvation." + str(2) + "." + str(j) + "." + x + ")" +
                  ", ylab=" + x + ", main=secenario " + str(j) + ", las=2)"
                )
    print("\n")

def run():

    # np.random.seed(0)
    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    plt.rcParams['figure.dpi'] = 300
    print("c(", end="")

    vec=list(Pneumocyte().x.values())
    for i in range(len(vec)):
        if i == len(vec) - 1:
            print(str(vec[i]), end="")
        else:
            print(str(vec[i]) + ", ", end="")
    print("),")

    components = list(Pneumocyte().x.keys())
    mat = np.average(QN(components, N_ITER, ORDERS, N_STATES, model=2, scenario=2, nutr_depr=0), axis=2)
    for vec in mat:
        print("c(", end="")
        for i in range(len(vec)):
            if i == len(vec) - 1:
                print(str(vec[i]), end="")
            else:
                print(str(vec[i]) + ", ", end="")
        print("),")

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
    fig.savefig('QN_plot.png')

    # print each value of components and the corresponding value in the last row of the matrix

if __name__ == '__main__':
    start_time = time.time()

    main()

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time} seconds")
    print("Complete")