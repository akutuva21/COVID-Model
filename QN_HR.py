import numpy as np

N_STATES = 2 # number of states for each component
MAX_STATES = N_STATES - 1 # maximum state value
N_ITER = 50 # number of iterations
ORDERS = 500 # number of models averaged
SIZE = 45 # number of samples in each boxplot

def AVG(*args):
    return np.mean(args)
def OR(*args):
    return np.mean(args)
def AND(*args):
    return np.min(args)
def NOT(value):
    return MAX_STATES - value

class Pneumocyte:
    MTORC_Scenario = -1
    SCENARIO = -1
    MTORC_knockout = False
    IFN_knockout = False

    components = ["ACE2", "ADAM_17", "AKT", "AMPK", "ANG_2", "ANG_2_T1R",
                    "Apoptosis", "BCL_2", "CASP1", "CASP8", "CASP9", "C_FLIP",
                    "CREB_1", "FADD", "FOXO3A", "HIF_1a", "Hypoxia", "IFN_a_b",
                    "IFNR", "IL1", "IL1R", "IL6", "IL6R", "IRF3", "IKKB_a_b",
                    "ISG", "MAPK_p38", "mTORC1", "mTORC2", "NFKB", "NLRP3",
                    "Nutr_Depr", "p53", "PI3K", "PTEN", "RIG1", "ROS", "SIL6R",
                    "STAT1", "STAT3", "TLR4", "TNF", "TNFR", "TSC2", "Viral_Repl", 
                    "Virus"]
    comp_idx = {comp: idx for idx, comp in enumerate(components)}
    
    def __init__(self, x=None):
        self.x = np.zeros(len(self.components)) if x is None else x
        self.x[self.comp_idx["Virus"]] = MAX_STATES
        self.x[self.comp_idx["ACE2"]] = MAX_STATES

    def compute_ACE2(self):
        return self.x[self.comp_idx["FOXO3A"]] - AVG(self.x[self.comp_idx["Virus"]], self.x[self.comp_idx["ADAM_17"]])
    def compute_ADAM_17(self):
        return AVG(self.x[self.comp_idx["ANG_2_T1R"]], self.x[self.comp_idx["HIF_1a"]])
    def compute_AKT(self):
        return AVG(self.x[self.comp_idx["mTORC2"]], self.x[self.comp_idx["FOXO3A"]], self.x[self.comp_idx["PI3K"]]) - self.x[self.comp_idx["PTEN"]]
    def compute_AMPK(self):
        return self.x[self.comp_idx["Nutr_Depr"]] - self.x[self.comp_idx["AKT"]]
    def compute_ANG_2(self):
        return NOT(self.x[self.comp_idx["ACE2"]]) # not expressed in-vitro
    def compute_ANG_2_T1R(self):
        return self.x[self.comp_idx["ANG_2"]]
    def compute_Apoptosis(self):
        return AVG(self.x[self.comp_idx["CASP8"]], self.x[self.comp_idx["CASP9"]])
    def compute_BCL_2(self):
        return AVG(self.x[self.comp_idx["NFKB"]], self.x[self.comp_idx["CREB_1"]]) - self.x[self.comp_idx["p53"]]
    def compute_CASP1(self):
        return self.x[self.comp_idx["NLRP3"]]
    def compute_CASP8(self):
        return AVG(self.x[self.comp_idx["FADD"]], self.x[self.comp_idx["p53"]]) - AVG(self.x[self.comp_idx["C_FLIP"]], self.x[self.comp_idx["FOXO3A"]])
    def compute_CASP9(self):
        return self.x[self.comp_idx["FOXO3A"]] - self.x[self.comp_idx["BCL_2"]]
    def compute_C_FLIP(self):
        return self.x[self.comp_idx["NFKB"]] - self.x[self.comp_idx["FOXO3A"]]
    def compute_CREB_1(self):
        return self.x[self.comp_idx["AKT"]]
    def compute_FADD(self):
        return self.x[self.comp_idx["TNFR"]]
    def compute_FOXO3A(self):
        return AVG(self.x[self.comp_idx["STAT3"]], self.x[self.comp_idx["MAPK_p38"]], self.x[self.comp_idx["AMPK"]]) - AVG(self.x[self.comp_idx["IKKB_a_b"]], self.x[self.comp_idx["AKT"]])
    def compute_HIF_1a(self):
        return AVG(self.x[self.comp_idx["ROS"]], self.x[self.comp_idx["Hypoxia"]])
    def compute_Hypoxia(self):
        return self.x[self.comp_idx["Hypoxia"]]
    def compute_IFN_a_b(self):
        if Pneumocyte.IFN_knockout == True:
            return 0
        if Pneumocyte.MTORC_Scenario > 0:
            if Pneumocyte.SCENARIO > 0:
                return AVG(self.x[self.comp_idx["IRF3"]], self.x[self.comp_idx["mTORC1"]]) - self.x[self.comp_idx["Viral_Repl"]]
            else:
                return AVG(self.x[self.comp_idx["IRF3"]], self.x[self.comp_idx["mTORC1"]])
        else:
            if Pneumocyte.SCENARIO > 0:
                return self.x[self.comp_idx["IRF3"]] - self.x[self.comp_idx["Viral_Repl"]]
            else:
                return self.x[self.comp_idx["IRF3"]]
    def compute_IFNR(self):
        return self.x[self.comp_idx["IFN_a_b"]]
    def compute_IL1(self):
        if Pneumocyte.MTORC_Scenario == 2:
            return AVG(self.x[self.comp_idx["CASP1"]], self.x[self.comp_idx["CASP8"]], self.x[self.comp_idx["NFKB"]], self.x[self.comp_idx["mTORC1"]]) 
        else:
            return AVG(self.x[self.comp_idx["CASP1"]], self.x[self.comp_idx["CASP8"]], self.x[self.comp_idx["NFKB"]])
    def compute_IL1R(self):
        return self.x[self.comp_idx["IL1"]]
    def compute_IL6(self):
        if Pneumocyte.MTORC_Scenario == 2:
            return AVG(self.x[self.comp_idx["NFKB"]], self.x[self.comp_idx["MAPK_p38"]], self.x[self.comp_idx["mTORC1"]]) 
        else:
            return AVG(self.x[self.comp_idx["NFKB"]], self.x[self.comp_idx["MAPK_p38"]])
    def compute_IL6R(self):
        return self.x[self.comp_idx["IL6"]]
    def compute_IRF3(self):
        return self.x[self.comp_idx["RIG1"]]
    def compute_IKKB_a_b(self):
        return AVG(self.x[self.comp_idx["TLR4"]], self.x[self.comp_idx["IL1R"]],
                   self.x[self.comp_idx["TNFR"]], self.x[self.comp_idx["AKT"]])
    def compute_ISG(self):
        if Pneumocyte.SCENARIO == 2:
            return self.x[self.comp_idx["STAT1"]] - self.x[self.comp_idx["Viral_Repl"]] 
        else:
            return self.x[self.comp_idx["STAT1"]]
    def compute_MAPK_p38(self):
        return AVG(self.x[self.comp_idx["ANG_2_T1R"]], self.x[self.comp_idx["TLR4"]], self.x[self.comp_idx["ROS"]])
    def compute_mTORC1(self):
        if Pneumocyte.MTORC_knockout == True:
            return 0
        return AVG(NOT(self.x[self.comp_idx["TSC2"]]), NOT(self.x[self.comp_idx["FOXO3A"]]))
    def compute_mTORC2(self):
        return AVG(self.x[self.comp_idx["PI3K"]], self.x[self.comp_idx["AMPK"]])
    def compute_NFKB(self):
        return AVG(self.x[self.comp_idx["IKKB_a_b"]], self.x[self.comp_idx["ROS"]]) - self.x[self.comp_idx["FOXO3A"]]
    def compute_NLRP3(self):
        return AVG(self.x[self.comp_idx["NFKB"]], self.x[self.comp_idx["RIG1"]])
    def compute_Nutr_Depr(self):
        return self.x[self.comp_idx["Nutr_Depr"]]
    def compute_p53(self):
        return AVG(self.x[self.comp_idx["Hypoxia"]], self.x[self.comp_idx["AMPK"]]) - self.x[self.comp_idx["Viral_Repl"]]
    def compute_PI3K(self):
        return AVG(self.x[self.comp_idx["TLR4"]], self.x[self.comp_idx["ROS"]],
                   self.x[self.comp_idx["IL6R"]])
    def compute_PTEN(self):
        return self.x[self.comp_idx["p53"]]
    def compute_RIG1(self):
        return self.x[self.comp_idx["Viral_Repl"]]
    def compute_ROS(self):
        return AVG(self.x[self.comp_idx["ANG_2_T1R"]],
                  self.x[self.comp_idx["MAPK_p38"]]) - self.x[self.comp_idx["FOXO3A"]]
    def compute_SIL6R(self):
        return AVG(self.x[self.comp_idx["ADAM_17"]], self.x[self.comp_idx["IL6R"]])
    def compute_STAT1(self):
        return self.x[self.comp_idx["IFNR"]]
    def compute_STAT3(self):
        return self.x[self.comp_idx["IL6R"]]
    def compute_TLR4(self):
        return self.x[self.comp_idx["Virus"]]
    def compute_TNF(self):
        if Pneumocyte.MTORC_Scenario == 2: # one of the cytokines
            return AVG(self.x[self.comp_idx["ADAM_17"]], self.x[self.comp_idx["NFKB"]], self.x[self.comp_idx["MAPK_p38"]], self.x[self.comp_idx["mTORC1"]]) 
        else:
            return AVG(self.x[self.comp_idx["ADAM_17"]], self.x[self.comp_idx["NFKB"]], self.x[self.comp_idx["MAPK_p38"]])
    def compute_TNFR(self):
        return self.x[self.comp_idx["TNF"]]
    def compute_TSC2(self):
        return AVG(self.x[self.comp_idx["p53"]], self.x[self.comp_idx["AMPK"]], self.x[self.comp_idx["HIF_1a"]]) - self.x[self.comp_idx["AKT"]]
    def compute_Viral_Repl(self):
        return AVG(self.x[self.comp_idx["Virus"]], self.x[self.comp_idx["mTORC1"]]) - self.x[self.comp_idx["ISG"]]
    def compute_Virus(self):
        return self.x[self.comp_idx["Virus"]]

    def set_oper(self, component):
        return lambda: getattr(self, f'compute_{component}')()

    def update_state(self, x, indices, update_funcs):
        """
        Updates the states for multiple indices using provided update functions in a vectorized manner.

        Args:
            x (np.ndarray): Array of current states.
            indices (np.ndarray): Array of indices to update.
            update_funcs (list): List of update functions corresponding to each index.

        Returns:
            np.ndarray: Updated states.
        """

        for idx, update_func in zip(indices, update_funcs):
            new_value = update_func()  # This will return the computed value, not modify self.x directly

            if new_value > x[idx]:
                x[idx] += 1
            elif new_value < x[idx]:
                x[idx] -= 1

            # Ensure the new state is within bounds
            x[idx] = max(0, min(MAX_STATES, x[idx]))

        return x

def QN(components, n_iter=N_ITER, orders=ORDERS, MTORC_Scenario=0, IFN_Scenario=0, nutr_depr=0, hypoxia=0, MTORC_knockout=False, IFN_knockout=False):
    """
    QN model.

    Args:
    components (list): List of components.
    n_iter (int): Number of iterations.
    orders (int): Number of models averaged.
    MTORC_Scenario (int): MTORC scenario.
    IFN_Scenario (int): IFN scenario.
    nutr_depr (int): Nutrient deprivation.
    hypoxia (int): Hypoxia.
    MTORC_knockout (bool): MTORC knockout.
    IFN_knockout (bool): IFN knockout.

    Returns:
    np.ndarray: Matrix of results.
    """

    cell = Pneumocyte()
    operations = {component: cell.set_oper(component) for component in components}
    
    for component in components:
        if component not in operations:
            raise NotImplementedError(f"Operation for '{component}' not defined in 'operations'")

    temp_mat = np.zeros((n_iter, len(components), orders))

    for order in range(orders):
        cell.x = np.zeros(len(components))
        cell.x[cell.comp_idx["Virus"]] = MAX_STATES
        cell.x[cell.comp_idx["ACE2"]] = MAX_STATES
        cell.x[cell.comp_idx["Nutr_Depr"]] = nutr_depr
        cell.x[cell.comp_idx["Hypoxia"]] = hypoxia
        Pneumocyte.MTORC_Scenario = MTORC_Scenario
        Pneumocyte.SCENARIO = IFN_Scenario
        Pneumocyte.MTORC_knockout = MTORC_knockout
        Pneumocyte.IFN_knockout = IFN_knockout

        for i in range(n_iter):
            indices = np.arange(len(components))
            np.random.shuffle(indices)

            update_funcs = [operations[components[idx]] for idx in indices]
            cell.x = cell.update_state(cell.x, indices, update_funcs)
            temp_mat[i, :, order] = cell.x

    return temp_mat

def main():
    '''
    Runs the experiment with different conditions and scenarios.
    '''
    nutr_cond = ["healthy", "starvation"]
    hyp_cond = ["normoxia", "hypoxia"]

    conditions = [
    {"use_hyp_cond": False, "use_nutr_cond": False, "MTORC_knockout": False, "IFN_knockout": False},
    {"use_hyp_cond": True,  "use_nutr_cond": False, "MTORC_knockout": False, "IFN_knockout": False},
    {"use_hyp_cond": False, "use_nutr_cond": True,  "MTORC_knockout": False, "IFN_knockout": False},
    {"use_hyp_cond": False, "use_nutr_cond": False, "MTORC_knockout": True,  "IFN_knockout": False},
    {"use_hyp_cond": False, "use_nutr_cond": False, "MTORC_knockout": False, "IFN_knockout": True}
    ]
    MTORC_Scenarios = [0, 1, 2]
    IFN_Scenarios = [0, 1, 2]

    components = ["ACE2", "ADAM_17", "AKT", "AMPK", "ANG_2", "ANG_2_T1R",
                  "Apoptosis", "BCL_2", "CASP1", "CASP8", "CASP9", "C_FLIP",
                  "CREB_1", "FADD", "FOXO3A", "HIF_1a", "Hypoxia", "IFN_a_b",
                  "IFNR", "IL1", "IL1R", "IL6", "IL6R", "IRF3", "IKKB_a_b",
                  "ISG", "MAPK_p38", "mTORC1", "mTORC2", "NFKB", "NLRP3",
                  "Nutr_Depr", "p53", "PI3K", "PTEN", "RIG1", "ROS", "SIL6R",
                  "STAT1", "STAT3", "TLR4", "TNF", "TNFR", "TSC2", "Viral_Repl",
                  "Virus"]
    comp = ["APO", "IFN", "IL6", "TNF", "Vir.rep"]
    # comp = components
    comp_idx = [6, 17, 21, 41, 44] # which columns to use (corresponds to comp)
    # comp_idx = list(range(0,len(components)))
    rows = [45] # which iteration / row to print
    n_trials = len(conditions)

    results_by_condition = []

    for trial in range(n_trials):
        cond = conditions[trial]
        nutr_depr = MAX_STATES if cond["use_nutr_cond"] else 0
        hypoxia = MAX_STATES if cond["use_hyp_cond"] else 0

        for MTORC_Scenario in MTORC_Scenarios:
            for IFN_Scenario in IFN_Scenarios:
                condition_name = f"{nutr_cond[nutr_depr]}_{hyp_cond[hypoxia]}_MTORC_{cond['MTORC_knockout']}_IFN_{cond['IFN_knockout']}_MTORC_Scenario_{MTORC_Scenario}_IFN_Scenario_{IFN_Scenario}"
                # print(f"Running: {condition_name}")

                results_matrix = np.zeros((SIZE, len(rows), len(comp_idx)))

                for sample in range(SIZE):
                    sim_results = QN(components, N_ITER, ORDERS, MTORC_Scenario=MTORC_Scenario, IFN_Scenario=IFN_Scenario, 
                                     nutr_depr=nutr_depr, hypoxia=hypoxia, MTORC_knockout=cond["MTORC_knockout"], 
                                     IFN_knockout=cond["IFN_knockout"])
                    averaged_results = np.average(sim_results, axis=2)
                    for row_idx, row in enumerate(rows):
                        results_matrix[sample, row_idx, :] = averaged_results[row - 1, comp_idx]

                results_by_condition.append({
                    "Nutrient_Condition": nutr_cond[nutr_depr],
                    "Hypoxia_Condition": hyp_cond[hypoxia],
                    "MTORC_Knockout": cond["MTORC_knockout"],
                    "IFN_Knockout": cond["IFN_knockout"],
                    "MTORC_Scenario": MTORC_Scenario,
                    "IFN_Scenario": IFN_Scenario,
                    "results_matrix": results_matrix
                })

                print_results(condition_name, results_matrix, comp, rows)
    save_results(results_by_condition, comp, rows)

def print_results(condition_name, results_matrix, components, rows):
    """
    Prints the simulation results.

    Args:
    condition_name (str): Name of the condition.
    results_matrix (np.ndarray): Matrix of results.
    components (list): List of components.
    rows (list): List of row indices.
    """
    for comp_idx, component in enumerate(components):
        for row_idx, row in enumerate(rows):
            values = results_matrix[:, row_idx, comp_idx]
            values_str = ", ".join(map(str, values))
            print(f"{condition_name}_{component}_{row} = c({values_str})")

def save_results(results_by_condition, comp, rows):
    """
    Saves the simulation results to a csv file, including the MTORC condition (healthy or starvation), hypoxia condition 
    (normoxia or hypoxia), the IFN Knockout, the MTORC knockout, and the model and scenario numbers.
    """

    with open("results.csv", "w") as file:
        file.write("Nutrient_Condition,Hypoxia_Condition,MTORC_Knockout,IFN_Knockout,MTORC_Scenario,IFN_Scenario,Component,Iteration,Value\n")

        for result in results_by_condition:
            for comp_idx, component in enumerate(comp):
                for row_idx, row in enumerate(rows):
                    values = result["results_matrix"][:, row_idx, comp_idx]
                    for i, value in enumerate(values):
                        file.write(f"{result['Nutrient_Condition']},{result['Hypoxia_Condition']},"
                                   f"{result['MTORC_Knockout']},{result['IFN_Knockout']},"
                                   f"{result['MTORC_Scenario']},{result['IFN_Scenario']},{component},{row},{value}\n")       
                
def run():
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import seaborn as sns

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
    mat = np.average(QN(components, N_ITER, ORDERS, N_STATES, MTORC_Scenario=2, IFN_Scenario=2, nutr_depr=0), axis=2)
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
    main()
    print("Complete")