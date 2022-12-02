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
        Number of samples to take for each iteration (default is 250)\
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
    
    temp_mat = np.zeros((n_iter + 1, len(components), orders))

    for k in np.arange(orders):
        comp = list(enumerate(components))
        x = {y: False for y in components}
        assert len(x) == len(components)
        x["Virus"] = True
        x["ACE2"] = True

        for j in np.arange(n_iter):
            np.random.shuffle(comp)
            indices, c = zip(*comp)
            for i in c:
                if i == "Virus":
                    x[i] = x["Virus"]
                elif i == "Viral_Repl":
                    x[i] = (x["Virus"]) and not x["ISG"]
                elif i == "ANG_2_T1R":
                    x[i] = x["ANG_2"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8193025/
                    # Ang II acts on angiotensin type 1 (AT1") receptors and
                    # activates the NADPH-oxidase complex producing superoxide
                    # and promoting cell pro-oxidative and pro-inflammatory
                    # responses
                elif i == "ADAM_17":
                    # A2_t1r activates ADAM 17, which promotes ACE2 shedding
                    # (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8185693/")
                    x[i] = x["ANG_2_T1R"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7224649/
                elif i == "TLR4":
                    x[i] = x["Virus"]
                    # Spike glycoprotein, the major infective surface protein
                    # of SARS-CoV-2 has been found as a ligand for human TLR4
                    # https://www.futuremedicine.com/doi/full/10.2217/fvl-2021-0249
                elif i == "RIG1":
                    x[i] = x["Viral_Repl"]
                    # AV activity of RIG-1 may inhibit of viral entry into the
                    # host cell by preventing the expression of ACE2
                    # https://www.news-medical.net/news/20210215/RIG-1-like-receptors-may-play-dominant-role-in-suppressing-SARS-CoV-2-infection.aspx
                elif i == "NFKB":
                    # x["ANG_2_T1R"] or x["PKC"] or x["RIG1"] or
                    x[i] = not x["IKKB a/b"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7648206/"
                elif i == "IRF3":
                    x[i] = x["RIG1"]  # and not x["Viral_Repl"]
                    # https://journals.asm.org/doi/10.1128/CMR.00299-20
                    # https://www.frontiersin.org/articles/10.3389/fcimb.2021.766922/full
                elif i == "STAT1":
                    x[i] = x["IFNR"]  # and not x["Viral_Repl"]
                    # After the infection, STAT1 activity is inhibited by
                    # the SARS-CoV-2 proteins, NSP1, and ORF6
                    # https://www.nature.com/articles/s41418-020-00633-7
                elif i == "STAT3":
                    # or (x["SIL6R"] and x["IL6"])
                    x[i] = x["IL6R"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7937040/
                elif i == "IL6R":
                    x[i] = x["IL6"]
                    # https://www.kegg.jp/pathway/map05171
                elif i == "SIL6R":
                    x[i] = x["ADAM_17"]  # or x["IL6"]
                elif i == "ISG":
                    x[i] = x["STAT1"]
                    # https://www.nature.com/articles/s41586-021-03234-7
                elif i == "C_FLIP":
                    x[i] = x["NFKB"] and not x["FOXO3A"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
                elif i == "NRLP3":
                    x[i] = x["NFKB"]
                    # https://www.nature.com/articles/ni.3772
                elif i == "CASP1":
                    x[i] = x["NRLP3"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6651423/
                elif i == "IFNR":
                    x[i] = x["IFN a/b"]
                    # need to confirm more
                    # https://www.frontiersin.org/articles/10.3389/fimmu.2020.606456/full
                elif i == "Bax_Bak":
                    x[i] = not x["BCL_2"] and x["tBid"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
                elif i == "CASP9":
                    x[i] = x["Bax_Bak"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
                elif i == "TNFR":
                    x[i] = x["TNF"]
                    # Do more research later
                    # https://www.frontiersin.org/articles/10.3389/fimmu.2020.585880/full
                elif i == "Pyroptosis":
                    x[i] = x["CASP1"]
                    # https://www.nature.com/articles/s41467-019-09753-2
                elif i == "IL1R":
                    x[i] = x["IL1"]
                elif i == "MLKL":
                    # and not (x["NRLP3"]")# or x["CASP1"]")
                    x[i] = x["RIPK1&3"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
                    # read more of this later
                elif i == "Necroptosis":
                    x[i] = x["MLKL"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
                elif i == "RIPK1&3":
                    x[i] = (x["RIG1"] or x["TLR4"]
                            or x["STAT1"]) and not x["CASP8"]
                    # read more later
                    # https://www.sciencedirect.com/science/article/pii/S1097276514008661
                elif i == "ANG_2":
                    # ACE2 converts Ang II into Ang-(1–7")
                    x[i] = x[i]
                elif i == "ANG_1_7":
                    x[i] = x["ACE2"] and x["ANG_2"]
                elif i == "ROS":
                    x[i] = x["ANG_2"] and not x["ANG_1_7"]  # not x["FOXO3A"]
                elif i == "RTK":
                    x[i] = x["Virus"]
                elif i == "AKT":
                    # maybe tlr4, not sure
                    x[i] = x["PI3K"] or x["mTORC2"] or x["TLR4"]
                elif i == "FOXO3A":
                    # x[i] = x["Virus"]
                    x[i] = (x["ROS"] or x["STAT3"]) and not (
                        x["IKKB a/b"] or x["AKT"])
                    # no direct relation
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8187014/"
                elif i == "TSC2":
                    x[i] = not x["AKT"] and not x["RTK"]
                elif i == "mTORC1":
                    x[i] = not x["TSC2"]
                elif i == "PKC":
                    x[i] = x["ANG_2_T1R"] or x["mTORC2"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9000463/
                elif i == "mTORC2":
                    x[i] = x["PI3K"]
                elif i == "CASP8":
                    # not sure about the or not x["AKT"]
                    x[i] = x["FADD"] or x["ROS"] or not x["AKT"]
                    # or not x["C_FLIP"]
                    # too many, will assume all is true
                    # only certain drugs inhibit casp 8
                    # https://www.nature.com/articles/1204926
                elif i == "FADD":
                    x[i] = x["TNFR"]
                    # Do more research on this later - previously was x["RIG1"]
                    # https://pubmed.ncbi.nlm.nih.gov/9430227/
                elif i == "Apoptosis":
                    # or x["CASP3"], previously was just 8
                    x[i] = x["CASP8"] or x["CASP9"]
                    # https://pubmed.ncbi.nlm.nih.gov/10200555/
                elif i == "BCL_2":
                    # previously and not x["AKT"]
                    x[i] = (x["NFKB"] or x["CREB_1"]) and not x["BAD"]
                elif i == "BAD":
                    # new
                    x[i] = not x["AKT"]
                elif i == "CREB_1":
                    x[i] = x["AKT"]
                elif i == "tBid":
                    # added ROS
                    x[i] = x["CASP8"] or x["ROS"]
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882451/
                elif i == "IFN a/b":
                    x[i] = (x["IRF3"] or x["NFKB"]) and not (
                        x["Viral_Proteins"])
                    # equivalent to Type 1 IFN
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7995242/
                elif i == "Cyt_R":
                    x[i] = x["TNF"] or x["IL6"]
                elif i == "PI3K":
                    x[i] = x["RTK"] or x["IFNR"] or x["Cyt_R"]
                elif i == "Autophagy":
                    # phenotype
                    x[i] = not x["mTORC1"]
                elif i == "Viral_Proteins":
                    # examples of viral proteins are S6K and 4E-BP1
                    x[i] = x["mTORC1"]
                
                elif i == "MAPK_p38":
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7228886/
                    # https://www.frontiersin.org/articles/10.3389/fphar.2021.631879/full
                    x[i] = (x["ANG_2_T1R"] or x["Virus"]) and not x["ANG_1_7"]
                elif i == "ACE2":
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7228886/
                    # PKC mediates ACE2 shedding from tubular cells
                    # or x["RIG1"]) #x["PKC"] and not x["ADAM_17"]
                    x[i] = not x["RIG1"] and not x["Virus"] and x["MAPK_p38"]
                    # x[i] = not x["Virus"]
                    # not sure about if there is a relation since Virus just
                    # relies on ACE2 to enter cells, the presence of ACE2
                    # promotes disease prog
                elif i == "IL6":
                    x[i] = x["NFKB"] or x["MAPK_p38"]
                elif i == "TNF":
                    x[i] = x["ADAM_17"] or x["NFKB"] or x["MAPK_p38"]
                    # https://www.kegg.jp/pathway/map05171
                elif i == "IL1":
                    x[i] = x["MLKL"] or x["NFKB"] or x["MAPK_p38"]
                    # Look into this more later
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
                elif i == "IKKB a/b":
                    # all phosphorylate the complex to release ikkb
                    # https://www.nature.com/articles/7290257
                    # https://pubmed.ncbi.nlm.nih.gov/16028365/
                    x[i] = not (x["TLR4"] or x["IL1R"]
                                or x["RIG1"] or x["TNFR"])
                else:
                    raise NotImplementedError(f'{i} was not executed (check names)')
            assert list(c) == [components[a] for a in indices]
            indices = np.sort(list(indices))
            temp_mat[j, :, k] = [x[components[a]] for a in indices]
        temp_mat[n_iter, :, k] = np.average(
            temp_mat[0:n_iter - 1, :, k], axis=0)
    return temp_mat


def main():
    np.random.seed(0)

    components = ["Virus", "Viral_Repl", "ACE2", "PKC", "ANG_2_T1R", "ANG_2",
                  "ANG_1_7", "ADAM_17", "SIL6R", "TLR4", "RIG1", "NFKB",
                  "IKKB a/b", "TNF", "IRF3", "STAT1", "STAT3", "IL6", "IL6R",
                  "ISG", "C_FLIP", "IFN a/b", "NRLP3", "CASP1", "FOXO3A",
                  "IFNR", "BCL_2", "tBid", "Bax_Bak", "CASP9", "ROS", "TNFR",
                  "FADD", "Pyroptosis", "IL1", "IL1R", "MLKL", "Necroptosis",
                  "RIPK1&3", "CASP8", "Apoptosis", "RTK", "PI3K", "AKT",
                  "TSC2", "mTORC1", "mTORC2", "BAD", "CREB_1", "Cyt_R",
                  "Autophagy", "Viral_Proteins", "MAPK_p38"]
    comp_edit = copy.deepcopy(components)
    n_iter = 25
    orders = 250
    # mat = np.zeros((n_iter + 1, len(components), orders))
    # consider simplification in the future, add negative feedback loops

    mat = np.average(BN(comp_edit, n_iter, orders), axis=2)

    yticklabels = [str(x) for x in 1 + np.arange(n_iter)]
    yticklabels.append("Average")

    cmap = cm.get_cmap('icefire_r')

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'truncated_%s' % cmap.name,
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    edited_comp = [x.replace("_", " ").replace("a/b", "α/β")
                   for x in components]
    ax = sns.heatmap(mat, cmap=truncate_colormap(cmap, 0.25, 0.75, n=200),
                     linewidths=.05, xticklabels=edited_comp,
                     yticklabels=yticklabels, vmin=0, vmax=1, alpha=1)
    ax.tick_params(axis='y', which='major', labelsize=10)

    # colorbar = ax.collections[0].colorbar
    # xmin, xmax, delta = np.min(mat), np.max(mat), 0.1
    # colorbar.set_ticks(np.arange(xmin, xmax + delta, delta))

    ax.set_xlabel('Model Component', fontsize=12)
    ax.set_ylabel('Iteration Number', fontsize=12)
    ax.set_title(
        f'Model Component Activation in COVID-19 with {orders} Samples',
        fontsize=14)

    fig = ax.get_figure()
    fig.set_size_inches(10, 7, forward=True)

    fig.tight_layout()
    # fig.show()
    # fig.savefig('plot.svg')
    fig.savefig('plot.png')


if __name__ == '__main__':
    main()
    print("Complete")
