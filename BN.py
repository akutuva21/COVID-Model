import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns

def BN(components, n_iter):
    temp_mat = np.zeros((n_iter + 1, len(components)))
    x = {y : False for y in components}
    x["Virus"] = True
    x["ACE2"] = True

    for j in range(n_iter):
        r = random.sample(list(x), len(x))
        for i in r:
            if i == "Virus":
                x[i] = x["Virus"]
            if i == "Viral_Repl":
                x[i] = (x["Virus"]) and not x["ISG"]
            if i == "ACE2":
                # PKC mediates ACE2 shedding from tubular cells
                x[i] = not x["RIG1"] and not x["Virus"] # or x["RIG1"]) #x["PKC"] and not x["ADAM_17"]
                # x[i] = not x["Virus"]
                # not sure about if there is a relation since Virus just relies on ACE2 to enter cells, the presence of ACE2 promotes disease prog
            if i == "ANG_2_T1R":
                x[i] = x["ANG_2"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8193025/
                # Ang II acts on angiotensin type 1 (AT1") receptors and activates the NADPH-oxidase complex producing superoxide and promoting cell pro-oxidative and pro-inflammatory responses
            if i == "ADAM_17":
                # A2_t1r activates ADAM 17, which promotes ACE2 shedding (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8185693/")
                x[i] = x["ANG_2_T1R"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7224649/
            if i == "TLR4":
                x[i] = x["Virus"]
                # Spike glycoprotein, the major infective surface protein of SARS-CoV-2 has been found as a ligand for human TLR4
                # https://www.futuremedicine.com/doi/full/10.2217/fvl-2021-0249
            if i == "RIG1":
                x[i] = x["Viral_Repl"]
                # antiviral activity of RIG-1 may comprise inhibition of viral entry into the host cell by preventing the expression of its receptor, ACE2
                # https://www.news-medical.net/news/20210215/RIG-1-like-receptors-may-play-dominant-role-in-suppressing-SARS-CoV-2-infection.aspx
            if i == "NFKB":
                x[i] = not x["IKKB a/b"] # x["ANG_2_T1R"] or x["PKC"] or x["RIG1"] or
                # should be good, may need to find what exactly inhibits NFKB
                # common drug therapy is inhibiting NFKB (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7648206/")
            if i == "IRF3":
                x[i] = x["RIG1"]# and not x["Viral_Repl"]
                # https://journals.asm.org/doi/10.1128/CMR.00299-20
                # SARS-CoV-2 membrane protein binds to importin karyopherin subunit alpha-6 (KPNA6") to inhibit interferon regulatory factor 3(IRF3") nuclear translocation
                # https://www.frontiersin.org/articles/10.3389/fcimb.2021.766922/full
            if i == "STAT1":
                x[i] = x["IFNR"] # and not x["Viral_Repl"]
                # After the infection, STAT1 activity is inhibited by the SARS-CoV-2 proteins, NSP1, and ORF6
                # https://www.nature.com/articles/s41418-020-00633-7
            if i == "STAT3":
                x[i] = x["IL6R"] #or (x["SIL6R"] and x["IL6"]) # IL6 is the main contributor to STAT3 - add il6 and il6r to the list of components
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7937040/
            if i == "IL6":
                x[i] = x["NFKB"]
            if i == "IL6R":
                x[i] = x["IL6"]
                # https://www.kegg.jp/pathway/map05171
            if i == "SIL6R":
                x[i] = x["ADAM_17"] #or x["IL6"]
            if i == "IKKB a/b":
                x[i] = not (x["TLR4"] or x["IL1R"] or x["RIG1"] or x["TNFR"]) # all phosphorylate the complex to release ikkb
            if i == "TNF":
                x[i] = x["ADAM_17"] or x["NFKB"]
                # https://www.kegg.jp/pathway/map05171
            if i == "ISG":
                x[i] = x["STAT1"]
                # https://www.nature.com/articles/s41586-021-03234-7
            if i == "C_FLIP":
                x[i] = x["NFKB"] and not x["FOXO3A"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
            if i == "INF a/b":
                x[i] = x["IRF3"]# and not x["Viral_Repl"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7995242/
            if i == "NRLP3":
                x[i] = x["NFKB"]
                # https://www.nature.com/articles/ni.3772
            if i == "CASP1":
                x[i] = x["NRLP3"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6651423/
            if i == "IFNR":
                x[i] = x["INF a/b"]
                # need to confirm more
                # https://www.frontiersin.org/articles/10.3389/fimmu.2020.606456/full
            if i == "TBid":
                x[i] = x["CASP8"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882451/
            if i == "Bax_Bak":
                x[i] = not x["BCL_2"] and x["TBid"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
            if i == "CASP9":
                x[i] = x["Bax_Bak"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
            if i == "TNFR":
                x[i] = x["TNF"]
                # Do more research later
                # https://www.frontiersin.org/articles/10.3389/fimmu.2020.585880/full
            if i == "FADD":
                x[i] = x["RIG1"]
                # Do more research on this later
                # https://pubmed.ncbi.nlm.nih.gov/9430227/
            if i == "Pyroptosis":
                x[i] = x["CASP1"]
                # https://www.nature.com/articles/s41467-019-09753-2
            if i == "IL1":
                x[i] = x["MLKL"] or x["NFKB"]
                # Look into this more later
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
            if i == "IL1R":
                x[i] = x["IL1"]
            if i == "MLKL":
                x[i] = x["RIPK1&3"] #and not (x["NRLP3"]")# or x["CASP1"]")
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
                # read more of this later
            if i == "Necroptosis":
                x[i] = x["MLKL"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
            if i == "RIPK1&3":
                x[i] = (x["RIG1"] or x["TLR4"] or x["STAT1"]) and not x["CASP8"]
                # read more later
                # https://www.sciencedirect.com/science/article/pii/S1097276514008661
            if i == "Apoptosis":
                x[i] = x["CASP8"]
                # https://pubmed.ncbi.nlm.nih.gov/10200555/
            
            ## NEW ###
            if i == "ANG_2":
                # ACE2 converts Ang II into Ang-(1–7")
                x[i] = x[i]
            if i == "ANG_1_7":
                x[i] = x["ACE2"] and x["ANG_2"]
            if i == "ROS":
                x[i] = x["ANG_2"] and not x["ANG_1_7"] #not x["FOXO3A"]
            if i == "RTK":
                x[i] = x["Virus"]
            if i == "PI3K":
                x[i] = x["RTK"]
            if i == "AKT":
                x[i] = x["PI3K"] or x["mTORC2"] or x["TLR4"] # maybe tlr4, not sure
            if i == "FOXO3A":
                #x[i] = x["Virus"]
                x[i] = (x["ROS"] or x["STAT3"]) and not (x["IKKB a/b"] or x["AKT"])
                # no direct relation (drug targetting in place - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8187014/")
            if i == "TSC2":
                x[i] = not x["AKT"] and not x["RTK"]
            if i == "mTORC1":
                x[i] = not x["TSC2"]
            if i == "PKC":
                x[i] = x["ANG_2_T1R"] or x["mTORC2"]
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9000463/
            if i == "mTORC2":
                x[i] = x["PI3K"]
            if i == "CASP8":
                x[i] = x["FADD"] or x["ROS"] or not x["AKT"] # not sure about the or not x["AKT"]
                # or not x["C_FLIP"]
                # too many, will assume all is true, only certain drugs inhibit casp 8
            if i == "BCL_2":
                x[i] = x["NFKB"] and not x["AKT"] # not sure about the and not x["AKT"]
                # https://www.nature.com/articles/1204926
        temp_mat[j, :]=[int(a) for a in x.values()]
    temp_mat[n_iter, :] = np.average(temp_mat[0:n_iter - 1,:],axis=0)
    return temp_mat

def main():
    random.seed(0)
    components = ["Virus","Viral_Repl","ACE2","PKC","ANG_2_T1R","ANG_2","ANG_1-7","ADAM_17","SIL6R","TLR4","RIG1","NFKB","IKKB a/b", "TNF","IRF3","STAT1","STAT3", "IL6", "IL6R","ISG","C_FLIP","INF a/b","NRLP3","CASP1","FOXO3A","IFNR","BCL_2","TBid","Bax_Bak","CASP9","ROS","TNFR","FADD","Pyroptosis","IL1","IL1R","MLKL","Necroptosis","RIPK1&3","CASP8","Apoptosis","RTK","PI3K","AKT","TSC2","mTORC1","mTORC2"]
    #components.sort()

    # find experiments that do this to confirm (grid search of genes array)

    n_iter = 25
    orders = 100
    mat = np.zeros((n_iter + 1, len(components)))
    # consider simplification in the future, add negative feedback loops

    for k in range(orders):
        mat += BN(components, n_iter)
    mat /= orders

    fig = plt.figure()
    fig.set_size_inches(10, 7, forward=True)
    yticklabels = [str(x) for x in range(1,n_iter + 1)]
    yticklabels.append("Average")

    cmap = cm.get_cmap('icefire_r')
    # def truncate_colormap_1(cmap, min1=0.0, max1=1.0, min2=1.0, max2=1.0, n=100:
    #     new_cmap = colors.LinearSegmentedColormap.from_list('truncated_%s' % cmap.name,
    #         cmap(np.concatenate((np.linspace(min1, max1, n), np.linspace(min2, max2, n)), axis = None)))
    #     return new_cmap
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list('truncated_%s' % cmap.name,
                                                            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    edited_comp = [x.replace("_", " ").replace("a/b", "α/β") for x in components]
    ax = sns.heatmap(mat, cmap = truncate_colormap(cmap,0.25, 0.75,n=200), linewidths=.05, xticklabels=edited_comp,
                     yticklabels=yticklabels, vmin=0, vmax=1, alpha = 0.7)
    ax.tick_params(axis='y', which='major', labelsize= 10)

    colorbar = ax.collections[0].colorbar
    xmin, xmax, delta = 0, 1, 0.1
    colorbar.set_ticks(np.arange(xmin, xmax + delta, delta))

    ax.set_xlabel('Model Component', fontsize=12)
    ax.set_ylabel('Iteration Number', fontsize=12)
    ax.set_title(f'Model Component Activation in COVID-19 with {orders} Samples', fontsize=14)

    plt.tight_layout()
    #plt.show()
    fig.savefig('plot.svg')
    fig.savefig('plot', format='png')

if __name__ == '__main__':
    main()
    print("Complete")