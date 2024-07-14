import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import ggplot, aes, geom_boxplot, position_dodge, theme_linedraw, theme_light, labs, geom_jitter, theme_bw, theme, element_text
import argparse

def generate_boxplots(filename='results.csv'):
    """
    Generates boxplots for the simulation results, specifically focusing on healthy and normoxia with no IFN and MTORC knockouts.

    Args:
    filename (str): Name of the file containing the simulation results.
    components (list): List of components.
    rows (list): List of row indices.
    """

    df = pd.read_csv(filename)
    df = df[(df["Nutrient_Condition"] == "healthy") & (df["Hypoxia_Condition"] == "normoxia") & 
            (df["MTORC_Knockout"] == False) & (df["IFN_Knockout"] == False)]

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    plt.rcParams['figure.dpi'] = 300

    sns.boxplot(x="Component", y="Value", data=df, ax=ax, color="white", linewidth=0.5, fliersize=1)
    sns.stripplot(x="Component", y="Value", data=df, ax=ax, color="black", size=1)

    ax.set_xlabel('Model Component', fontsize=12)
    ax.set_ylabel('Value', fontsize=12)
    ax.set_title('Model Component Activation in COVID-19', fontsize=14)

    plt.show()

def generate_boxplots_ggplot(filename = 'results.csv'):
    """
    Generates boxplots for the simulation results using ggplot.
    """

    # Read and filter the data
    df = pd.read_csv(filename)
    df = df[(df["Nutrient_Condition"] == "healthy") & 
            (df["Hypoxia_Condition"] == "normoxia") & 
            (df["MTORC_Knockout"] == False) & 
            (df["IFN_Knockout"] == False)]

    # Create the plot
    p = (ggplot(df, aes(x="Component", y="Value"))
        + geom_boxplot(color="black", fill="white", size=0.5, outlier_size=1)
        + geom_jitter(color="black", size=1, width=0.1)
        + labs(x="Model Component", y="Value", title="Model Component Activation in COVID-19")
        + theme_bw()
        + theme(
            figure_size=(9, 7),
            dpi=300,
            axis_text_x=element_text(size=12),
            axis_text_y=element_text(size=12),
            plot_title=element_text(size=14)
        ))

    # Show the plot
    p.show()

def generate_APO_test(filename='results.csv', component='APO', mtorc_scenario=-1, ifn_scenario=-1):
    """
    Generates boxplots for the simulation results using ggplot using filtered data.
    """
    df = pd.read_csv(filename)
    default_names = df.iloc[0, 0:4].values
    
    df = df[(df["Component"] == component)]
    if mtorc_scenario != -1:
        df = df[(df["MTORC_Scenario"] == mtorc_scenario)]
    if ifn_scenario != -1:
        df = df[df["IFN_Scenario"] == ifn_scenario]

    names = [""] * len(df)

    for i in range(len(df)):
        if (df.iloc[i, 0:4].values == default_names).all():
            names[i] = "WT"
        elif df.iloc[i,0] == "starvation":
            names[i] = "Nutrient_Deprivation"
        elif df.iloc[i,1] == "hypoxia":
            names[i] = "Hypoxia"
        elif df.iloc[i,2] == True:
            names[i] = "MTORC_Knockout"
        elif df.iloc[i,3] == True:
            names[i] = "IFN_Knockout"
        else:
            names[i] = "test"

    Scenario = [f"{row[4]}_{row[5]}" for row in df.itertuples(index=False, name=None)]

    # Ensure the lengths of names and scenario match Apoptosis
    if len(df) != len(names) != len(Scenario):
        raise ValueError("Length of Apoptosis data does not match the length of names")

    # Create the DataFrame
    data = pd.DataFrame({
        'variety': Scenario,
        'names': names,
    })

    data = pd.concat([data, df.reset_index(drop=True)], axis=1)
    comp = {
        "APO": "Apoptosis",
        "IFN": "IFN_a_b",
        "IL6": "IL6",
        "TNF": "TNF",
        "Vir.rep": "Viral Replication"
    }
    comp_name = comp.get(component)

    # Create the plot
    p = (ggplot(data, aes(x='names', y='Value', fill='Scenario'))
         + geom_boxplot(width=0.5, position=position_dodge(0.75), color="black")  # Boxplots with black borders and dodge for spacing
         # + scale_fill_manual(values=bw_palette)  # Use the black and white palette
         + theme_linedraw()
         + theme_light()
         + labs(x="Treatment Type", y="Degree of Activation", title=f"{comp_name} Activation in COVID-19"))

    # Show the plot
    p.save("boxplot.png", width=10, height=5, units='in')
    p.show()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Generate test plots.")
    parser.add_argument('--filename', type=str, default='results.csv', help='Path to the CSV file.')
    parser.add_argument('--component', type=str, default='APO', help='Component to filter.')
    parser.add_argument('--mtorc_scenario', type=int, default=-1, help='MTORC Scenario filter value.')
    parser.add_argument('--ifn_scenario', type=int, default=-1, help='IFN Scenario filter value.')
    
    args = parser.parse_args()
    generate_APO_test(filename=args.filename, component=args.component, 
                      mtorc_scenario=args.mtorc_scenario, ifn_scenario=args.ifn_scenario)