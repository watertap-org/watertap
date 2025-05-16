import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pyomo.environ as pyo
from pyomo.network import Arc


def plot_network(m, stream_table, path_to_save=None, figsize=(8, 6)):
    """
    Plot the network of the flowsheet.

    Args:
        m (pyomo.core.base.PyomoModel): The Pyomo model to plot.
        stream_table (pandas.DataFrame): The stream table to plot.
        path_to_save (str): The path to save the plot.

    Returns:
        None. Saves the plot to the path_to_save.
    """
    plt.figure(figsize=figsize)
    G = nx.DiGraph()

    # specify types of nodes
    included_types = [
        "_ScalarMixer",
        "_ScalarCSTR",
        "_ScalarCSTR_Injection",
        "_ScalarClarifier",
        "_ScalarSeparator",
        "_ScalarProduct",
        "_ScalarAD",
        "_ScalarFeed",
        "_ScalarPressureChanger",
        "ScalarParam",
        "_ScalarDewateringUnit",
        "_ScalarElectroNPZO",
        "_ScalarThickener",
        "OrderedScalarSet",
        "_ScalarTranslator_ADM1_ASM2D",
        "_ScalarTranslator_ASM2d_ADM1",
    ]

    for unit in m.fs.component_objects(pyo.Block, descend_into=True):
        if type(unit).__name__ in included_types:
            node_name = unit.name.split("fs.")[-1]
            G.add_node(node_name)

    # edges with mass concentrations as labels
    for arc in m.fs.component_objects(Arc, descend_into=True):
        source_name = arc.source.parent_block().name.split("fs.")[-1]
        dest_name = arc.destination.parent_block().name.split("fs.")[-1]

        G.add_edge(source_name, dest_name)

        column_mapping = {
            "thickener outlet": "thickener",
            "ADM-ASM translator outlet": "translator_asm2d_adm1",
            "dewater outlet": "dewater",
            "electroNP treated": "electroNP",
            "electroNP byproduct": "electroNP",
            "Treated water": "Treated",
            "Sludge": "Sludge",
            "Feed": "FeedWater",
        }

        # Check if source_key is one of the values in column_mapping, and if so, update source_key to the corresponding key
        source_key = next(
            (k for k, v in column_mapping.items() if v == source_name), source_name
        )

        # Construct the stream table key from source and destination (using original names)
        stream_key = f"{source_name}__{dest_name}".replace("_", "_")
        if stream_key in stream_table.columns:
            # Handle missing values in stream table
            try:
                S_N2 = np.round(
                    float(stream_table.loc["Mass Concentration S_N2", stream_key]), 3
                )
            except (ValueError, TypeError, KeyError):
                S_N2 = "-"
            try:
                S_NO3 = np.round(
                    float(stream_table.loc["Mass Concentration S_NO3", stream_key]), 3
                )
            except (ValueError, TypeError):
                S_NO3 = "-"
            try:
                S_NH4 = np.round(
                    float(stream_table.loc["Mass Concentration S_NH4", stream_key]), 3
                )
            except (ValueError, TypeError):
                S_NH4 = "-"
            try:
                S_PO4 = np.round(
                    float(stream_table.loc["Mass Concentration S_PO4", stream_key]), 3
                )
            except (ValueError, TypeError):
                S_PO4 = "-"
            try:
                S_O2 = np.round(
                    float(stream_table.loc["Mass Concentration S_O2", stream_key]), 3
                )
            except (ValueError, TypeError):
                S_O2 = "-"
            try:
                X_AUT = np.round(
                    float(stream_table.loc["Mass Concentration X_AUT", stream_key]), 3
                )
            except (ValueError, TypeError):
                X_AUT = "-"
            try:
                X_PP = np.round(
                    float(stream_table.loc["Mass Concentration X_PP", stream_key]), 3
                )
            except (ValueError, TypeError):
                X_PP = "-"
            try:
                X_PAO = np.round(
                    float(stream_table.loc["Mass Concentration X_PAO", stream_key]), 3
                )
            except (ValueError, TypeError):
                X_PAO = "-"
            # add volumetric flowrate "Volumetric Flowrate"
            try:
                V_flow = np.round(
                    float(stream_table.loc["Volumetric Flowrate", stream_key]), 3
                )
            except (ValueError, TypeError):
                V_flow = "-"
            label = f"N2: {S_N2}\nO2: {S_O2}\nNO3: {S_NO3}\nNH4: {S_NH4}\nPO4: {S_PO4}\nX_AUT: {X_AUT}\nX_PP: {X_PP}\nX_PAO: {X_PAO}\nV_flow: {V_flow}"
            G.edges[source_name, dest_name]["label"] = label

    pos = nx.kamada_kawai_layout(
        G, weight="weight", scale=3, center=None, dim=2
    )  # Increased scale for more spread
    # Adjust node positions for better layout
    position_adjustments = {
        "MX1": [-0.2, 0.2],
        "MX3": [0.2, -0.1],
        "MX4": [-0.9, -0.1],
        "R1": [0, 0.2],
        "R2": [-0.2, 0.3],
        "R5": [0, 0.3],
        "R6": [0, 0.3],
        "R7": [0, -0.2],
        "SP1": [0.1, 0],
        "SP2": [-0.3, 0.2],
        "CL": [0.0, 0.2],
        "CL2": [-0.1, 0.0],
        "P1": [-0.2, 0.0],
        "thickener": [-0.6, 0.2],
        "dewater": [-0.2, -0.1],
        "translator_adm1_asm2d": [-0.1, -0.1],
        "AD": [-0.1, 0.0],
        "translator_asm2d_adm1": [-0.4, 0],
        "FeedWater": [0.2, 0.1],
    }

    for node, [dx, dy] in position_adjustments.items():
        if node in pos:
            pos[node][0] += dx
            pos[node][1] += dy

    # Add line breaks at underscores in node labels
    labels = {node: node.replace("_", "\n") for node in G.nodes()}
    nx.draw(
        G,
        pos,
        with_labels=True,
        labels=labels,
        node_size=500,
        node_color="skyblue",
        font_size=10,
        font_weight="bold",
        arrows=True,
    )

    edge_labels = nx.get_edge_attributes(G, "label")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)

    plt.title("BSM2 + electroNP Flowsheet")
    if path_to_save is not None:
        plt.savefig(path_to_save, dpi=300)
    plt.show()
