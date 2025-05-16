"""
Tools for interactive network visualization using Dash.
"""

import socket
import threading
import time
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import networkx as nx
from dash.exceptions import PreventUpdate
from idaes.core.util.tables import stream_table_dataframe_to_string
from pyomo.network import Arc
import dash_cytoscape as cyto
import idaes.core.base.unit_model as unit_model_base
import os
import json


def find_available_port(start_port=8050, max_port=9000):
    """Find an available port starting from start_port up to max_port"""
    for port in range(start_port, max_port):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", port))
                return port
            except OSError:
                continue
    raise RuntimeError(f"No available ports found between {start_port} and {max_port}")


def run_dash_app(app, port):
    """Run Dash app in a separate thread"""
    if app is None:
        print("Error: Cannot run None app")
        return None

    def run():
        try:
            app.run(debug=False, port=port, host="127.0.0.1")
        except Exception as e:
            print(f"Error running Dash app: {str(e)}")
            import traceback

            traceback.print_exc()

    thread = threading.Thread(target=run)
    thread.daemon = True  # Thread will be killed when main program exits
    thread.start()
    time.sleep(2)  # Give the server more time to start
    return thread


def create_dash_app(model, stream_table=None, show_labels=True):
    """Create a Dash app for interactive network visualization using dash-cytoscape, with save/load layout."""
    try:
        app = dash.Dash(__name__)
        layout_file = os.path.join(os.path.dirname(__file__), "network_layout.json")
        # Create network graph
        G = nx.DiGraph()
        unit_nodes = {}
        for obj in model.fs.component_objects():
            if isinstance(obj, unit_model_base.UnitModelBlockData):
                G.add_node(obj.name)
                unit_nodes[obj.name] = obj
        for arc in model.fs.component_objects(Arc):
            source_unit = arc.source.parent_block().name
            dest_unit = arc.destination.parent_block().name
            if source_unit in unit_nodes and dest_unit in unit_nodes:
                G.add_edge(source_unit, dest_unit)
        if len(G.nodes()) == 0:
            print("Warning: No nodes found in the network!")
            return None
        # Try to load saved layout
        saved_positions = None
        if os.path.exists(layout_file):
            try:
                with open(layout_file, "r") as f:
                    saved_positions = json.load(f)
            except Exception as e:
                print(f"Could not load saved layout: {e}")
        # Use saved positions if available, else spring layout
        if saved_positions:
            pos = {
                node: (p["x"], p["y"])
                for node, p in saved_positions.items()
                if node in G.nodes()
            }
            # Fallback for any new nodes
            for node in G.nodes():
                if node not in pos:
                    pos[node] = (0, 0)
        else:
            pos = nx.spring_layout(G, k=1.5, iterations=200, seed=42)
        cy_elements = []
        for node in G.nodes():
            x, y = pos[node]
            if saved_positions:
                # Use saved positions as-is
                cy_elements.append(
                    {
                        "data": {"id": node, "label": node},
                        "position": {"x": float(x), "y": float(y)},
                        "grabbable": True,
                    }
                )
            else:
                # Scale spring layout for better spread
                cy_elements.append(
                    {
                        "data": {"id": node, "label": node},
                        "position": {"x": float(x) * 500, "y": float(y) * 500},
                        "grabbable": True,
                    }
                )
        for edge in G.edges():
            cy_elements.append(
                {
                    "data": {
                        "source": edge[0],
                        "target": edge[1],
                        "label": f"{edge[0]}→{edge[1]}",
                    }
                }
            )
        stylesheet = [
            {
                "selector": "node",
                "style": {
                    "content": "data(label)" if show_labels else "",
                    "width": 50,
                    "height": 50,
                    "background-color": "#1f77b4",
                    "color": "black",
                    "font-size": 16,
                    "text-valign": "center",
                    "text-halign": "center",
                    "z-index": 10,
                },
            },
            {
                "selector": "edge",
                "style": {
                    "curve-style": "bezier",
                    "target-arrow-shape": "triangle",
                    "line-color": "#888",
                    "target-arrow-color": "#888",
                    "width": 3,
                    "font-size": 10,
                    "color": "#888",
                    "text-background-color": "#fff",
                    "text-background-opacity": 1,
                    "text-background-shape": "roundrectangle",
                    "text-rotation": "autorotate",
                },
            },
        ]
        layout_children = [
            html.H1("Process Network Visualization"),
            cyto.Cytoscape(
                id="cytoscape-network",
                elements=cy_elements,
                layout={"name": "preset"},
                style={"width": "100%", "height": "700px", "background": "#f8f9fa"},
                stylesheet=stylesheet,
                userZoomingEnabled=True,
                userPanningEnabled=True,
                boxSelectionEnabled=True,
                autoungrabify=False,
                autolock=False,
                minZoom=0.2,
                maxZoom=2,
            ),
            html.Button("Save Layout", id="save-layout-btn", n_clicks=0),
            html.Div(
                id="save-layout-msg", style={"marginTop": "10px", "color": "green"}
            ),
        ]
        if stream_table is not None:
            layout_children.append(
                html.Div(
                    id="stream-table",
                    children=[
                        html.H3("Stream Table"),
                        html.Pre(stream_table_dataframe_to_string(stream_table)),
                    ],
                )
            )
        app.layout = html.Div(layout_children)

        # Add callback to save layout
        @app.callback(
            dash.dependencies.Output("save-layout-msg", "children"),
            [dash.dependencies.Input("save-layout-btn", "n_clicks")],
            [dash.dependencies.State("cytoscape-network", "elements")],
        )
        def save_layout(n_clicks, elements):
            if n_clicks and elements:
                # Only save node positions
                node_positions = {
                    el["data"]["id"]: el["position"]
                    for el in elements
                    if "position" in el
                }
                try:
                    with open(layout_file, "w") as f:
                        json.dump(node_positions, f, indent=2)
                    return "Layout saved!"
                except Exception as e:
                    return f"Error saving layout: {e}"
            return ""

        return app
    except Exception as e:
        print(f"Error creating Dash app: {str(e)}")
        import traceback

        traceback.print_exc()
        return None
