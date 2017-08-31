"""
plotly plot, use python-igraph
"""
import os
import sys
import time
import json
import igraph as ig
# import plotly.plotly as py
import plotly.offline as py
from plotly.graph_objs import Figure, Data, Font, Annotation, Annotations, \
    Margin, XAxis, YAxis, ZAxis, Scene, Layout, Line, Scatter3d, Marker


def plotly_plot_network(file_dir):
    """
    plotly test
    """
    data = []
    with open(os.path.join(file_dir, "output", "miserables.json"), 'r') as f_h:
        data = json.loads(f_h.read())

    n_n = len(data['nodes'])
    # Define the list of edges and the Graph object from edges
    n_l = len(data['links'])
    edges = [(data['links'][k]['source'], data['links'][k]['target'])
             for k in range(n_l)]

    g_o = ig.Graph(edges, directed=False)

    labels = []
    group = []
    for node in data['nodes']:
        labels.append(node['name'])
        group.append(node['group'])

    # Get the node positions, set by the Kamada-Kawai layout for 3D graphs
    layt = g_o.layout('kk', dim=3)

    # Set data for the Plotly plot of the graph
    x_n = [layt[k][0] for k in range(n_n)]  # x-coordinates of nodes
    y_n = [layt[k][1] for k in range(n_n)]  # y-coordinates
    z_n = [layt[k][2] for k in range(n_n)]  # z-coordinates
    x_e = []
    y_e = []
    z_e = []
    for e in edges:
        # x-coordinates of edge ends
        x_e += [layt[e[0]][0], layt[e[1]][0], None]
        y_e += [layt[e[0]][1], layt[e[1]][1], None]
        z_e += [layt[e[0]][2], layt[e[1]][2], None]

    trace1 = Scatter3d(x=x_e,
                       y=y_e,
                       z=z_e,
                       mode='lines',
                       line=Line(color='rgb(125,125,125)', width=1),
                       hoverinfo='none'
                       )
    trace2 = Scatter3d(x=x_n,
                       y=y_n,
                       z=z_n,
                       mode='markers',
                       name='actors',
                       marker=Marker(symbol='dot',
                                     size=6,
                                     color=group,
                                     colorscale='Viridis',
                                     line=Line(
                                         color='rgb(50,50,50)', width=0.5)
                                     ),
                       text=labels,
                       hoverinfo='text'
                       )

    axis = dict(showbackground=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''
                )
    layout = Layout(
        title="3D visualization of Network",
        width=1000,
        height=1000,
        showlegend=False,
        scene=Scene(
            xaxis=XAxis(axis),
            yaxis=YAxis(axis),
            zaxis=ZAxis(axis),
        ),
        margin=Margin(
            t=100
        ),
        hovermode='closest',
        annotations=Annotations([
            Annotation(
                showarrow=False,
                text="N-Propane: <a href='https://github.com/AdamPI314/catalytic_cycle'>[1] catalytic cycle</a>",
                xref='paper',
                yref='paper',
                x=0,
                y=0.1,
                xanchor='left',
                yanchor='bottom',
                font=Font(
                    size=14
                )
            )
        ]),)

    fig_data = Data([trace1, trace2])
    fig = Figure(data=fig_data, layout=layout)

    # help(py.plot)
    # plot(figure_or_data, show_link=True, link_text='Export to plot.ly', validate=True, output_type='file', include_plotlyjs=True,
    #      filename='temp-plot.html', auto_open=True, image=n_none, image_filename='plot_image', image_width=800, image_height=600, config=n_none)

    py.plot(fig, show_link=False, auto_open=False,
            filename=os.path.join(file_dir, "output", "n_les-Miserables.html"))

    return


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)
    I_TIME = time.time()

    plotly_plot_network(FILE_DIR)    

    E_TIME = time.time()
    print("Time elapsed: {0:.5} seconds".format(E_TIME - I_TIME))
