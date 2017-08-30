import os
import time
import json
from urllib.request import urlopen
import numpy as np
import igraph as ig
# import plotly.plotly as py
import plotly.offline as py
from plotly.graph_objs import *

if __name__ == '__main__':
    I_TIME = time.time()

    data = []
    f_h = urlopen(
        "https://raw.githubusercontent.com/plotly/datasets/master/miserables.json")
    data = json.loads(f_h.read())

    print(data.keys())
    N = len(data['nodes'])
    print(N)
    # Define the list of edges and the Graph object from Edges
    L = len(data['links'])
    Edges = [(data['links'][k]['source'], data['links'][k]['target'])
             for k in range(L)]

    G = ig.Graph(Edges, directed=False)

    labels = []
    group = []
    for node in data['nodes']:
        labels.append(node['name'])
        group.append(node['group'])

    # Get the node positions, set by the Kamada-Kawai layout for 3D graphs
    layt = G.layout('kk', dim=3)
    print(np.shape(layt))

    # Set data for the Plotly plot of the graph
    Xn = [layt[k][0] for k in range(N)]  # x-coordinates of nodes
    Yn = [layt[k][1] for k in range(N)]  # y-coordinates
    Zn = [layt[k][2] for k in range(N)]  # z-coordinates
    Xe = []
    Ye = []
    Ze = []
    for e in Edges:
        # x-coordinates of edge ends
        Xe += [layt[e[0]][0], layt[e[1]][0], None]
        Ye += [layt[e[0]][1], layt[e[1]][1], None]
        Ze += [layt[e[0]][2], layt[e[1]][2], None]

    trace1 = Scatter3d(x=Xe,
                       y=Ye,
                       z=Ze,
                       mode='lines',
                       line=Line(color='rgb(125,125,125)', width=1),
                       hoverinfo='none'
                       )
    trace2 = Scatter3d(x=Xn,
                       y=Yn,
                       z=Zn,
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
        title="Network of coappearances of characters in Victor Hugo's novel<br> Les Miserables (3D visualization)",
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
                text="Data source: <a href='http://bost.ocks.org/mike/miserables/miserables.json'>[1] miserables.json</a>",
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

    data = Data([trace1, trace2])
    fig = Figure(data=data, layout=layout)

    # plot(figure_or_data, show_link=True, link_text='Export to plot.ly', validate=True, output_type='file', include_plotlyjs=True,
    #      filename='temp-plot.html', auto_open=True, image=None, image_filename='plot_image', image_width=800, image_height=600, config=None)

    py.plot(fig, show_link=False, auto_open=False,
            filename='./Les-Miserables.html')

    E_TIME = time.time()
    print("hello")
    print(E_TIME - I_TIME)
