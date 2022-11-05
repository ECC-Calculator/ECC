import os
import pandas as pd
import plotly_express as px

def plot(x_coordinates, y_coordinates):
    if not os.path.exists("images"):
        os.mkdir("images")

    df = pd.DataFrame(dict(y = y_coordinates, x = x_coordinates))
    fig = px.scatter(df, x = "x", y = "y", title = "Points on EC", labels = dict(y="Y-axis", x="X-axis"))
    '''
    #for smaller graph with difference of 1 between coordinates
    fig.update_xaxes(dtick = 1)
    fig.update_yaxes(dtick = 1)
    '''
    # fig.show()
    fig.write_image("images/fig1.png")
    