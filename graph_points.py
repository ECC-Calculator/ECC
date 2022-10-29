import pandas as pd
import plotly_express as px

def plot(x_coordinates, y_coordinates):
    df = pd.DataFrame(dict(y = y_coordinates, x = x_coordinates))
    graph = px.scatter(df, x = "x", y = "y", title = "Points on EC", labels = dict(y="Y-axis", x="X-axis"))
    graph.show()
    