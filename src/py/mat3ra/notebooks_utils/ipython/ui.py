import json
import os
import uuid
from typing import List

import ipywidgets as widgets
from IPython.display import HTML, display
from pandas import DataFrame
from pandas.io.formats.style import Styler
from pydantic import BaseModel

from ..settings import use_interactive_JSON_viewer

_WEB_DIR = os.path.join(os.path.dirname(__file__), "web")


def display_JSON(obj, interactive_viewer: bool = use_interactive_JSON_viewer, level: int = 2) -> None:
    """
    Displays JSON, either interactively or via a text dump to Stdout.

    The interactive viewer is based on https://github.com/mljar/mercury/blob/main/mercury/widgets/json.py.

    Args:
        obj (dict): Object to display as nicely-formatted JSON
        interactive_viewer (bool): Whether to use the interactive viewer or not
        level (int): The level to which the JSON should be expanded by default
    """
    if interactive_viewer:
        if isinstance(obj, (dict, list)):
            json_str = json.dumps(obj)
        else:
            json_str = obj

        id = str(uuid.uuid4())

        with open(os.path.join(_WEB_DIR, "renderjson.css")) as fp:
            css = fp.read()

        with open(os.path.join(_WEB_DIR, "renderjson.js")) as fp:
            js = fp.read()

        display(HTML(f'<style>{css}</style><div id="{id}"></div>'))
        display(
            HTML(
                f"<script>{js} "
                f"renderjson.set_show_to_level({str(level)}); "
                f'renderjson.set_icons("▸","▾"); '
                f'document.getElementById("{id}").appendChild(renderjson({json_str}))</script>'
            )
        )
    else:
        print(json.dumps(obj, indent=4))


def dataframe_to_html(df: DataFrame, text_align: str = "center") -> Styler:
    """
    Converts Pandas dataframe to HTML.

    Args:
        df (pd.DataFrame): Pandas dataframe.
        text_align (str): text align. Defaults to center.
    """
    styles = [
        dict(selector="th", props=[("text-align", text_align)]),
        dict(selector="td", props=[("text-align", text_align)]),
    ]
    return df.style.set_table_styles(styles)


class MaterialViewProperties(BaseModel):
    repetitions: List[int] = [1, 1, 1]
    rotation: str = "0x,0y,0z"
    title: str = "Material"


def create_image_widget(image_data, format="png", object_fit="contain"):
    """
    Creates an Image widget with specified layout settings.

    Args:
        image_data (bytes): The image data to be displayed.
        format (str): The format of the image, default is 'png'.
        object_fit (str): CSS object-fit property value, default is 'contain'.

    Returns:
        widgets.Image: A configured image widget.
    """
    image = widgets.Image(value=image_data, format=format)
    image.layout.object_fit = object_fit
    return image


def create_responsive_image_grid(image_tuples, max_columns=3):
    """
    Create a responsive image grid. Limits the grid to a maximum of three columns.

    Args:
        image_tuples (list): List of tuples where each tuple contains an image and a title.
        max_columns (int): Maximum number of columns in the grid.
    """
    items = [
        widgets.VBox(
            [
                widgets.Label(value=title, layout=widgets.Layout(height="30px", align_self="center")),
                create_image_widget(image, object_fit="contain"),
            ],
            layout=widgets.Layout(align_items="center", padding="0px 0px 10px 0px"),
        )
        for image, title in image_tuples
    ]

    column_width = f"minmax(100px, {100 / max_columns}%)"
    grid = widgets.GridBox(
        items,
        layout=widgets.Layout(
            grid_template_columns=f"repeat({max_columns}, {column_width})",
            grid_gap="10px",
            width="100%",
        ),
    )
    return grid


def get_viewer_html(div_id, width, height=None, title="Viewer", custom_styles=""):
    """
    Generate HTML container for a viewer.

    Args:
        div_id: Unique ID for the container div
        width: Width in pixels
        height: Height in pixels (optional, if not provided only width is set)
        title: Title to display above the viewer
        custom_styles: Additional inline CSS styles
    """
    if height is not None:
        size_style = f"width:{width}px; height:{height}px;"
    else:
        size_style = f"width:{width}px;"

    return f"""
    <h2>{title}</h2>
    <div id="{div_id}" style="{size_style} {custom_styles}"></div>
    """


def get_viewer_js(
    data_json,
    div_id,
    bundle_url,
    render_function,
    data_var_name="data",
    extra_config_json=None,
    css_url=None,
):
    """
    Generate JavaScript to load and render a viewer bundle.

    Args:
        data_json: JSON string of data to render
        div_id: Container div ID
        bundle_url: URL to the JS bundle
        render_function: Name of the window function to call
        data_var_name: Variable name for the data
        extra_config_json: Optional extra config as JSON string
        css_url: Optional CSS file URL to load
    """
    extra_config_arg = f", {extra_config_json}" if extra_config_json else ""
    css_loader = (
        f"""
    document.head.insertAdjacentHTML(
        'beforeend',
        '<link rel="stylesheet" href="{css_url}"/>');
    """
        if css_url
        else ""
    )

    return f"""
    const {data_var_name}={data_json};
    const container = document.getElementById('{div_id}');
    (async function() {{
        await import('{bundle_url}');
        window.{render_function}({data_var_name}, container{extra_config_arg});
    }})();
    {css_loader}
    """
