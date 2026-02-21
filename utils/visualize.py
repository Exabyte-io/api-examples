import io
import json
import os
import time
import uuid
from enum import Enum
from typing import Dict, List, Optional, Tuple, Union

import ipywidgets as widgets
from ase.build import make_supercell
from ase.io import write
from IPython.display import HTML, Javascript, display
from mat3ra.made.material import Material
from mat3ra.made.tools.convert import to_ase
from mat3ra.utils.array import convert_to_array_if_not
from pandas import DataFrame
from pandas.io.formats.style import Styler
from pydantic import BaseModel

from utils import settings


class ViewersEnum(str, Enum):
    wave = "wave"
    ase = "ase"


def dataframe_to_html(df: DataFrame, text_align: str = "center") -> Styler:
    """
    Converts Pandas dataframe to HTML.
    See https://pandas.pydata.org/pandas-docs/stable/style.html for more information about styling.

    Args:
        df (pd.DataFrame): Pandas dataframe.
        text_align (str): text align. Defaults to center.
    """
    styles = [
        dict(selector="th", props=[("text-align", text_align)]),
        dict(selector="td", props=[("text-align", text_align)]),
    ]
    return df.style.set_table_styles(styles)


def display_JSON(
    obj: Union[dict, list], interactive_viewer: bool = settings.use_interactive_JSON_viewer, level: int = 2
) -> None:
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

        web_dir = os.path.join(os.path.dirname(__file__), "web")

        with open(os.path.join(web_dir, "renderjson.css")) as fp:
            css = fp.read()

        with open(os.path.join(web_dir, "renderjson.js")) as fp:
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


def get_material_image(material: Material, title: str, rotation="0x,0y,0z", repetitions=[1, 1, 1]):
    """
    Returns an image of the material structure with the specified title.

    Args:
        material (Material): Material object to visualize.
        title (str): Title of the image.
        rotation (str): Rotation of the image, in degrees around the x, y, and z axes (e.g., "-90x,90y,0z").
        repetitions (list): Repetitions alongside a,b,c lattice vectors.

    Returns:
        tuple: Tuple containing the image and the title.
    """

    ase_atoms = to_ase(material)
    # Create supercell for visualization
    supercell_matrix = [[repetitions[0], 0, 0], [0, repetitions[1], 0], [0, 0, repetitions[2]]]
    material_repeat = make_supercell(ase_atoms, supercell_matrix)
    text = f"{ase_atoms.get_chemical_formula()} - {title} - rotation: {rotation}"

    # Write image to a buffer to display in HTML
    buf = io.BytesIO()
    write(buf, material_repeat, format="png", rotation=rotation)
    buf.seek(0)
    return buf.read(), text


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
    Create a responsive image grid that can display images from a specified folder.
    Ensures images are displayed at their true sizes and limits the grid to a maximum of three columns.

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


class MaterialViewProperties(BaseModel):
    repetitions: List[int] = [1, 1, 1]
    rotation: str = "0x,0y,0z"
    title: str = "Material"


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
        render_function: Name of the window function to call (e.g., 'renderThreeDEditor', 'renderResults')
        data_var_name: Variable name for the data (e.g., 'materialConfig', 'results')
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


def get_wave_viewer(material, div_id, width, height, title):
    size = min(width, height)
    html = get_viewer_html(
        div_id=div_id,
        width=size,
        height=size,
        title=title,
        custom_styles="border:1px solid #333;",
    )
    js = get_viewer_js(
        data_json=material.to_json(),
        div_id=div_id,
        bundle_url="https://exabyte-io.github.io/wave.js/main.js",
        render_function="renderThreeDEditor",
        data_var_name="materialConfig",
        css_url="https://exabyte-io.github.io/wave.js/main.css",
    )
    return html, js


def render_wave(material, properties, width=600, height=600):
    timestamp = time.time()
    div_id = f"wave-{timestamp}"
    html, js = get_wave_viewer(material, div_id, width, height, properties.title)
    display(HTML(html))
    display(Javascript(js))


def render_wave_grid(
    materials: List[Material],
    list_of_properties_configs: List[MaterialViewProperties],
    width=400,
    height=400,
    max_columns=3,
):
    html_items = []
    js_items = []
    timestamp = time.time()

    for i, material in enumerate(materials):
        properties_config = list_of_properties_configs[i]
        div_id = f"wave-{i}-{timestamp}"
        html, js = get_wave_viewer(material, div_id, width, height, properties_config.title)
        html_items.append(widgets.HTML(html))
        js_items.append(Javascript(js))

    grid = widgets.GridBox(
        html_items,
        layout=widgets.Layout(
            grid_template_columns=f"repeat({max_columns}, 1fr)",
            grid_gap="10px",
            width="100%",
        ),
    )
    display(grid)
    for js in js_items:
        display(js)


def process_material_entry(
    material_entry: Union[Material, Dict], default_properties: MaterialViewProperties
) -> Tuple[Material, MaterialViewProperties]:
    """
    Process the material entry and return the material and properties.
    Args:
        material_entry: Material or a dictionary containing the material and properties.
        default_properties: Default properties to use if not specified in the material entry.

    Returns:
        Tuple[Material, MaterialViewProperties]: Material and properties.

    """
    if isinstance(material_entry, Material):
        material = material_entry
        properties = default_properties
    elif isinstance(material_entry, dict) and "material" in material_entry:
        material = material_entry["material"]
        properties = MaterialViewProperties(
            title=material_entry.get("title", default_properties.title),
            repetitions=material_entry.get("repetitions", default_properties.repetitions),
            rotation=material_entry.get("rotation", default_properties.rotation),
        )
    else:
        raise ValueError("Invalid material entry")
    return material, properties


def visualize_materials(
    materials: Union[List[Material], List[Dict[str, Union[Material, dict]]]],
    repetitions: Optional[List[int]] = [1, 1, 1],
    rotation: Optional[str] = "0x,0y,0z",
    title: Optional[str] = "Material",
    viewer: ViewersEnum = ViewersEnum.ase,
) -> None:
    """
    Visualize the material(s) in the output cell.
    Args:
        materials: Mist of Materials or a list of dictionaries:
        {"material": Material, "title": str,"repetitions": List[int], "rotation": str}.
        repetitions (Optional[List[int]]): Repetitions alongside a, b, c lattice vectors.
        rotation (Optional[str]): Rotation of the image, in degrees around the x, y, and z axes (e.g., "-90x,90y,0z").
        title (Optional[str]): Title of the image.
        viewer (ViewersEnum): Viewer to use for visualization. "ase" or "wave".

    Returns:
        None
    """
    materials = convert_to_array_if_not(materials)

    if not materials:
        print("No materials to visualize.")
        return

    default_properties = MaterialViewProperties(
        title=title if title is not None else MaterialViewProperties.title,
        repetitions=repetitions if repetitions is not None else MaterialViewProperties.repetitions,
        rotation=rotation if rotation is not None else MaterialViewProperties.rotation,
    )

    if viewer == ViewersEnum.wave:
        wave_materials = []
        wave_properties_list = []
        for material_entry in materials:
            material, material_properties = process_material_entry(material_entry, default_properties)
            wave_materials.append(material)
            wave_properties_list.append(material_properties)
        if len(wave_materials) == 1:
            # Render single material in the wave viewer, larger size and hotkeys working
            render_wave(wave_materials[0], properties=wave_properties_list[0])
        else:
            render_wave_grid(materials=wave_materials, list_of_properties_configs=wave_properties_list)

    else:
        items = []
        for material_entry in materials:
            material, properties = process_material_entry(material_entry, default_properties)
            if material:
                image_data, image_title = get_material_image(
                    material, title=properties.title, rotation=properties.rotation, repetitions=properties.repetitions
                )
                items.append((image_data, image_title))

        display(create_responsive_image_grid(items))


def visualize_workflow(workflow, level: int = 2) -> None:
    """
    Visualize a workflow by displaying its JSON configuration.

    Args:
        workflow: Workflow object with a to_dict() method
        level: Expansion level for the JSON viewer (default: 2)

    Returns:
        None
    """
    workflow_config = workflow.to_dict()
    display_JSON(workflow_config, level=level)


def visualize_properties(results, width=900, title="Properties", extra_config=None):
    """
    Visualize properties using a Prove viewer.

    Args:
        results: List[dict] of property JSON objects (or a single dict).
        width: Container width in pixels.
        title: Title displayed above the viewer.
        extra_config: Optional dict with materials, components, callbacks, etc.
    """
    if isinstance(results, dict):
        results = [results]

    # Only keep the properties we can display
    results = [r for r in results if "value" in r or "values" in r]

    timestamp = time.time()
    div_id = f"prove-{timestamp}"
    results_json = json.dumps(results)
    extra_config_json = json.dumps(extra_config) if extra_config else "undefined"

    html = get_viewer_html(
        div_id=div_id,
        width=width,
        title=title,
        custom_styles="border:1px solid #ddd; padding:12px; background:#fff; color:#111;",
    )
    js = get_viewer_js(
        data_json=results_json,
        div_id=div_id,
        bundle_url="https://exabyte-io.github.io/prove/main.js",
        render_function="renderResults",
        data_var_name="results",
        extra_config_json=extra_config_json,
    )

    display(HTML(html))
    display(Javascript(js))
