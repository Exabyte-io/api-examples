import io
from typing import List, Union

import ipywidgets as widgets
from ase.build import make_supercell
from ase.io import write
from IPython.display import display
from mat3ra.made.material import Material
from mat3ra.made.tools.convert import to_ase
from mat3ra.utils.array import convert_to_array_if_not


def get_material_visualization_html(
    material: Material, title: str, rotation="0x,0y,0z", number_of_repetitions: int = 1
):
    """
    Returns an image of the material structure with the specified title.

    Args:
        material (Material): Material object to visualize.
        title (str): Title of the image.
        rotation (str): Rotation of the image.
        number_of_repetitions (int): Number of repetitions of the unit cell.

    Returns:
        tuple: Tuple containing the image and the title.
    """

    material = to_ase(material)
    # Set the number of unit cell repetition for the structure
    n = number_of_repetitions
    material_repeat = make_supercell(material, [[n, 0, 0], [0, n, 0], [0, 0, 1]])
    text = f"{material.symbols} - {title}"

    # Write image to a buffer to display in HTML
    buf = io.BytesIO()
    write(buf, material_repeat, format="png", rotation=rotation)
    buf.seek(0)
    return (buf.read(), text)


def create_responsive_image_grid(image_tuples, max_columns=3):
    """
    Create a responsive image grid that can display images from a specified folder.

    Args:
        image_tuples (list): List of tuples where each tuple contains an image and a title.
        max_columns (int): Maximum number of columns in the grid.
    """
    items = [
        widgets.VBox(
            [
                widgets.Label(value=title, layout=widgets.Layout(height="30px", overflow="hidden")),
                widgets.Image(value=image, format="png", width="auto", height="auto"),
            ]
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


def visualize(
    materials: Union[List[Material], Material],
    title: str = "Material",
    number_of_repetitions: int = 1,
    rotation="0x,0y,0z",
):
    materials = convert_to_array_if_not(materials)
    items = [
        get_material_visualization_html(
            material, title=f"{title} {i}", number_of_repetitions=number_of_repetitions, rotation=rotation
        )
        for i, material in enumerate(materials)
    ]

    display(create_responsive_image_grid(items))
