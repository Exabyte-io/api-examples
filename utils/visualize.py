import io
from typing import Dict, List, Optional, Union

import ipywidgets as widgets
from ase.build import make_supercell
from ase.io import write
from IPython.display import display
from mat3ra.made.material import Material
from mat3ra.made.tools.convert import to_ase
from mat3ra.utils.array import convert_to_array_if_not
from pydantic import BaseModel


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
        )  # Adjust padding as needed
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


def visualize_materials(
    materials_to_view: Union[List[Material], List[Dict[str, Union[Material, dict]]]],
    repetitions: Optional[List[int]] = [1, 1, 1],
    rotation: Optional[str] = "0x,0y,0z",
    title: Optional[str] = "Material",
) -> None:
    """
    Visualize the material(s) in the output cell.
    Args:
        materials_to_view: Mist of Materials or a list of dictionaries:
        {"material": Material, "title": str,"repetitions": List[int], "rotation": str}.
        repetitions (Optional[List[int]]): Repetitions alongside a, b, c lattice vectors.
        rotation (Optional[str]): Rotation of the image, in degrees around the x, y, and z axes (e.g., "-90x,90y,0z").
        title (Optional[str]): Title of the image.

    Returns:
        None
    """
    materials_to_view = convert_to_array_if_not(materials_to_view)

    if not materials_to_view:
        print("No materials to visualize.")
        return

    default_properties = MaterialViewProperties(
        title=title if title is not None else MaterialViewProperties.title,
        repetitions=repetitions if repetitions is not None else MaterialViewProperties.repetitions,
        rotation=rotation if rotation is not None else MaterialViewProperties.rotation,
    )

    items = []
    for material_entry in materials_to_view:
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
            print("Invalid material entry:", material_entry)
            continue

        image_data, image_title = get_material_image(
            material, title=properties.title, rotation=properties.rotation, repetitions=properties.repetitions
        )
        items.append((image_data, image_title))

    display(create_responsive_image_grid(items))
