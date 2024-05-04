import base64
import io
from typing import List, Union

from ase.build import make_supercell
from ase.io import write
from IPython.display import HTML, display
from mat3ra.made.material import Material
from mat3ra.made.tools.convert import to_ase
from mat3ra.utils.array import convert_to_array_if_not


def get_material_visualization_html(
    material: Material, title: str, rotation: str = "0x", number_of_repetitions: int = 1
):
    """
    Returns an HTML string with a Base64-encoded image for visualization,
    including the name of the file, positioned horizontally.
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
    img_str = base64.b64encode(buf.read()).decode("utf-8")
    html_str = f"""
    <div style="display: inline-block; margin: 10px; vertical-align: top;">
        <p>{text}</p>
        <img src="data:image/png;base64,{img_str}" alt="{title}" />
    </div>
    """
    return html_str


def visualize(materials: Union[List[Material], Material], title: str = "Material", number_of_repetitions: int = 1):
    materials = convert_to_array_if_not(materials)
    for i, material in enumerate(materials):
        material_image = get_material_visualization_html(
            material, title=f"{title} {i}", number_of_repetitions=number_of_repetitions
        )
        display(HTML(material_image))
