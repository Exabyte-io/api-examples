import ipywidgets as widgets
from IPython.display import display, clear_output

# Parameters for the
# as defined in: https://pymatgen.org/pymatgen.analysis.interfaces.html#pymatgen.analysis.interfaces.zsl.ZSLGenerator
ZSL_GENERATOR_PARAMS = {
    "MAX_AREA": 400,
    "MAX_AREA_TOL": 0.09,
    "MAX_LENGTH_TOL": 0.03,
    "MAX_ANGLE_TOL": 0.01,
    "STRAIN_TOL": 10e-6,
}

INTERFACE_PARAMS = {
    "SUBSTRATE_INDEX": 0,
    "LAYER_INDEX": 1,
    "SUBSTRATE_MILLER": (1, 1, 1),
    "SUBSTRATE_THICKNESS": 1,
    "LAYER_MILLER": (0, 0, 1),
    "LAYER_THICKNESS": 1,
    "DISPLACEMENT_X": 0.0,
    "DISPLACEMENT_Y": 0.0,
    "DISPLACEMENT_Z": 3.0,
}

PLOT_SETTINGS = {
    "X_MIN": 0.01,  # percentage
    "X_MAX": 100,  # percentage
    "Y_MIN": 1,  # number of atoms
    "Y_MAX": 1000,  # number of atoms
    "X_SCALE": "log",
    "Y_SCALE": "log",
}


def get_zsl_generator_params():
    return ZSL_GENERATOR_PARAMS


def get_interface_params():
    return INTERFACE_PARAMS


def get_plot_settings():
    return PLOT_SETTINGS


def create_widgets(params):
    widgets_dict = {}
    for key, val in params.items():
        if isinstance(val, tuple):
            widgets_dict[key] = widgets.HBox(
                [
                    widgets.FloatText(value=v, description=f"{key}_{i+1}", continuous_update=False)
                    for i, v in enumerate(val)
                ]
            )
        else:
            widgets_dict[key] = widgets.FloatText(value=val, description=key, continuous_update=False)
    return widgets_dict


def display_form(zsl_widgets, interface_widgets, update_func):
    for widget in zsl_widgets.values():
        display(widget)
    for widget in interface_widgets.values():
        display(widget)
    update_button = widgets.Button(description="Update Parameters")
    update_button.on_click(update_func)
    display(update_button)


def update_params(button):
    for key, widget in zsl_widgets.items():
        ZSL_GENERATOR_PARAMS[key] = widget.value
    for key, widget in interface_widgets.items():
        if isinstance(widget, widgets.HBox):
            INTERFACE_PARAMS[key] = tuple(w.value for w in widget.children)
        else:
            INTERFACE_PARAMS[key] = widget.value
    clear_output(wait=True)
    print("Parameters updated:")
    print("ZSL_GENERATOR_PARAMS:", ZSL_GENERATOR_PARAMS)
    print("INTERFACE_PARAMS:", INTERFACE_PARAMS)


zsl_widgets = create_widgets(ZSL_GENERATOR_PARAMS)
interface_widgets = create_widgets(INTERFACE_PARAMS)
