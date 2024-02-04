from IPython.display import display, Javascript
import json
import time


def set_data(key, value):
    """
    This function takes a Python object, serializes it to JSON, and sends it to the host environment
    through a JavaScript function defined in the JupyterLite extension `data_bridge`.

    Args:
        materials (object): The Python object to be sent to the host environment.
    """
    serialized_data = json.dumps({[key]: value})
    js_code = f"""
    (function() {{
        window.sendDataToHost({serialized_data})
        console.log({serialized_data})
    }})();
    """

    display(Javascript(js_code))
    print("Status: {key} sent to host.")


def get_data(key):
    """
    This function requests materials from the host environment through a JavaScript function defined in the JupyterLite
    extension `data_bridge`. The materials are then returned to the Python environment.
    """
    js_code = """
    (function() {
        if (window.requestDataFromHost) {
            window.requestDataFromHost();
            
        } else {
            console.error('requestDataFromHost function is not defined on the window object.');
        }
    })();
    """

    display(Javascript(js_code))
    time.sleep(2)  # JS postMessage is asynchronous, so we need to wait for the response from JS host
    print("Status: {key} received")
