import os
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    readme = f.read()

if "COLAB_JUPYTER_IP" in os.environ:  # running in Colab
    requirements_file = "requirements-colab.txt"
else:
    requirements_file = "requirements.txt"

with open(requirements_file, "r") as f:
    requirements = [line for line in f.read().splitlines() if not line.startswith("#")]


setup(
    name="mat3ra-api-examples",
    version="0.1.0",
    description="Mat3ra API Examples",
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    package_data={"examples": ["settings.json"]},
    entry_points={
        "console_scripts": [
            "notebook-path = examples.utils.initialize_settings:get_notebook_path",
        ],
    },
    install_requires=requirements,
)
