from .ipython.packages.install import install_packages_python
from .primitive.environment import is_pyodide_environment
from .pyodide.packages.install import install_packages_pyodide


async def install_packages(notebook_name_pattern: str, config_file_path: str = "", verbose: bool = True):
    """
    Install the packages listed in config.yml for the given notebook name pattern.

    Usage in notebooks:
        from mat3ra.notebooks_utils.packages import install_packages
        await install_packages("my_notebook")

    Args:
        notebook_name_pattern (str): Pattern matched against notebook names in config.yml.
        config_file_path (str): Path to config.yml; empty string uses the JupyterLite default (/drive/config.yml).
        verbose (bool): Whether to print install progress.
    """
    if is_pyodide_environment():
        await install_packages_pyodide(notebook_name_pattern, verbose)
    else:
        install_packages_python(notebook_name_pattern, verbose)
