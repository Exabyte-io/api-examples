def install_packages_python(pkgs: str, verbose: bool = True):
    """
    Install a package in a standard Python environment.

    Args:
        pkgs (str): The name of the package to install.
        verbose (bool): Whether to print the name of the installed package.
    """
    # NOTE: in a regular Python environment packages should be installed via pip,
    # not programmatically from config.yml. Direct user to do so.
    print('To install packages, run `pip install ".[all]"` in the terminal')
