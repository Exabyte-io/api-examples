# Exabyte RESTful API Examples

This repository contains examples for performing most-common tasks in the Exabyte.io platform through its RESTful API. Examples are presented in [Jupyter Notebook](http://jupyter.org/) format.

# Setup

1. Clone the repository.
    
    ```bash
    git clone git@github.com:Exabyte-io/exabyte-api-examples.git
    ```

2. Install [virtualenv](https://virtualenv.pypa.io/en/stable/) using [pip](https://pip.pypa.io/en/stable/):

    ```bash
    pip install virtualenv
    ```

3. Create virtual environment and install required packages:

    ```bash
    cd exabyte-api-examples
    virtualenv .env
    source .env/bin/activate
    pip install -r requirements.txt
    ```

4. Run Jupyter. This will open a notebook in your browser:

    ```bash
    cd examples
    jupyter notebook --config=config.py
    ```

# Usage

1. Examples require an existing account with Exabyte.io platform. Register [here](https://platform.exabyte.io/register) to obtain one.

2. Open [settings](examples/settings.ipynb) and adjust it to provide authentication parameters. See the [corresponding example](examples/system/get_authentication_params.ipynb) for how to obtain the authentication parameters with your username and password.

3. Open the desired example notebook, adjust it as necessary and run it.


# Contribute

If you would like to add new examples or adjust existing ones, please consider the following points.

1. Put examples into the corresponding directories by domain.

2. Walk the readers through the examples by providing step-by-step explanation similar to our examples, e.g. [this](material/get_materials_by_formula.ipynb).

3. We use post-save hooks to automatically convert notebooks to python scripts. See [jupyter notebook config](config.py) for more information.
 
4. As it is difficult to review the notebooks on GitHub we [automatically](.gitattributes) add iPython notebooks to [Git LFS](https://git-lfs.github.com/).
