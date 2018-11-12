# Exabyte API Examples

This repository contains examples for performing most common tasks in the Exabyte.io platform through its RESTful application programming interface (REST API). Examples are grouped together by domain (eg. "materials") and are presented in a self-documented format inside [Jupyter](http://jupyter.org/). [Navigate](examples/) to a notebook page inside this repository to view its content online.

## Setup

Follow the steps below in order to setup and view the Jupyter notebooks:

1. Clone repository:
    
    ```bash
    git clone git@github.com:Exabyte-io/exabyte-api-examples.git
    ```

2. Install [virtualenv](https://virtualenv.pypa.io/en/stable/) using [pip](https://pip.pypa.io/en/stable/) if not already present:

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

4. Run Jupyter and open a notebook in a browser:

    ```bash
    cd examples
    jupyter notebook --config=config.py
    ```

## Usage

In order to run or edit the examples:

1. Assert an existing Exabyte.io account. Examples require an account to run. New users can register [here](https://platform.exabyte.io/register) to obtain one.

2. Open [settings](examples/settings.ipynb) and adjust it to provide the API authentication parameters. See the [corresponding example](examples/system/get_authentication_params.ipynb) to learn how to obtain the authentication parameters.

3. Open the desired example notebook, adjust it as necessary and run.


## Contribute

This is an open-source repository and we welcome contributions for other use cases. The original set of examples is only meant to demonstrate the capabilities and can be extended.

We suggest forking this repository and introducing the adjustments there. The changes in the fork can further be considered for merging into this repository as it is commonly used on Github. This process is explained in more details online [here](https://gist.github.com/Chaser324/ce0505fbed06b947d962), for example.
 
If you would like to add new examples or adjust existing ones, please consider the following points:

1. Put examples into the corresponding directories by domain.

2. Walk the readers through the examples by providing step-by-step explanation similar to our examples, e.g. [this](examples/material/get_materials_by_formula.ipynb).

> NOTE: We use post-save hooks to automatically convert notebooks to python scripts. See [jupyter notebook config](config.py) for more information. As it is difficult to review the notebooks on GitHub we [automatically](.gitattributes) add iPython notebooks to [Git LFS](https://git-lfs.github.com/).

## Links

1. Exabyte.io RESTful API, description in the online documentation: [link](https://docs.exabyte.io/rest-api/overview/)
