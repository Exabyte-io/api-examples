# Exabyte API Examples

## Contents of this Repository

Below, we list the contents of this repository, in roughly the order that a user might want to go through it in order to learn how our API works.

| Folder            | Notebook                                | Description |
| ------------------|-----------------------------------------| ----------- |
| [Examples/System](examples/system/)     | [Get Authentication Params](examples/system/get_authentication_params.ipynb)               | Demonstrates how to programatically find your user ID and access token, which is to [authenticate](https://docs.exabyte.io/rest-api/authentication/) for many portions of the Exabyte API.
| [Examples/Workflow](examples/workflow/) | [Get Workflows](examples/workflow/get_workflows.ipynb)                           | Walks through how to [query](https://docs.exabyte.io/rest-api/query-structure/) the Exabyte API to programatically search for workflows. In this example, we search for workflows that calculate the total energy of a material.
| [Examples/Material](examples/material/) | [Get Materials by Formula](examples/material/get_materials_by_formula.ipynb)                | Shows how [queries](https://docs.exabyte.io/rest-api/query-structure/) can be made to search for materials stored on your account by their formula. In this example, we search for a system containing Si.
| [Examples/Material](examples/material/) | [Create Material](examples/material/create_material.ipynb)                         | Gives an overview of how materials can be generated in [JSON format](https://docs.exabyte.io/materials/data/) and uploaded to your user account. In this example, we create an FCC Si crystal and upload it.
| [Examples/Material](examples/material/) | [Import Materials from Materials Project](examples/material/import_materials_from_materialsproject.ipynb) | Demonstrates how materials can be imported from [Materials Project](https://materialsproject.org/about), if their Materials Project ID is known. In this example, we import monoclinic and hexagonal SiGe cells.
| [Examples/Material](examples/material/) | [Import Materials from Poscar](examples/material/import_materials_from_poscar.ipynb)            | Provides an example of how materials can be imported directly from Poscar files (a common chemical file format best-known [for its use in VASP](https://www.vasp.at/wiki/index.php/Input)). In this example, we import the unit cell of SiGe.
| [Examples/Job](examples/job/)            | [Create and Submit Job](examples/job/create_and_submit_job.ipynb)                   | Shows how to use the Exabyte API to [create jobs](https://docs.exabyte.io/jobs/data/) and run them on our cluster. In this example, we run a DFT calculation to get the total energy of an FCC Si unit cell using Quantum Espresso.
| [Examples/Job](examples/job/)            | [Run Simulations and Extract Properties](examples/job/run-simulations-and-extract-properties.ipynb)  | Leads you through the process of copying a [bank workflow](https://docs.exabyte.io/workflows/bank/) to your account and using it to automatically calculate the [properties](https://docs.exabyte.io/properties/overview/) of multiple materials. In this example, we determine the [band gap](https://docs.exabyte.io/properties-directory/non-scalar/band-gaps/) of Si and Ge.
| [Examples/Job](examples/job/)            | [ML - Train Model Predict Properties](examples/job/ml-train-model-predict-properties.ipynb)     | Walks you through automated dataset generation and the training/prediction of material properties using [machine learning](https://docs.exabyte.io/software-directory/overview/#machine-learning). In this example, we calculate the band gaps of Si and SiGe, and using various materials properties as descriptors, train a model to predict their band gaps. Finally, we use this trained model to predict the band gap of Ge.



## Setup

NOTE: tested with python version 3.8.6, please assert that the virtual environment is created with it.

Follow the steps below in order to setup and view the Jupyter notebooks:

0. [Install git-lfs](https://help.github.com/articles/installing-git-large-file-storage/) [[3](#links)] in order to get access to the source code and notebook files.

1. Clone repository:

    ```bash
    git clone git@github.com:Exabyte-io/exabyte-api-examples.git
    ```

    In case for some reason git-lfs was not installed at the time of cloning, the files can be pulled after installing git-lfs, through `git lfs pull`.

    Related to this, please be aware that as the `.ipynb` and `.poscar` files are stored on git-lfs, they are not part of the zip archive downloaded through GitHub's web interface.

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
    jupyter lab --config=config.py
    ```

## Usage

In order to run or edit the examples:

1. Assert an existing Exabyte.io account. Examples require an account to run. New users can register [here](https://platform.exabyte.io/register) to obtain one.

2. Open [settings](examples/settings.py) and adjust it to provide the API authentication parameters. See the [corresponding example](examples/system/get_authentication_params.ipynb) to learn how to obtain the authentication parameters.

3. Open the desired example notebook, adjust it as necessary and run.


## Contribute

This is an open-source repository and we welcome contributions for other use cases. The original set of examples is only meant to demonstrate the capabilities and can be extended.

We suggest forking this repository and introducing the adjustments there. The changes in the fork can further be considered for merging into this repository as it is commonly used on Github. This process is explained in more details elsewhere online [[4](#links)].
 
If you would like to add new examples or adjust existing ones, please consider the following:

1. Put examples into the corresponding directories by domain.

2. Walk the readers through the examples by providing step-by-step explanation similar to [this](examples/material/get_materials_by_formula.ipynb).

3. We use post-save hooks to automatically convert notebooks to python scripts. See [config](examples/config.py) file for more information. In order to facilitate code review, we exclude notebook sources from version control and store them in Git LFS [[3](#links)]. Please follow this convention.

## Links

1. Exabyte.io RESTful API, description in the online documentation: [link](https://docs.exabyte.io/rest-api/overview/)
2. Jupyter.org, official website: [link](http://jupyter.org/)
3. Git Large File Storage, official website: [link](https://git-lfs.github.com/)
4. GitHub Standard Fork & Pull Request Workflow, online explanation: [link](https://gist.github.com/Chaser324/ce0505fbed06b947d962) 
