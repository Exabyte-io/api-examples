# Exabyte API Examples

This repository explains how to perform some of the most common tasks in the Exabyte.io platform through its RESTful application programming interface (REST API) [[1](#links)] by virtue of examples. Examples are grouped together by domain (eg. "materials") and are presented in a self-documented format inside Jupyter notebooks [[2](#links)]. In order to view the content online, [navigate](examples/) to a notebook page inside this repository.

## Setup

NOTE: testied with python version 3.8.6.

Follow the steps below in order to setup and view the Jupyter notebooks:

0. [Install git-lfs](https://help.github.com/articles/installing-git-large-file-storage/) [[3](#links)] in order to get access to the source code and notebook files.

1. Clone repository:

    Mac / Unix:
    
    ```bash
    git clone git@github.com:Exabyte-io/exabyte-api-examples.git
    ```
   
   Windows:
   
   On windows machines (including the CygWin shell), symlinks are not cloned by default. Thus, the following command must be used:
   
   ```powershell
   git clone -c core.symlinks=True git@github.com:Exabyte-io/exabyte-api-examples.git
   ```

    In case for some reason git-lfs was not installed at the time of cloning, the files can be pulled after subsequent installation later, through `git lfs pull`.

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

2. Open [settings](examples/settings.py) and adjust it to provide the API authentication parameters. See the [corresponding example](examples/system/get_authentication_params.ipynb) to learn how to obtain the authentication parameters.

3. Open the desired example notebook, adjust it as necessary and run.


## Contribute

This is an open-source repository and we welcome contributions for other use cases. The original set of examples is only meant to demonstrate the capabilities and can be extended.

We suggest forking this repository and introducing the adjustments there. The changes in the fork can further be considered for merging into this repository as it is commonly used on Github. This process is explained in more details elsewhere online [[4](#links)].
 
If you would like to add new examples or adjust existing ones, please consider the following:

1. Put examples into the corresponding directories by domain.

2. Walk the readers through the examples by providing step-by-step explanation similar to [this](examples/material/get_materials_by_formula.ipynb).

3. We use post-save hooks to automatically convert notebooks to python scripts. See [config](config.py) file for more information. In order to facilitate code review, we exclude notebook sources from version control and store them in Git LFS [[3](#links)]. Please follow this convention.

## Links

1. Exabyte.io RESTful API, description in the online documentation: [link](https://docs.exabyte.io/rest-api/overview/)
2. Jupyter.org, official website: [link](http://jupyter.org/)
3. Git Large File Storage, official website: [link](https://git-lfs.github.com/)
4. GitHub Standard Fork & Pull Request Workflow, online explanation: [link](https://gist.github.com/Chaser324/ce0505fbed06b947d962) 
