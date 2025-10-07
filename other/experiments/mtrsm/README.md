## Usage

### With standalone scripts

Install dependencies (using pyenv):

```bash
pyenv local 3.10.13
python -m venv .venv-3.10.13
source .venv-3.10.13/bin/activate
pip install mattersim
```

Run the script(s):

```bash
source .venv-3.10.13/bin/activate
python <script_name.py>
```

### With JupyterLab

Install dependencies:

```bash
pyenv local 3.10.13
python -m venv .venv-3.10.13
source .venv-3.10.13/bin/activate
pip install "mat3ra-api-examples[forcefields]"
```

As below. With `pyenv`.

```bash
pyenv local 3.10.13
python -m venv .venv-3.10.13
source .venv-3.10.13/bin/activate
pip install mattersim
```

See the included notebook.
