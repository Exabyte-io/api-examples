from setuptools import setup, find_packages


setup(
    name="utils",
    packages=find_packages(),
    package_data={"utils": ["settings.json"]},
    entry_points={
        "console_scripts": [
            "notebook-info = utils.initialize_settings:cli",
        ],
    },
)
