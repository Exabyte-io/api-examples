from setuptools import find_namespace_packages, setup

setup(
    name="exabyte_api_examples_utils",
    version="0.0.1",
    author='Exabyte Inc.',
    author_email='info@exabyte.io',
    url="https://github.com/Exabyte-io/exabyte-api-examples",
    description="Additional Utilities for Exabyte API Examples",
    long_description=open("exabyte_api_examples_utils/README.md").read(),
    long_description_content_type="text/markdown",
    packages=find_namespace_packages(include=["exabyte_api_examples_utils*"]),
    # packages=['exabyte_api_examples_utils'],
    # package_dir={'exabyte_api_examples_utils': 'utils'},
    package_data={"": ["**/*.poscar", "requirements.txt"]},
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords=[
        "materials science",
        "machine learning",
    ],
    install_requires=[
        "tabulate",
        "pymatgen",
    ],
    extras_require={
        "test": ["pytest", "pytest-cov"],
    },
)