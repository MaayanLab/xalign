import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="xalign",
    version="0.1.74",
    author="Alexander Lachmann",
    author_email="alexander.lachmann@mssm.edu",
    description="Alignment in a python wrapper.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/maayanlab/xalign",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    package_data={
        "xalign": ["data/*"]
    },
    include_package_data=True,
    install_requires=[
        'pandas',
        'numpy',
        'scikit-learn',
        'progress',
        'loess',
        'tqdm',
        'statsmodels',
        'mygene',
    ],
    python_requires='>=3.6',
)