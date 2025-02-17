from setuptools import find_packages, setup

with open("README.md", "r") as f:
    description = f.read()

setup(
    name="fractoolbox",
    version="0.0.0",
    description="A toolbox for structural geology, borehole image analysis, and geomechanics",
    # package_dir={"": "fractoolbox"},
    packages=find_packages(),
    long_description=description, 
    long_description_content_type="text/markdown",
    url="https://github.com/ICWallis/fractoolbox",
    author="Irene Wallis",
    author_email="irene@cubicearth.nz",
    license="Apatche 2.0",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache 2.0 License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
    install_requires=[],
    extras_require={
    #    "dev": ["pytest>=7.0", "twine>=4.0.2"],
    },
    python_requires=">=3.8",
)