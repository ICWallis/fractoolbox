# from setuptools import setup, find_packages

# setup(
#     name='fractoolbox',
#     version='0.1.1',
#     packages=find_packages(),
#     install_requires=[
#         # Add your package dependencies here
#     ],
#     # entry_points={
#     #     'console_scripts': [
#     #         # Add command line scripts here
#     #     ],
#     # },
#     author='Irene Wallis',
#     author_email='irene@cubicearth.nz',
#     description='A toolbox for structural geology, borehole image analysis, and geomechanics',
#     #long_description=open('README.md').read(),
#     #long_description_content_type='text/markdown',
#     #url='https://github.com/yourusername/sharkie',
#     classifiers=[
#         'Programming Language :: Python :: 3',
#         'License :: OSI Approved :: MIT License',
#         'Operating System :: OS Independent',
#     ],
#     python_requires='>=3.6',
# )


from setuptools import find_packages, setup

# with open("app/README.md", "r") as f:
#     long_description = f.read()

setup(
    name="fractoolbox",
    version="0.0.0",
    description="A toolbox for structural geology, borehole image analysis, and geomechanics",
    # package_dir={"": "fractoolbox"},
    packages=find_packages(),
    long_description= "Long description TBC", # long_description,
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