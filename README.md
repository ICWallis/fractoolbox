# fractoolbox

Python tools for structural geology and borehole image analysis that includes data handling, frequency and geometric analysis, and reservoir geomechanics.  

fractoolbox is in rapid development so the content of this repo will change regularly and may not be stable. We welcome use of the fractoolbox libraries and notebooks but only on an 'as is' basis and without any warranty implied. 

The Juypter Notebooks (.ipynb) in the main directory are step-by-step walkthroughs of fractoolbox and other related libraries (such as mplstereonet). Some also describe the theory underlying the tools and point readers to relevant academic literature. If you are new to fractoolbox, these are a good place to start.

A series of Jupter Notebooks are included in the directory _published-case-study_ that describe fractoolbox methods applied to data from geothermal wells in Indonesia and New Zealand. These notebooks are designed for in-person training with geothermal operators. The paper that accompanies these Juypter Notebooks will be presented at this year's New Zealand Geothermal Workshop.     

>Wallis, I.C., Rowland, J. V. and Dempsey, D. E., Allan, G., Sidik, R., Martikno, R., McLean, K., Sihotang, M., Azis, M. and Baroek, M. 2020 (submitted) Approaches to imaging feedzone diversity with case studies from Sumatra, Indonesia, and the Taupō Volcanic Zone, New Zealand. New Zealand Geothermal Workshop: Waitangi, New Zealand.

The Juypter Notebooks in fractoolbox can be viewed as static pages in your browser window by clicking on them here. Alternatively, they may be explored in the interactive environment of mybinder.org by clicking on the button below. Binder allows you to edit the code and re-run the notebook to generate a new output, all without downloading or installing anything on your computer. Because fractoolbox is regularly changing, the binder instance may take a few minutes to load. 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ICWallis/fractoolbox/master)

### New to python? Start here.

If you're new to python, download the 'Python for Geoscientists' course from [here](https://github.com/ddempsey/python_for_geoscientists) because it includes a quick-start guide and will set you on the pythonic path.  

We recommend using [Anaconda](https://www.anaconda.com/) to manage Python packages and the Juypter Notebook installation. It comes with many of packages/libraries required by fractoolbox and is a simple way to get started. 

fractoolbox is made with Python 3.6

There are two libraries used in fractoolbox that are not already built into Annaconda. Once you you have installed Annaconda, these libraries may be installed using _anaconda prompt_ by typing:

> $ pip install mplstereonet

> https://github.com/joferkington/mplstereonet

Due to compatibility with resent versions of matplotlib, ensure you have mplstereonet version 0.6 or later.

> $ pip install seaborn

> https://seaborn.pydata.org/

When mplstereonet is used after seaborn, there may be an issue with mplstereonet drawing plot elements like the grid. If this occurs, use sns.reset_defaults() just prior to building a stereonet with mplstereonet. 

Now you're set up, use the green 'Code' button on this page to download a copy of this repository. Unzip and have fun exploring your own fracture data.

The 'fractoolbox' folder needs to be in the same folder as your code when you call it with...

> import fractoolbox as ftb

The fractoolbox module has not yet been deployed as a package so it can't be installed using pip or conda. 

### License

Apache 2.0 

https://choosealicense.com/licenses/apache-2.0/

### Development

To date, the content of fractoolbox has been generated during doctoral research undertaken by [Irene Wallis](https://www.cubicearth.nz/), which is supervised by [David Dempsey](https://sites.google.com/view/dempsey-research-group/home) and [Julie (JR) Rowland](http://www.science.auckland.ac.nz/people/profile/j-rowland)  and includes code/math contributions from [Evert Durán](https://unidirectory.auckland.ac.nz/profile/e-quintero). 

Feedback and contributions welcome.