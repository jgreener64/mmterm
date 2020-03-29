import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mmterm",
    version="0.1",
    author="Joe G Greener",
    author_email="j.greener@ucl.ac.uk",
    description="View proteins and trajectories in the terminal",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jgreener64/mmterm",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Environment :: Console :: Curses",
        "Natural Language :: English",
    ],
    license="MIT",
    keywords="protein structure model trajectory view terminal",
    scripts=["bin/mmterm"],
    install_requires=["numpy", "biopython", "drawille"],
)
