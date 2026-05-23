import setuptools
from accomp import __version__, name

# with open("README.md", "r") as fh:
#     long_description = fh.read()

setuptools.setup(
    name = name,
    version = __version__,
    author = "Yayaati Chachan",
    author_email = "ychachan@ucsc.edu",
    description = "A package to relate composition of planets to their building blocks.",
    packages = setuptools.find_packages(),
    classifiers = (
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ),
    include_package_data=True,
    zip_safe = False,
    install_requires = [
        "numpy",
        "scipy",
        "emcee", "dynesty", "corner",
        "future", "nose", "setuptools", "configparser",
        "matplotlib"]
)
