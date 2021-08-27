import setuptools

__version__ = "v1.0.4a"
__author__ = "S. Winkel"


with open("README.md", "r") as fh:
    long_description = fh.read()

requires = [
    "cycler",
    "kiwisolver",
    "matplotlib",
    "numpy",
    "pandas",
    "Pillow",
    "pyparsing",
    "pyproj",
    "pyshp",
    "python-dateutil",
    "pytz",
    "six",
    "scipy",
    "tk",
    "tabulate",
    "XlsxWriter",
]

setuptools.setup(
    name="breakwater",
    version=__version__,
    author=__author__,
    author_email="pybreakwater@gmail.com",
    description="A package for the conceptual design of breakwaters",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Sander-w/breakwater/",
    packages=setuptools.find_packages(include=["breakwater", "breakwater.*"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "License :: Free for non-commercial use",
    ],
    python_requires=">=3.6",
    install_requires=requires,
    keywords="conceptual design hydraulic engineering breakwaters",
    include_package_data=True,
    data_files=[("", ["LICENSE", "README.md"])],
)
