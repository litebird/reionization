from setuptools import find_packages, setup

setup(
    name="REIO_forecast",
    description="Likelihood for reionization forecast",
    long_description_content_type="text/markdown",
    author="Matthieu Tristram",
    url="",
    license="GNU license",
    classifiers=[
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    packages=find_packages(),
    python_requires=">=3.5",
    install_requires=["astropy", "cobaya>=3.0"],
    package_data={"REIO_forecast": ["PLK.yaml","SO.yaml","S4.yaml","LB.yaml"]},
)
