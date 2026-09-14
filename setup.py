from setuptools import setup, find_packages

setup(
    name="seventdt",
    version="0.2.0",
    description="Empirical exposure time calculator for the 7DT telescope",
    author="GWUniverse @ Seoul National University",
    python_requires=">=3.10",
    packages=find_packages(),
    package_data={"seventdt": ["depth_model_*.json"]},
    install_requires=[
        "numpy",
        "astropy",
        "astroplan",
        "matplotlib",
        "tqdm",
    ],
    extras_require={
        "fit": ["scikit-learn"],
    },
)
