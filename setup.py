"""Setup script for symmer-magic package."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="symmer-magic",
    version="0.1.5",
    author="Sam Alterman",
    author_email="samalterman@gmail.com",
    description="SRE calculations for the Symmer quantum chemistry package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https:/github.com/samalterman/symmer-magic",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "openfermion>=1.5.0",
        "symmer>=0.0.10",
        "sympy>=1.10",
        'joblib>=1.3.2'
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=3.0",
            "black>=22.0",
            "flake8>=4.0",
            "mypy>=0.950",
        ],
        "notebooks": [
            "jupyter>=1.0",
            "matplotlib>=3.5.0",
            "ipykernel>=6.0",
        ],
    }
)