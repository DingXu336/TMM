#!/usr/bin/env python3
"""
Setup configuration for TMM (Transfer Matrix Method) package.

This package provides Transfer Matrix Method calculations for optical properties
of layered structures, commonly used in physical chemistry and optics research.

Author: Ding Xu (Physical Chemistry Researcher)
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="tmm-optics",
    version="1.0.0",
    author="Ding Xu",
    author_email="dingxu@example.com",  # Replace with actual email
    description="Transfer Matrix Method for optical calculations in layered structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dingxu/tmm-optics",  # Replace with actual repository
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "black>=22.0",
            "flake8>=5.0",
            "mypy>=1.0",
            "sphinx>=5.0",
            "sphinx-rtd-theme>=1.0",
        ],
        "gui": [
            "tkinter",  # Usually included with Python
        ],
    },
    entry_points={
        "console_scripts": [
            "tmm-gui=tmm.gui.main:main",
            "tmm-convert=tmm.utils.converter:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
) 