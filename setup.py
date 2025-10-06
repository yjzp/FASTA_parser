"""
Setup script для установки пакета fasta-parser
"""

from setuptools import setup, find_packages
import os

# Читаем README файл
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Читаем версию из __init__.py
version = {}
with open("fasta_parser/__init__.py") as fp:
    exec(fp.read(), version)

setup(
    name="fasta-parser",
    version=version['__version__'],
    author=version['__author__'],
    author_email="biology@university.edu",
    description="Python библиотека для работы с биологическими последовательностями в формате FASTA",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/username/fasta-parser",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        # Никаких внешних зависимостей для основного функционала
    ],
    extras_require={
        "dev": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
            "pytest>=6.0.0",
            "pytest-cov>=2.0.0",
            "black>=21.0.0",
            "flake8>=3.8.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "fasta-demo=examples.demo:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)