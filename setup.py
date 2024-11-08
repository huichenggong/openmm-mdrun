from setuptools import setup, find_packages
import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

setup(
    name='Omdrun',
    version=get_version("Omdrun/__init__.py"),
    author='Chenggong Hui',
    author_email='chenggong.hui@mpinat.mpg.de',
    description='A wrapper on openmm for basic mdrun + GCNCMC',
    packages=find_packages(),
    install_requires=["openmm"],
    python_requires='>=3.7',
    entry_points={
        'console_scripts': [
            'openmm_mdrun=Omdrun.openmm_mdrun:main',
            'openmm_GCNCMC=Omdrun.openmm_GCNCMC:main',
        ],
    },
    classifiers=['Programming Language :: Python :: 3',],
)