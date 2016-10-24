from setuptools import setup
import unittest


def test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='*_test.py')
    return test_suite

setup(
    name='g4predict',
    version='0.1',
    description=(
        'Predict putative G Quadruplexes using an extension'
        'of the Quadparser method'
    ),
    author='Matthew Parker',
    entry_points={
        'console_scripts': [
            'g4predict = g4funcs.g4predict:main'
        ]
    },
    packages=[
        'g4funcs',
    ],
    install_requires=[
        'regex>=2016.3.2'
    ],
    test_suite='setup.test_suite'
)
