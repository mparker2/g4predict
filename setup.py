from setuptools import setup

setup(
    name='g4predict',
    version='0.1',
    description=(
        'Predict putative G Quadruplexes using an extension'
        'of the Quadparser method'
    ),
    author='Matthew Parker',
    scripts=[
        'scripts/g4predict.py'
    ],
    packages=[
        'g4funcs',
    ],
    install_requires=[
        'regex>=2016.3.2'
    ]
)
