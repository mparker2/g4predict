## g4predict

### Predict putative G Quadruplexes using an extension of the Quadparser method

<img src="./g4folding.gif" width="300">

Requirements:

* Linux (merge and filter_overlapping methods require calling sort behind the scenes)
* Python 3.5 (have not tested with Python 2.7 but it may work!)
* regex>=2016.3.2

regex module should be installed automatically by setuptools

TO INSTALL:

    tar xf g4predict.tar.gz
    cd g4predict
    python setup.py build
    python setup.py install
