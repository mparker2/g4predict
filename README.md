## g4predict

### Predict putative G Quadruplexes using an extension of the Quadparser method

Requirements:

* Linux (program calls sort behind the scenes)
* Python 3.5 (have not tested with Python 2.7 but it probably won't work!)
* regex>=2016.3.2

these should be installed automatically by setuptools

TO INSTALL:

    tar xf g4predict.tar.gz
    cd g4predict
    python setup.py build
    python setup.py install

tested on:

Scientific Linux release 6.8 (Carbon)
Ubuntu 16.04.1 LTS (Xenial Xerus)

(both using conda)