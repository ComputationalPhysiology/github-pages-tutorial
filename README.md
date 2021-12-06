# GitHub actions

An example repo for how to set up github actions. 
Read more about github actions at https://docs.github.com/en/actions


## Example: Caesar cipher
In this example we only need to install pytest

```
name: testing Caesar cipher   # Name of the Ci

on: [push]  # When to run the CI

jobs:
  test:  # Name of job
    name: Test on ubuntu with python 3.9
    runs-on: ubuntu-latest  # OS to run on
    
    steps:
      - uses: actions/checkout@v2  # Pre-made action to check out repo
      - name: Set up Python
        uses: actions/setup-python@v2   # Pre-made action to set up python 
        with:
          python-version: 3.9
      - name: Install dependencies
        run: |
          python -m pip install pip --upgrade
          python -m pip install pytest pytest-cov
      - name: Test with pytest
        run: |
          python -m pytest test_caesar_cipher.py --cov=caesar_cipher --cov-report term-missing
```

### Testing multiple python versions

```
name: testing Caesar cipher   # Name of the Ci

on: [push] 

jobs:
  test:
    name: Test on ubuntu with python ${{ matrix.python-version }}
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.7, 3.8, 3.9]

    steps:
      - uses: actions/checkout@v2
      - name: Set up Python ${{ matrix.python-version }}
        uses: actions/setup-python@v2
        with:
          python-version: ${{ matrix.python-version }}
      - name: Install dependencies
        run: |
          python -m pip install pip --upgrade
          python -m pip install pytest pytest-cov
      - name: Test with pytest
        run: |
          python -m pytest test_caesar_cipher.py --cov=caesar_cipher --cov-report term-missing
```


### Test different python versions on different operating systems


```
name: testing Caesar cipher   # Name of the Ci

on: [push]

jobs:
  test:
    name: Test on ${{ matrix.os }} with python ${{ matrix.python-version }}
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
        python-version: [3.7, 3.8, 3.9]

    steps:
      - uses: actions/checkout@v2
      - name: Set up Python ${{ matrix.python-version }}
        uses: actions/setup-python@v2
        with:
          python-version: ${{ matrix.python-version }}
      - name: Install dependencies
        run: |
          python -m pip install pip --upgrade
          python -m pip install pytest pytest-cov
      - name: Test with pytest
        run: |
          python -m pytest test_caesar_cipher.py --cov=caesar_cipher --cov-report term-missing
```

## Test FEniCS code


```
name: Test laplace

on: [push]

jobs:
  test:
    name: Test
    runs-on: ubuntu-latest
    container:
      image: quay.io/fenicsproject/stable:latest

    steps:
      - uses: actions/checkout@v2

      - name: Install dependencies
        run: |
          python -m pip install pip --upgrade
          python -m pip install pytest pytest-cov
      - name: Test with pytest
        run: |
          python -m pytest test_laplace.py --cov=laplace --cov-report term-missing
