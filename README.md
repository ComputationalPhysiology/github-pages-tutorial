[![testing Caesar cipher](https://github.com/finsberg/2021-12-06-github-actions/actions/workflows/main.yml/badge.svg)](https://github.com/finsberg/2021-12-06-github-actions/actions/workflows/main.yml)

[![Test laplace](https://github.com/finsberg/2021-12-06-github-actions/actions/workflows/test-laplace.yml/badge.svg)](https://github.com/finsberg/2021-12-06-github-actions/actions/workflows/test-laplace.yml)

# GitHub actions and Github Pages

An example repo for how to set up github actions and GitHub pages. 
Read more about github actions at https://docs.github.com/en/actions


## GitHub actions

### Example: Caesar cipher
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

#### Testing multiple python versions

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


#### Test different python versions on different operating systems


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

### Test FEniCS code


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
```

## GitHub pages

Github pages allows you to turn a Github repo into a website.
This is a good way to host documentation for your github project. 
You can read the full documentation here: https://docs.github.com/en/pages

Other popular alternatives are https://readthedocs.org and https://www.netlify.com.

The website will be hosted at https://username.github.io/myrepo
where `username` is your github username or organization name and `myrepo` it the name of the `repo`.

For example I have create a repo called `ap_features` at https://github.com/ComputationalPhysiology/ap_features and its documentation is hosted on https://computationalphysiology.github.io/ap_features

### Building a website for a repo

If you want to do create a website for your repo you need to do the following
1. Create a clean branch called `gh-pages`
   ```
   git checkout --orphan gh-pages # Creates a clean branch
   git rm -rf .  # Delete everything in the folder (or do 'git rm -rf --dry-run .' so see which files that will be deleted)
   ```

2. Add the html for your website, commit it to the branch and push it
3. Go to Setting->Pages, select `gh-pages` as source

#### Example

Try creating a simple html file called `index.html` with the following content
```html
<h1>Hello from Github pages</h1>
```
Now commit it, push it and update the settings.


#### Deleting the website

To delete the website you simply delete the `gh-pages` branch

### Using a static site generator


#### Jekyll

A popular static site generator that works very well for GitHub pages is Jekyll. You can read more about Jekyll here: https://docs.github.com/en/pages/setting-up-a-github-pages-site-with-jekyll


#### Sphinx
Sphinx is a popular static site generator that is used to generate documentation for python projects, see e.g https://github.com/finsberg/sphinx-tutorial for a sphinx-tutorial that was given on an earlier tools meetup.

- Go back to the `master` branch
  ```
  git checkout master
  ```
- Install `sphinx`
  ```
  python -m pip install sphinx
  ```
- Create a directory called `docs`, go into it and run Sphinx quick-start
  ```
  mkdir docs
  cd docs
  sphinx-quickstart
  ```
- Hit `make html` to make the html
- Take a look at the generated website by opening `_build/html/index.html`
- Open `index.rst` and try to add some content and build the docs again.
- Try to change the `html_theme` from `alabaster` to `classic` (more themes can be found here: https://sphinx-themes.org)
- Open `conf.py` and update `extensions` to be the following
    ```python
    extensions = [
            "sphinx.ext.autodoc",
            "sphinx.ext.mathjax",
            "sphinx.ext.viewcode",
            "sphinx.ext.napoleon",
            "sphinx.ext.githubpages",
        ]
    ```
    Also add the following line on the top
    ```python
    from pathlib import Path
    import sys
    here = Path(__file__).parent
    sys.path.insert(0, here.parent.as_posix())
    ```
- Run `sphinx-apidoc` doc to generate documentation from your code
    ```
    # From /docs folder
    sphinx-apidoc -o . ..
    ```
    Now try to rebuild the documentation
    Note that if you are using e.g `dolfin` and you don't have `dolfin` installed you can still build the documentation, but you need to mock the import first. To to this add the following to your `conf.py`
    ```python
    autodoc_mock_imports = ["dolfin"]
    ```
    Also note that it you change your docstring, then you need to also run `sphinx-apidoc` again so it might be a good idea to automate this process. 
    Finally you should add the files to your `index.rst`, e.g if you do
    ```
    $ sphinx-apidoc -o . ..
    Creating file ./caesar_cipher.rst.
    Creating file ./test_caesar_cipher.rst.
    Creating file ./modules.rst.
    ```
    then add the following to `index.rst`
    ```
    .. toctree::
       :maxdepth: 2
       :caption: Contents:
    
       caesar_cipher
       test_caesar_cipher
       modules
    ```
    (note the indentation here)


If you want to write your documentation in Markdown in stead of reStructuredText, you can install `myst-parser` and add it as an extension to the `conf.py`, see https://www.sphinx-doc.org/en/master/usage/markdown.html

### Automating the build using github actions

Currently all the html is now in `docs/_build/html` so we need to somehow move this to the `gh-pages` branch to the root folder. Luckily for us there is an existing Github action for us that we can use

Delete the files generated from `sphinx-apidoc` as well as the `dock/_build` folder and add the `docs` folder as to your `master` branch (it might be a good idea to add these files and folders to your `.gitignore` file).

Next, create a file called `.github/workflows/docs.yml` with the following content

```
name: github pages

on:
  push:
    branches:
      - master

jobs:
  deploy:
    runs-on: ubuntu-18.04
    steps:
      - uses: actions/checkout@v2

      - name: Setup Python
        uses: actions/setup-python@v2
        with:
          python-version: "3.8"

      - name: Upgrade pip
        run: |
          python3 -m pip install --upgrade pip
      - name: Install dependencies
        run: |
          python -m pip install sphinx
      - name: Build docs
        run: |
          sphinx-apidoc -o docs .
          cd docs && make html

      - name: Deploy
        uses: peaceiris/actions-gh-pages@v3
        with:
          github_token: ${{ secrets.GITHUB_TOKEN }}
          publish_dir: ./docs/_build/html
```

### Building a website for your user or organization

It is possible to also build a website for your Github user or organization. For example, I have also created a website at https://computationalphysiology.github.io which is indented to describe who we are (but it is work in progress).

In this case you need to create a repository called `username.github.io`, see e.g https://github.com/ComputationalPhysiology/computationalphysiology.github.io. Note that in this case, the html must be at the `master` or `main` branch.