# How to create a Bioconda release

1. Create a `meta.yaml` file similar to [this example](meta.yaml).
2. Open a PR on [Bioconda recipes](https://github.com/bioconda/bioconda-recipes)
3. Follow the [contribution manual](https://bioconda.github.io/contributor/index.html) to replicate the [contribution workflow](https://bioconda.github.io/contributor/workflow.html)

See [this PR](https://github.com/bioconda/bioconda-recipes/pull/20701) for an example of how it's done.

# How to be notified of PR requests

Depending on how your Bioconda release was created, whenever you create a new GitHub release, the Bioconda Bot will pick up the changes and automatically create a pull request (PR).

To properly be notified of such PR's and to be able to contribute to them, you will need to:

1. Join the [Bioconda gitter lobby](https://gitter.im/bioconda/Lobby) and ask nicely to become a member.
2. Make sure that the yaml file lists you as a recipe maintainer.

# How to resolve a Bioconda PR issues.

Eg: [This PR](https://github.com/bioconda/bioconda-recipes/pull/21112) failed. Let's see what the issue might be.

1. Check the [Circle-CI](https://github.com/bioconda/bioconda-recipes/pull/21112/checks?check_run_id=534673776) results to see where the build fails
2. Click on the appropriate link [like this one](https://circleci.com/gh/bioconda/bioconda-recipes/100487?utm_campaign=vcs-integration-link&utm_medium=referral&utm_source=github-checks-link) and notice that code changes now require python>=3.7
3. Clone the appropriate branch, make changes, and then push:
```bash
$ git clone https://github.com/bioconda/bioconda-recipes.git -b bump/cmash
$ cd bioconda-recipes/recipes/cmash
<edit the yaml file>
$ git commit -a  # Make sure to add what changes you made and the PR issue number to the commit message
$ git push

```

# How to test build locally to make sure that you don't run into a Bioconda PR issue.

This is still a work in progress, have not been able to get anything to install properly, but this is how Bioconda is doing it on their back-end:
```
#!/bin/bash -eo pipefail
bioconda-utils build recipes config.yml \
  --docker --mulled-test \
  --git-range master HEAD
```


## Make sure you have Docker installed

If you don't have Docker installed, do the following:
```
sudo apt-get update
sudo apt install docker.io
sudo systemctl start docker
sudo systemctl enable docker
```

## Create bioconda conda environment

1. Add the bioconda channels to your conda installation
```bash
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```
2. Install bioconda utils
```bash
$ conda create -y -n bioconda python=3.6  # seems to fail with python>=3.7
$ conda activate bioconda
$ conda install -y bioconda-utils
```
## Install Bioconda CircleCI

This is the continuous integration checker that the Bioconda bot is using to check if your recipe builds.

1. Install the CircleCI docker image: `sudo docker pull bioconda/bioconda-utils-build-env:latest`
2. `sudo docker run -t --name test bioconda/bioconda-utils-build-env`
2. `cd` to where you made your changes. eg. `cd bioconda-recipes/recipes/cmash` 
 
