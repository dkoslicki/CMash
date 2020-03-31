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

This will allow you to add the tag `please review & merge` to request a review __after__ all tests are passing as indicated by the BiocondaBot.

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

Then wait for the BiocondaBot to run the continuous integration tests, and see what else might pop up.

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
 

# As close as I've gotten it to working
https://bioconda.github.io/recipes/bioconda-utils/README.html

https://github.com/bioconda/bioconda-utils/blob/master/Dockerfile

https://quay.io/repository/biocontainers/bioconda-utils?tab=tags
```
# or 0.16.14--py_0
# or sudo docker run -it bioconda/bioconda-utils-build-env:latest /bin/bash  # has yum
sudo docker pull quay.io/biocontainers/bioconda-utils:0.16.14--py_0
sudo docker run -it quay.io/biocontainers/bioconda-utils:0.16.14--py_0 /bin/bash
#conda update -n base conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y conda-build
git clone https://github.com/bioconda/bioconda-recipes.git
cd bioconda-recipes/recipes
conda-build cmash  # just check if it builds
git checkout -b local_test  # make a local branch to make changes to your meta.yaml file
vi /bioconda-recipes/recipes/cmash/meta.yaml  # make your edits
git config --global user.email "dmkoslicki@gmail.com"
git config --global user.name "dkoslicki"
git commit -a  # make all your changes to the local branch DO NOT PUSH
cd /bioconda-recipes
bioconda-utils build recipes config.yml --mulled-test --git-range master HEAD --force  # can't find the package manager on this docker image, so foregoe the --docker option that's actually being used.

vi /bioconda-recipes/recipes/cmash/meta.yaml
```

# Message sent to Bioconda team

tldr: How does one set up a local testing environment for Bioconda tests run by the Bioconda Bot?

tldr2: I want my students to be able to test locally so they don't need to push/make PR's to bioconda-recipes and wait for the Bioconda Bot to run the tests, then manually dig through the errors on GitHub.

tldr3: happy to open an issue instead of trying to discuss this here. Also happy to update [the docs](https://bioconda.github.io/contributor/building-locally.html) once a solution is found.

I'm trying to set up a local testing environment so I can check my recipes pass all the Circle-CI tests before bumping versions/making changes/etc to the Bioconda repo. The documentation [here](https://bioconda.github.io/contributor/building-locally.html), is unfortunately not reproducible. Examples:

1. Circle CI not in docker container
```bash
$ docker pull bioconda/bioconda-utils-build-env:latest
$ sudo docker run -it bioconda/bioconda-utils-build-env:latest /bin/bash  # missing this command and entry point in the above linked documentation
[root@xyz /]# circleci build
bash: circleci: command not found
```
Found a hidden folder `.circleci` in `bioconda-recipes`, but the `setup.sh` fails due to a missing `common.sh`. It [appears](https://circleci.com/gh/bioconda/bioconda-recipes/100722#config/containers/0) to come from some curl-ing of raw github user content, but even after `curl -s https://raw.githubusercontent.com/bioconda/bioconda-common/master/common.sh > .circleci/common.sh`, the circleci `setup.sh` fails.

2. Bootstrap method fails
```bash
$ sudo docker run -it bioconda/bioconda-utils-build-env:latest /bin/bash
[root@xyz /]# ./bootstrap.py /tmp/miniconda
ERROR conda.core.link:_execute(502): An error occurred while installing package 'conda-forge::aiofiles-0.4.0-py_1001'.
FileNotFoundError(2, "No such file or directory: '/tmp/miniconda/miniconda/bin/python3.7'")
```
`bootstrap.py` appears to be installing and using python3.6.

Tried using `bioconda-utils`  locally as [indicated in the readme](https://bioconda.github.io/recipes/bioconda-utils/README.html). But `conda install bioconda-utils` fails to run properly when installed locally:
```bash
$ conda install bioconda-utils  # using conda 4.8.3 on Ubuntu 18.04, installs without error
$ bioconda-utils -h
<snip>
ModuleNotFoundError: No module named 'gitdb.utils.compat'
```
Maybe related to [this other issue](https://github.com/gitpython-developers/GitPython/issues/983)? But no success when trying to downgrade gitdb and/or gitdb2 (earlier version not on conda, still no success using pip to pin the versions, same `ModuleNotFoundError`.

So going back to the docker container:

1. Using `bioconda-utils-build-env` with docker still fails
```bash
$ sudo docker run -it bioconda/bioconda-utils-build-env:latest /bin/bash
[root@xyz /]# clone, make branch, do trivial edit to package meta.yaml that I know passes all tests, commit them, etc.
[root@xyz /]# bioconda-utils build recipes config.yml --docker --mulled-test --git-range master HEAD --force
<snip>
No such file or directory: 'docker': 'docker'  # have fun trying to install docker on CentOS 6.10 due to docker requiring kernel>=3.10
```

2. Using `bioconda-utils-build-env` _without_ docker still fails
```bash
$ sudo docker run -it bioconda/bioconda-utils-build-env:latest /bin/bash
[root@xyz /]# clone, make local branch, do trivial edit to package meta.yaml that I know passes all tests, builds with conda-built, etc. and then commit changes
[root@xyz /]# bioconda-utils build recipes config.yml --mulled-test --git-range master HEAD --force
<snip>
(ERR) /opt/conda/conda-bld/cmash_1585274295882/_test_env_placehold_placehold_placehold_<snip>_placehold_placehold_pl/bin/MakeStreamingDNADatabase.py
BIOCONDA ERROR COMMAND FAILED (exited with 1)
```
Ok, so this appears to be related to issue [#21121](https://github.com/bioconda/bioconda-recipes/issues/21121), so let's pull the most up-to-date docker from [quay](https://quay.io/repository/biocontainers/bioconda-utils?tab=tags) and try again:

```bash
$ sudo docker run -it quay.io/biocontainers/bioconda-utils:0.16.14--py_0
# clone, make branch, etc.
[root@xyz /]# bioconda-utils build recipes config.yml --docker --mulled-test --git-range master HEAD --force
FileNotFoundError: [Errno 2] No such file or directory: '/usr/local/conda-bld/conda_build_config_0_-e_conda_build_config.yaml'
#ok, so now without docker (since this docker image appears to have no package manager)
[root@xyz /]# bioconda-utils build recipes config.yml --mulled-test --git-range master HEAD --force
# BIOCONDA ERROR
```
In more detail, the build worked, but the `mulled-test` failed due to:
```
/usr/local/bin/involucro -v=2 -f /usr/local/lib/python3.7/site-packages/galaxy/tools/deps/mulled/invfile.lua -set CHANNELS='conda-forge,file:///usr/local/conda-bld,bioconda,defaults' -set TARGETS='cmash=0.5.0=py_0' -set REPO='quay.io/biocontainers/cmash:0.5.0--py_0' -set BINDS='build/dist:/usr/local/,/usr/local/conda-bld:/usr/local/conda-bld' -set CONDA_IMAGE='quay.io/dpryan79/mulled_container:latest' -set TEST='MakeStreamingDNADatabase.py -h && MakeStreamingPrefilter.py -h && StreamingQueryDNADatabase.py --help && python -c "import CMash"' build-and-test
```
It failed since this appears to be calling a docker image `quay.io/dpryan79/mulled_container:latest` by @dpryan79, but as previously mentioned _there is no docker in this container and no way to install it_.

Let's step back one and then try an earlier version `0.16.13--py_0`:
```bash
$ sudo docker run -it quay.io/biocontainers/bioconda-utils:0.16.13--py_0
# clone, make branch, etc.
[root@xyz /]# bioconda-utils build recipes config.yml --docker --mulled-test --git-range master HEAD --force
FileNotFoundError: [Errno 2] No such file or directory: '/usr/local/conda-bld/conda_build_config_0_-e_conda_build_config.yaml'
#so once more without docker
[root@xyz /]# bioconda-utils build recipes config.yml --mulled-test --git-range master HEAD --force
# same placehold_placehold_..._placehold_placehold error
```

Soooo... exactly how are people managing to do local tests?


