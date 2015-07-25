#!/bin/bash -e

# Settings
REPO_PATH=git@github.com:A2-Collaboration-dev/ant.git
HTML_PATH=doxygen/html
COMMIT_USER="Documentation Builder"
COMMIT_EMAIL="neiser@kph.uni-mainz.de"
CHANGESET=$(git rev-parse --verify HEAD)

# Get a clean version of the HTML documentation repo.
rm -rf ${HTML_PATH}
mkdir -p ${HTML_PATH}
git clone -b gh-pages "${REPO_PATH}" --single-branch ${HTML_PATH}

# rm all the files through git to prevent stale files.
cd ${HTML_PATH}
git rm -rf .
cd -

# Generate the HTML documentation.
cd build
make doxygen
cd -

# Create and commit the documentation repo.
cd ${HTML_PATH}
rm -f *.map *.md5
git add .
git config user.name "${COMMIT_USER}"
git config user.email "${COMMIT_EMAIL}"
git commit -m "Automated documentation build for changeset ${CHANGESET}."
git push origin gh-pages
cd -
