#!/bin/bash

set -e

# the secret password is not set for pull requests
# and we also don't want the doc to be updated for that.. so that fits.

if [[ $TRAVIS_PULL_REQUEST != 'false' ]]; then exit; fi
if [[ $TRAVIS_BRANCH != 'master' ]]; then exit; fi

openssl aes-256-cbc -K $encrypted_b7ce407a834b_key -iv $encrypted_b7ce407a834b_iv -in ci/travis_rsa.enc -out ci/travis_rsa -d
chmod 0600 ci/travis_rsa
cp ci/travis_rsa ~/.ssh/id_rsa 


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
git reset --hard 7742183a2e40171ab3c3d4aaaf1f9cedde566b37
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
git push -f origin gh-pages
cd -
