#!/bin/bash
# source this file to get asked for your name and
# your email address, then the appropiate git env is setup

echo "Please enter your full name: "
read GIT_NAME
echo "Please enter email address: "
read GIT_EMAIL

export GIT_COMMITTER_NAME=$GIT_NAME
export GIT_COMMITTER_EMAIL=$GIT_EMAIL
export GIT_AUTHOR_NAME=$GIT_NAME
export GIT_AUTHOR_EMAIL=$GIT_EMAIL
