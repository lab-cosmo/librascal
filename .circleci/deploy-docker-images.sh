#!/usr/bin/env bash
set -u

ALL_ENVS=("gcc-5" "gcc-9" "clang-4" "clang-9")
VERSION=$(cat VERSION)

if [[ "${VERSION}" == "" ]]
then
    echo "missing the VERSION file. Are you running this from .circleci?"
    exit 1
fi

for env in ${ALL_ENVS[@]}
do
    cd $env
    docker build -t rascal-ci-$env .
    docker tag rascal-ci-$env cosmoepfl/rascal-ci-$env:$VERSION
    cd -
done

for env in ${ALL_ENVS[@]}
do
    docker push cosmoepfl/rascal-ci-$env:$VERSION
done
