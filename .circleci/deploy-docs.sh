#!/usr/bin/env bash
set -xu

git config --global user.email "$GH_EMAIL"
git config --global user.name "$GH_NAME"

git clone $CIRCLE_REPOSITORY_URL --branch gh-pages gh-pages
rm -rf gh-pages/*

cp -a build/docs/html/. gh-pages/.

mkdir -p gh-pages/.circleci && cp -a .circleci/. gh-pages/.circleci/.

cd gh-pages
git add -A
git commit -m "Automated deployment to GitHub Pages: ${CIRCLE_SHA1}" --allow-empty

cd -
