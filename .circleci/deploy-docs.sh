#!/usr/bin/env bash
set -xu

git config --global user.email "$GH_EMAIL"
git config --global user.name "$GH_NAME"

git clone $CIRCLE_REPOSITORY_URL --branch gh-pages gh-pages
rm -rf gh-pages/*

cp -a build/docs/html/. gh-pages/.

# disable CircleCI on gh-pages branch
mkdir -p gh-pages/.circleci
cat << EOF > gh-pages/.circleci/config.yml
version: 2

jobs:
  build:
    docker:
      # Very small docker image
      - image: alpine:3
    steps:
      - run: echo "done"

workflows:
  version: 2
  build_and_test:
    jobs:
      - build
EOF

cd gh-pages
git add -A
git commit -m "Automated deployment to GitHub Pages: ${CIRCLE_SHA1}" --allow-empty

cd -
