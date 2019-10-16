# docker image used for Circle CI

In order to cut down on the run time of circle CI, a custom docker image is
used. It is based on Ubuntu 18.04, with rascal specific dependencies installed.
See the [Dockerfile] for more information.

## Updating the image

1. make your changes to the Dockerfile
2. build the image locally
```bash
cd .circleci
docker build -t rascal-ci .
```
3. tag the image (replace <VERSION> with the right value) and upload to docker hub
```bash
docker tag rascal-ci cosmoepfl/rascal-ci:<VERSION>
docker push cosmoepfl/rascal-ci:<VERSION>
```
