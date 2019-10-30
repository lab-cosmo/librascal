# docker image used for Circle CI

In order to cut down on the run time of circle CI, custom docker images are
used, one for for each compiler. It is based on Ubuntu 18.04, with rascal
specific dependencies installed.  See the `<compiler>/Dockerfile` for more
information.

## Updating the images

1. update the VERSION file
2. make your changes to one of the `Dockerfile`s
3. build and upload the images to dockerhub using the `deploy-docker-images.sh`
   script. You will need to login to the dockerhub
   [cosmoepfl](https://cloud.docker.com/u/cosmoepfl/repository/list)
   organization.
4. update the config.yml file to reference the new VERSION

You can then test you changes locally using `circleci local execute --job
<job_name>`. Make sure to commit your changes before, since the local repository
is cloned in the `checkout` step.
