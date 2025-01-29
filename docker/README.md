# G2C Project Dockerfiles

Each subdirectory herein contains one Dockerfile corresponding to one Docker image used in this project.  

Example build command executed from this directory:
```
  docker build \
    -f g2c_pipeline/Dockerfile \
    --platform linux/x86_64 \
    --progress plain \
    --tag vanallenlab/g2c_pipeline:$TAG \
    ./
```
