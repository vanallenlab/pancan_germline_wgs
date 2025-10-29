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

Note that the main G2C analysis image, `g2c_analysis`, has its environment pre-packaged in `g2c_analysis_base`. Thus, if you are trying to update `g2c_analysis`, you should usually first rebuild `g2c_analysis_base`.  

For example, to rebuild the `g2c_analysis` image to reflect the repo hash `26ce17a`, you would run the following on an arm64 (Apple M1/M2) machine:  

```
# 1. Declare the git hash as an environment variable for convenience
export g2c_hash=26ce17a

# 2. Build the base image, which installs all dependencies but results in a huge (>8GB) image
docker build \
  -f g2c_analysis_base/Dockerfile \
  --platform linux/x86_64 \
  --progress plain \
  --tag vanallenlab/g2c_analysis_base:$g2c_hash \
  --build-arg G2C_TAG="$g2c_hash" \
  ./

# 3. Build the analysis docker, which contains a compressed frozen version of the analysis environment from the base image
docker build \
  -f g2c_analysis/Dockerfile \
  --platform linux/x86_64 \
  --progress plain \
  --tag vanallenlab/g2c_analysis:$g2c_hash \
  --build-arg BASE_HASH="$g2c_hash" \
  --build-arg G2C_TAG="$g2c_hash" \
  ./

4. Push only the analysis docker image to Dockerhub for deployment
docker push vanallenlab/g2c_analysis:$g2c_hash
```