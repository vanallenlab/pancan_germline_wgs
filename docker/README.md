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