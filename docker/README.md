The container to access auspice view on a server can be executed as follows, provided that the Dockerfile is in the current directory. The
user will need to replace the path to the auspice folder during docker run with the appropriate directory path:
 
```
  docker build -t auspice . 
  docker run -p 4000:4000 -v path_to_auspice_folder/:/app auspice auspice view --datasetDir .
```
