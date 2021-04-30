The container to access auspice view on a server can be executed as follows, provided that the Dockerfile is in the current directory: 
```
  docker build -t auspice . 
  # retrieve correct docker image id
  docker run -dit #image id
  docker ps #retrieve container id
  docker exec #container_id auspice view --datasetDir path_to_auspice_ncov/
```
