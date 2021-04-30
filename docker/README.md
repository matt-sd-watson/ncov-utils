```
function test() {
  docker build -t auspice . 
  # retrieve correct docker image id
  docker run -dit #image id
  docker ps #retrieve container id
  docker exec #container_id auspice view --datasetDir path_to_auspice_ncov/
}
```
