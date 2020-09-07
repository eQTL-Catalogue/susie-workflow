### Build the Docker container
sudo docker build -t eqtlcatalogue/susie-finemapping:v20.08.1 .

### Push to DockerHub
docker push eqtlcatalogue/susie-finemapping:v20.08.1

### Build a local copy of the Singularity container
singularity build susie-finemapping.img docker://kauralasoo/susie-finemapping:latest
