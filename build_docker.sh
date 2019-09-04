### Build the Docker container
sudo docker build -t kauralasoo/susie-finemapping .

### Push to DockerHub
docker push kauralasoo/susie-finemapping

### Build a local copy of the Singularity container
singularity build susie-finemapping.img docker://kauralasoo/susie-finemapping:latest
