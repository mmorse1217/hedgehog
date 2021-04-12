image-build: Dockerfile
	docker build -f Dockerfile -t hedgehog .
container-create:
	echo "making container...";
	docker create -it -v `pwd`:/hedgehog --name hedgehog-dev hedgehog
container-start:
	echo "starting container...";
	docker start hedgehog-dev
container-stop:
	echo "stopping container...";
	docker stop hedgehog-dev
container-exec:
	echo "exec-ing into container...";
	docker exec -it hedgehog-dev bash
