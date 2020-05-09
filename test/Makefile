# Makefile to build and run docker
# Inpsired by: https://gist.github.com/mpneuried/0594963ad38e68917ef189b4e6a269db
SHELL=/bin/bash
APP_NAME=fv3atm_test

# HELP
# This will output the help for each task
# thanks to https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help
help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

.DEFAULT_GOAL := help

# DOCKER TASKS

.PHONY: build
build: ## Build the container
	docker build -t $(APP_NAME) .

.PHONY: build-nc
build-nc: ## Build the container without caching
	docker build --no-cache -t $(APP_NAME) .

.PHONY: run
run: ## Run container (no mounts)
	docker run -it --rm --mount type=bind,source=`pwd`/work,target=/work --name="$(APP_NAME)" $(APP_NAME)

.PHONY: enter
enter:
	docker exec -it $(APP_NAME) /bin/bash

.PHONY: stop
stop: ## Stop and remove a running container
	-docker stop $(APP_NAME); docker rm $(APP_NAME)

.PHONY: clean
clean: ## remove rundir
	-/bin/rm -rf ./work/rundir

.PHONY: space
space: ## remove all containers/images/volumes not currently in use (DANGER ZONE!!!)
	-docker system prune -a

