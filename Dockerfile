# Take inspiration from
# https://github.com/jrnold/docker-stan/blob/master/rstan/Dockerfile

FROM rocker/tidyverse:4.1.3

RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils ed libnlopt-dev \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/

# Global site-wide config -- neeeded for building packages
RUN mkdir -p $HOME/.R/ \
    && echo "CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -flto -ffat-lto-objects  -Wno-unused-local-typedefs \n" >> $HOME/.R/Makevars \
    && echo "CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -flto -ffat-lto-objects  -Wno-unused-local-typedefs -Wno-ignored-attributes -Wno-deprecated-declarations\n" >> $HOME/.R/Makevars

# Install rstan
# Instead, we can install rstan by restoring package library with renv (remove dependency on EczemaPredPOEM first)
RUN install2.r --error --deps TRUE \
    rstan \
    rstantools \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# Create a user variable
ENV USER=rstudio

# Create project directory and set it as working directory
WORKDIR /home/$USER/EczemaTreat

# Install R packages to local library using renv
COPY [".Rprofile", "renv.lock", "./"]
COPY renv/activate.R ./renv/
RUN chown -R rstudio . \
 && sudo -u rstudio R -e 'renv::restore(confirm = FALSE, exclude = "TanakaData")'
