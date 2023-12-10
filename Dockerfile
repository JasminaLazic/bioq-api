
FROM rstudio/plumber

# install the linux libraries needed for plumber
RUN apt-get update -qq && apt-get install -y \
  libssl-dev \
  libcurl4-gnutls-dev \
  libgdal-dev \
  libproj-dev \
  libgeos-dev

# create the application folder
RUN mkdir -p ~/application

# copy everything from the current directory into the container
COPY "/" "application/"
WORKDIR "application/data/" 

# open port 80 to traffic
EXPOSE 80

# install plumber
RUN R -e "install.packages('terra', dependencies = TRUE, repos = 'https://cloud.r-project.org/')"
RUN R -e "install.packages('plumber', dependencies = TRUE, repos = 'https://cloud.r-project.org/')"

# OPTIONAL: Install any additional R packages you may need for your model crate to run
RUN R -e "install.packages(c( 'data.table', 'plyr', 'dplyr', 'reshape2', 'rredlist', 'vegan', 'jsonlite', 'FD', 'igraph', 'bipartite', 'tidyr', 'ggplot2', 'terra', 'raster', 'fundiversity', 'randomForest', 'DT'), repos='https://cran.rstudio.com/')"


# when the container starts, start the main.R script
CMD ["BioQApi.R"]