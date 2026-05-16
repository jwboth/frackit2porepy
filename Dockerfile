FROM git.iws.uni-stuttgart.de:4567/tools/frackit/ubuntu18.04:latest

ENV DEBIAN_FRONTEND=noninteractive

# 1) Ensure Frackit is installed/updated in the image
RUN bash ./install_or_update.sh

# 2) Copy your application sources and config file into Frackit's appl tree
COPY appl/frackit2porepy/ frackit/appl/frackit2porepy/

# 3) Copy the appl/CMakeLists.txt that registers the application
COPY appl/CMakeLists.txt frackit/appl/CMakeLists.txt

# 4) Explicitly build the target to fail the image build early
RUN cd frackit/build && cmake .
RUN cd frackit/build/appl/frackit2porepy && make -j"$(nproc)" frackit2porepy

# 5) Copy the config file to a shared location that can be mounted from outside the container. This allows users to modify the config without rebuilding the image.
COPY shared/config.toml frackit/shared/config.toml

# Copy entrypoint
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

VOLUME ["/frackit/shared"]
ENTRYPOINT ["/entrypoint.sh"]
CMD []