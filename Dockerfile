FROM eclipse-temurin:17-jre

# Install dependencies for X11 and IGV
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    xvfb \
    libxext6 \
    libxrender1 \
    libxtst6 \
    libxi6 \
    libxrandr2 \
    libgtk-3-0 \
    && rm -rf /var/lib/apt/lists/*

# Download and install IGV
WORKDIR /opt
# Using a specific version of IGV (2.16.2) known to work well
RUN wget https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_2.16.2.zip \
    && unzip IGV_2.16.2.zip \
    && mv IGV_2.16.2 igv \
    && rm IGV_2.16.2.zip

# Add IGV to path
ENV PATH="/opt/igv:${PATH}"

# Create a wrapper script to run IGV with xvfb
# Note: IGV 2.16+ uses a different launcher script usually, let's just point to igv.sh provided by them or invoke java directly if needed.
# The zip usually contains 'igv.sh'. Let's check if we can just use that.
# But we need xvfb-run wrapper.
RUN echo '#!/bin/bash\n\
xvfb-run --auto-servernum --server-num=1 --server-args="-screen 0 1024x768x24" java -Xmx4g --module-path="/opt/igv/lib" --module=org.igv/org.broad.igv.ui.Main "$@"' > /usr/local/bin/igv_headless.sh \
    && chmod +x /usr/local/bin/igv_headless.sh

ENTRYPOINT ["/bin/bash"]
