FROM continuumio/miniconda3

RUN apt update
RUN apt-get install -y xvfb
RUN apt-get install -qqy mesa-utils libgl1-mesa-glx  libglib2.0-0
RUN apt-get install -qqy libfontconfig1 libxrender1 libdbus-1-3 libxkbcommon-x11-0 libxi6
RUN apt-get install -qqy libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0
RUN apt-get install -qqy libxcb-xinerama0 libxcb-xinput0 libxcb-xfixes0 libxcb-shape0

RUN pip install gctree
RUN conda install -c bioconda phylip
COPY . /gctree
WORKDIR /gctree
CMD . tests/test.sh
