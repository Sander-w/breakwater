# use python package as base
FROM python:3-slim
ARG GITHUB_TOKEN
ARG PYPI_URL=https://${GITHUB_TOKEN}@pypi.data.vanoord.com/

# upgrade packages
RUN apt update && apt install -y python3-tk git procps
RUN pip install --upgrade pip

# install python requirements
WORKDIR /breakwater
ADD ./requirements.txt /breakwater/requirements.txt
RUN pip install -r requirements.txt --extra-index-url https://$GITHUB_TOKEN@pypi.data.vanoord.com/

# install python package
ADD . /breakwater
RUN pip install -e . --extra-index-url https://$GITHUB_TOKEN@pypi.data.vanoord.com/

EXPOSE 8888

CMD ["sh", "-c", "tail -f /dev/null"]