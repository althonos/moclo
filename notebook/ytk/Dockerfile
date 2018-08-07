FROM python:alpine

ADD *.py *.ipynb nb/
ADD misc/* misc/

RUN apk add --no-cache --update $(cat misc/libs.txt) \
 && apk add --no-cache --update --virtual=.build-dependencies $(cat misc/build.txt) \
 && pip install --no-cache-dir -r misc/requirements.txt \
 && jupyter nbextension enable --py widgetsnbextension \
 && jupyter serverextension enable --py jupyterlab \
 && jupyter labextension install @jupyter-widgets/jupyterlab-manager qgrid \
 && /misc/glib-patch.sh \
 && apk del .build-dependencies

ENV LANG=C.UTF-8
EXPOSE 8888
CMD jupyter lab \
    --ip=* \
    --port=8888 \
    --no-browser \
    --NotebookApp.token='' \
    --notebook-dir=/nb \
    --allow-root \
    /nb/YTK\ Protocol.ipynb
