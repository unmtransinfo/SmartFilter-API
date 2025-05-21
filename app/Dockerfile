FROM python:3.12

WORKDIR /app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

RUN mkdir -p flasgger_static && \
    cp -r /usr/local/lib/python3.12/site-packages/flasgger/ui3/static/* flasgger_static

ENV APP_PORT=8000
CMD echo "RUNTIME: APP_PORT is set to ${APP_PORT}" && gunicorn --bind "0.0.0.0:${APP_PORT}" --reload app:app