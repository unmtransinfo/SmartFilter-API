import os
from os import environ

FLASK_ENV = environ.get("FLASK_ENV")
APP_NAME = environ.get("APP_NAME")
APP_PORT = environ.get("APP_PORT")
APP_URL = environ.get("APP_URL") or "localhost"
URL_PREFIX = environ.get("URL_PREFIX") or ""

class Config:
    SQLALCHEMY_DATABASE_URI = os.getenv("DATABASE_URL", "sqlite:///pains.db")
    SQLALCHEMY_TRACK_MODIFICATIONS = False
