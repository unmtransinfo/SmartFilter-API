# models.py
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class PainsRun(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    timestamp = db.Column(db.DateTime, server_default=db.func.now())
    inputs  = db.Column(db.JSON, nullable=False)  
    passed  = db.Column(db.JSON, nullable=False)  
    failed  = db.Column(db.JSON, nullable=False)  
