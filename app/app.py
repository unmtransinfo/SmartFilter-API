import yaml
from blueprints.version import register_routes
from config import Config
from flasgger import LazyJSONEncoder, Swagger
from flask import Flask
from flask_cors import CORS
from werkzeug.middleware.proxy_fix import ProxyFix

def _load_api_spec() -> dict:
    with open("swagger_template.yml", "r") as file:
        api_spec = yaml.safe_load(file)
    return api_spec

def create_app():
    app = Flask(__name__)
    app.wsgi_app = ProxyFix(app.wsgi_app, x_proto=1, x_host=1)
    app.json_encoder = LazyJSONEncoder

    app.config.from_object(Config)


    # load swagger template
    swagger_template = _load_api_spec()
    VERSION = swagger_template["info"]["version"]
    VERSION_URL_PREFIX = f"/api/v{VERSION}"

    # setup swagger config
    URL_PREFIX = app.config.get("URL_PREFIX", "")
    IN_PROD = app.config.get("FLASK_ENV", "") == "production"
    print("IN_PROD value:", IN_PROD)
    swagger_config = {
    "headers": [],
    "specs": [
        {
            "endpoint": "apispec_1",
            "route": f"/{URL_PREFIX}/apidocs/apispec_1.json",
        }
    ],
    "static_url_path": f"/{URL_PREFIX}/flasgger_static/",
    "swagger_ui": True,
    "specs_route": f"/{URL_PREFIX}/apidocs/",
    "ui_params": {
           "url_prefix": f"/{URL_PREFIX}/" if IN_PROD else "",
    },
    }
    swagger_template["swaggerUiPrefix"] = f"/{URL_PREFIX}/" if IN_PROD else ""
    # setup swagger and register routes
    swagger = Swagger(app, config=swagger_config, template=swagger_template)
    register_routes(app, IN_PROD, VERSION_URL_PREFIX)
    return app


app = create_app()
CORS(app)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=app.config.get("APP_PORT"), debug=False)

