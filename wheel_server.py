# This script is used to serve the static files from the 'dist' directory for installation of python wheels.
# Install deps with `pip install ".[dev]"`, build tool `pip install build` and run `python -m build`
# Run the server with `python wheel_server.py`
# The server will be available at http://localhost:8080, use micropip in pyodide to install wheel from that URL, i.e:
# await micropip.install("http://localhost:8080/mat3ra_api_examples-2024.3.30.post2-py3-none-any.whl", deps=False)

import os
from http.server import HTTPServer, SimpleHTTPRequestHandler


class CORSHTTPRequestHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "X-Requested-With, Content-Type")
        return super(CORSHTTPRequestHandler, self).end_headers()


if __name__ == "__main__":
    port = 8080
    bind_addr = "localhost"
    directory = "./dist"  # make sure this is the correct relative path to your 'dist' directory

    os.chdir(directory)  # Change the current working directory to the specified 'directory'

    httpd = HTTPServer((bind_addr, port), CORSHTTPRequestHandler)
    print(f"Serving at http://{bind_addr}:{port}")
    httpd.serve_forever()
