# This script is used to serve the static files from the 'dist' directory for installation of python wheels.
# 1. Install deps with `pip install ".[dev]"`,
# 2. Install build tool with `pip install build`
# 3. Run `python -m build`
# You should see something alongside:
# ```Successfully built mat3ra_api_examples-dev9+g7c6e8d9.tar.gz and
# mat3ra_api_examples-dev9+g7c6e8d9-py3-none-any.whl
# ```
# 4. Copy the wheel file name
# 5. Run the server with `python wheel_server.py`
# The server will be available at http://localhost:8080,
# 6. use micropip in pyodide to install wheel from that URL `<server_url>/<wheel_file_name>`, i.e:
# await micropip.install("http://localhost:8080/mat3ra_api_examples-dev9+g7c6e8d9-py3-none-any.whl", deps=False)

import os
import socket
from http.server import HTTPServer, SimpleHTTPRequestHandler


class CORSHTTPRequestHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "X-Requested-With, Content-Type")
        return super(CORSHTTPRequestHandler, self).end_headers()


def check_port(host, port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex((host, port)) == 0


if __name__ == "__main__":
    port = 8080
    bind_addr = "localhost"
    directory = "./dist"  # make sure this is the correct relative path to your 'dist' directory

    os.chdir(directory)  # Change the current working directory to the specified 'directory'

    while check_port(bind_addr, port):
        print(f"Port {port} is already in use. Trying with port {port + 1}.")
        port += 1

    httpd = HTTPServer((bind_addr, port), CORSHTTPRequestHandler)
    print(f"Serving at http://{bind_addr}:{port}")
    httpd.serve_forever()
