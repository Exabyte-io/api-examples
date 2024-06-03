import argparse
import glob
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


def inform_user(port):
    whl_files = glob.glob("*.whl")
    file = whl_files[0] if whl_files else None
    url_str = f"http://localhost:{port}/{file}"
    print("Copy URL to use in notebook or `config.yml`: ", url_str)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Python wheel server.")
    parser.add_argument("--port", type=int, default=8080, help="Port to run the server on.")
    parser.add_argument("--dir", type=str, default="./dist", help="Directory to serve.")
    args = parser.parse_args()

    port = args.port
    bind_addr = "localhost"
    directory = args.dir  # Change this line

    os.chdir(directory)  # Change the current working directory to the specified 'directory'

    while check_port(bind_addr, port):
        print(f"Port {port} is already in use. Trying with port {port + 1}.")
        port += 1

    httpd = HTTPServer((bind_addr, port), CORSHTTPRequestHandler)
    print(f"Serving at http://{bind_addr}:{port}")
    inform_user(port)
    httpd.serve_forever()
