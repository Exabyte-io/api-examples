import json
import argparse

from endpoints.login import LoginEndpoint

HOST = 'platform.exabyte.io'
PORT = 443
SECURE = True


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--username', required=True, help='username')
    parser.add_argument('-p', '--password', required=True, help='password')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    login_endpoint = LoginEndpoint(HOST, PORT, args.username, args.password, secure=SECURE)
    auth_params = login_endpoint.login()
    print json.dumps(auth_params, indent=4)
