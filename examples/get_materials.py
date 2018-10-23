import json
import argparse

from endpoints.login import LoginEndpoint
from endpoints.materials import MaterialEndpoints

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
    response = login_endpoint.login()
    account_id = response["X-Account-Id"]
    auth_token = response["X-Auth-Token"]

    materials_endpoint = MaterialEndpoints(HOST, PORT, account_id=account_id, auth_token=auth_token, secure=SECURE)
    materials = materials_endpoint.list()
    print json.dumps(materials, indent=4)
