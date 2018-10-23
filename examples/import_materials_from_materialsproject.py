import json
import argparse

from endpoints.login import LoginEndpoint
from endpoints.materials import MaterialEndpoints


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-H', '--host', default="platform.exabyte.io", help='RESTful API hostname')
    parser.add_argument('-P', '--port', type=int, default=443, help='RESTful API port')
    parser.add_argument('-S', '--insecure', action="store_true", default=False, help='Whether to use SSL')
    parser.add_argument('-u', '--username', required=True, help='Your Exabyte username')
    parser.add_argument('-p', '--password', required=True, help='Your Exabyte password')
    parser.add_argument('-k', '--key', required=True, help='materialsproject key')
    parser.add_argument('-m', '--material', dest="material_ids", action="append", required=True, help='material ID')
    parser.add_argument('-t', '--tag', dest="tags", action="append", help='material tag')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    login_endpoint = LoginEndpoint(args.host, args.port, args.username, args.password, secure=not args.insecure)
    response = login_endpoint.login()
    account_id = response["X-Account-Id"]
    auth_token = response["X-Auth-Token"]

    materials_endpoint = MaterialEndpoints(args.host, args.port, account_id, auth_token, secure=not args.insecure)
    materials = materials_endpoint.import_from_materialsproject(args.key, args.material_ids, args.tags or [])
    print json.dumps(materials, indent=4)
