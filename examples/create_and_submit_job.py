import json
import argparse

from endpoints.jobs import JobEndpoints
from endpoints.login import LoginEndpoint

HOST = 'platform.exabyte.io'
PORT = 443
SECURE = True

config = {
    "_material": {
        "_id": ""
    },
    "workflow": {
        "_id": ""
    },
    "name": "TEST JOB"
}


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--username', required=True, help='username')
    parser.add_argument('-p', '--password', required=True, help='password')
    parser.add_argument('-m', '--material', required=True, help='material ID')
    parser.add_argument('-w', '--workflow', required=True, help='workflow ID')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    login_endpoint = LoginEndpoint(HOST, PORT, args.username, args.password, secure=SECURE)
    response = login_endpoint.login()
    account_id = response["X-Account-Id"]
    auth_token = response["X-Auth-Token"]

    config['_material']['_id'] = args.material
    config['workflow']['_id'] = args.workflow

    job_endpoints = JobEndpoints(HOST, PORT, account_id=account_id, auth_token=auth_token, secure=SECURE)
    job = job_endpoints.create(config)
    job_endpoints.submit(job['_id'])
    print json.dumps(job_endpoints.get(job['_id']), indent=4)
