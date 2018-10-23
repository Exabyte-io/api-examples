import json
import argparse

from endpoints.login import LoginEndpoint
from endpoints.materials import MaterialEndpoints

HOST = 'platform.exabyte.io'
PORT = 443
SECURE = True

MATERIAL = {
    "name": "TEST MATERIAL",
    "basis": {
        "elements": [
            {
                "id": 1,
                "value": "Si"
            },
            {
                "id": 2,
                "value": "Si"
            }
        ],
        "coordinates": [
            {
                "id": 1,
                "value": [
                    0,
                    0,
                    0
                ]
            },
            {
                "id": 2,
                "value": [
                    0.25,
                    0.25,
                    0.25
                ]
            }
        ],
        "units": "crystal",
        "name": "basis"
    },
    "lattice": {
        "type": "FCC",
        "a": 3.867,
        "b": 3.867,
        "c": 3.867,
        "alpha": 60,
        "beta": 60,
        "gamma": 60,
        "units": {
            "length": "angstrom",
            "angle": "degree"
        },
        "vectors": {
            "a": [
                3.867,
                0,
                0
            ],
            "b": [
                1.9335000000000004,
                3.348920236434424,
                0
            ],
            "c": [
                1.9335000000000004,
                1.1163067454781415,
                3.1573922784475164
            ],
            "name": "lattice vectors",
            "alat": 1,
            "units": "angstrom"
        }
    },
    "tags": [
        "REST API"
    ]
}


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
    material = materials_endpoint.create(MATERIAL)
    print json.dumps(material, indent=4)
