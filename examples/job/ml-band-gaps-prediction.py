#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse

from tabulate import tabulate

from endpoints.jobs import JobEndpoints
from endpoints.login import LoginEndpoint
from endpoints.materials import MaterialEndpoints
from examples.utils import wait_for_jobs_to_finish
from endpoints.raw_properties import RawPropertiesEndpoints


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-H', '--host', default="platform.exabyte.io", help='RESTful API hostname')
    parser.add_argument('-P', '--port', type=int, default=443, help='RESTful API port')
    parser.add_argument('-u', '--username', required=True, help='Your Exabyte username')
    parser.add_argument('-p', '--password', required=True, help='Your Exabyte password')
    parser.add_argument('-tmi', '--train-material-id', dest="train_material_ids", action="append", required=True, help='train material ID')
    parser.add_argument('-pmi', '--predict-material-id', dest="predict_material_id", required=True, help='predict material ID')
    parser.add_argument('-oi', '--owner-id', dest="owner_id", help='owner ID')
    parser.add_argument('-wi', '--workflow-id', dest="workflow_id", help='train workflow ID')
    parser.add_argument('-pi', '--project-id', dest="project_id", help='project id')
    parser.add_argument('-j', '--job-prefix', dest="job_prefix", default="job", help='job name prefix')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    # params
    owner_id = args.owner_id
    job_prefix = args.job_prefix
    project_id = args.project_id
    workflow_id = args.workflow_id
    train_material_ids = args.train_material_ids
    predict_material_id = args.predict_material_id

    endpoint_options = LoginEndpoint.get_endpoint_options(args.host, args.port, args.username, args.password)

    # endpoints
    job_endpoints = JobEndpoints(**endpoint_options)
    material_endpoints = MaterialEndpoints(**endpoint_options)
    raw_property_endpoints = RawPropertiesEndpoints(**endpoint_options)

    # create ML Train job and submit it
    job_name = "-".join((job_prefix, "train"))
    job_config = job_endpoints.get_config(train_material_ids, workflow_id, project_id, owner_id, job_name,
                                          is_multi_material=True)
    job = job_endpoints.create(job_config)
    job_endpoints.submit(job["_id"])

    # wait for job to finish
    wait_for_jobs_to_finish(job_endpoints, [job])

    # extract predict workflow
    unit_flowchart_id = job["workflow"]["subworkflows"][0]["units"][4]["flowchartId"]
    predict_workflow = raw_property_endpoints.get_property(job["_id"], unit_flowchart_id, "workflow:ml_predict")["data"]

    # create ML Predict job and submit it
    job_name = "-".join((job_prefix, "predict"))
    job_config = job_endpoints.get_config([predict_material_id], predict_workflow["_id"], project_id, owner_id, job_name)
    job = job_endpoints.create(job_config)
    job_endpoints.submit(job["_id"])

    # wait for job to finish
    wait_for_jobs_to_finish(job_endpoints, [job])

    # extract band gaps
    predict_material = material_endpoints.get(predict_material_id)

    unit_flowchart_id = job["workflow"]["subworkflows"][0]["units"][3]["flowchartId"]
    predicted_properties = raw_property_endpoints.get_property(job["_id"], unit_flowchart_id, "predicted_properties")

    predict_material_properties = predicted_properties["data"]["values"][predict_material["exabyteId"]]
    band_gaps = next((v for v in predict_material_properties if v["name"] == "band_gaps"))

    band_gap_direct = next((v for v in band_gaps["values"] if v["type"] == "direct"), None)["value"]
    band_gap_indirect = next((v for v in band_gaps["values"] if v["type"] == "indirect"), None)["value"]

    # print results
    header = ["DIRECT-GAP", "INDIRECT-GAP"]
    rows = [[band_gap_direct, band_gap_indirect]]
    print tabulate(rows, header, tablefmt='grid', stralign='center')

