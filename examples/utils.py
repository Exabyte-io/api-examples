#!/usr/bin/env python
# coding: utf-8

# This module defines a set of functions which are used in other examples. 

# In[ ]:


import time
from tabulate import tabulate


def get_jobs_statuses_by_ids(endpoint, job_ids):
    """
    Get jobs statues by their IDs.
    
    Args:
        endpoint (endpoints.jobs.JobEndpoints): an instance of job endpoints class
        job_ids (list): list of job IDs to get the status for
    
    Returns:
        list: list of job statuses
    """
    jobs = endpoint.list({"_id": {"$in": job_ids}}, {"fields": {"status": 1}})
    return [job["status"] for job in jobs]


def wait_for_jobs_to_finish(endpoint, job_ids, pulling_interval=60):
    """
    Wait for jobs to finish and print their statuses.
    A job is considered finished if it is not in "pre-submission", "submitted", or, "active" status.
    
    Args:
        endpoint (endpoints.jobs.JobEndpoints): an instance of job endpoints class
        job_ids (list): list of job IDs to wait for
        pulling_interval (int): job pulling interval in seconds. Defaults to 60.
    """
    while True:
        statuses = get_jobs_statuses_by_ids(endpoint, job_ids)

        errored_jobs = len([status for status in statuses if status == "error"])
        active_jobs = len([status for status in statuses if status == "active"])
        finished_jobs = len([status for status in statuses if status == "finished"])
        submitted_jobs = len([status for status in statuses if status == "submitted"])

        headers = ["SUBMITTED-JOBS", "ACTIVE-JOBS", "FINISHED-JOBS", "ERRORED-JOBS"]
        row = [submitted_jobs, active_jobs, finished_jobs, errored_jobs]
        print tabulate([row], headers, tablefmt='grid', stralign='center')

        if all([status not in ["pre-submission", "submitted", "active"] for status in statuses]): break
        time.sleep(pulling_interval)

