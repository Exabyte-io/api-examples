#!/usr/bin/env python
# coding: utf-8

# This module defines a set of functions which are used in other examples. 

# In[ ]:


import time
from tabulate import tabulate


def get_jobs_statuses_by_ids(endpoint, job_ids):
    """
    Get jobs statues by their IDs.
    """
    jobs = endpoint.list({"_id": {"$in": job_ids}}, {"fields": {"status": 1}})
    return [job["status"] for job in jobs]


def wait_for_jobs_to_finish(endpoint, job_ids, job_pull_interval=60):
    """
    Wait for jobs to finish. 
    A job is considered finished if it is not in "pre-submission", "submitted", or, "active" status.
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
        time.sleep(job_pull_interval)

