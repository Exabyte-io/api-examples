import time
from tabulate import tabulate


def get_jobs_in_state(jobs, state):
    return [job for job in jobs if job["status"] == state]


def print_status(jobs):
    error_jobs = get_jobs_in_state(jobs, "error")
    active_jobs = get_jobs_in_state(jobs, "active")
    finished_jobs = get_jobs_in_state(jobs, "finished")
    submitted_jobs = get_jobs_in_state(jobs, "submitted")
    header = ['SUBMITTED-JOBS', 'ACTIVE-JOBS', 'FINISHED-JOBS', 'ERROR-JOBS']
    row = [[len(submitted_jobs), len(active_jobs), len(finished_jobs), len(error_jobs)]]
    print tabulate(row, header, tablefmt='grid', stralign='center')
    print "\n"


def wait_for_jobs_to_finish(job_endpoints, jobs):
    while True:
        all_jobs = job_endpoints.list(query={"_id": {"$in": [job["_id"] for job in jobs]}}, projection={"status": 1})
        print_status(all_jobs)
        if all([job["status"] not in ["pre-submission", "submitted", "active"] for job in all_jobs]): break
        time.sleep(10)
