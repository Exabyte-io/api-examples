# This module defines a set of common functions which are used in other examples.
import datetime
import json
import os
import time
import urllib.request
import uuid
from typing import List, Union

from exabyte_api_client.endpoints.bank_workflows import BankWorkflowEndpoints
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.raw_properties import RawPropertiesEndpoints
from IPython.display import HTML, display
from pandas import DataFrame
from pandas.io.formats.style import Styler
from tabulate import tabulate

from . import settings

# GENERIC UTILITIES


def update_json_file_kwargs(path_to_json_file: str = "settings.json", **kwargs) -> None:
    """
    This function updates settings.json for a given kwargs if kwargs
    contains variables different from those already in json

    Args:
        path_to_json_file (str): the path to the json file to be updated
        **kwargs (dict): A dict of keyword arguments

    Returns:
        None
    """

    # 1. Assert the json file is where we think it is
    assert os.path.isfile(path_to_json_file)

    # 2. Load settings.json
    with open(path_to_json_file) as settings_json_file:
        variables = json.load(settings_json_file)

    # 3. Update json file if kwargs contains new variables
    if kwargs != variables:
        updated_variables = {**variables, **kwargs}
        with open(path_to_json_file, "w") as settings_json_file:
            json.dump(updated_variables, settings_json_file, indent=4)


def save_files(job_id: str, job_endpoint: JobEndpoints, filename_on_cloud: str, filename_on_disk: str) -> None:
    """
    Saves a file to disk, overwriting any files with the same name as filename_on_disk

    Args:
        job_id (str): ID of the job
        job_endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        filename_on_cloud (str): Name of the file on the server
        filename_on_disk (str): Name the file will be saved to

    Returns:
        None
    """
    files = job_endpoint.list_files(job_id)
    for file in files:
        if file["name"] == filename_on_cloud:
            file_metadata = file

    # Get a download URL for the CONTCAR
    signed_url = file_metadata["signedUrl"]

    # Download the contcar to memory
    server_response = urllib.request.urlopen(signed_url)

    # Write it to disk
    with open(filename_on_disk, "wb") as outp:
        outp.write(server_response.read())


# JOB UTILITIES


def get_jobs_statuses_by_ids(endpoint: JobEndpoints, job_ids: List[str]) -> List[str]:
    """
    Gets jobs statues by their IDs.

    Args:
        job_endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        job_ids (list): list of job IDs to get the status for

    Returns:
        list: list of job statuses
    """
    jobs = endpoint.list({"_id": {"$in": job_ids}}, {"fields": {"status": 1}})
    return [job["status"] for job in jobs]


def wait_for_jobs_to_finish(endpoint: JobEndpoints, job_ids: list, poll_interval: int = 10) -> None:
    """
    Waits for jobs to finish and prints their statuses.
    A job is considered finished if it is not in "pre-submission", "submitted", or, "active" status.

    Args:
        job_endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        job_ids (list): list of job IDs to wait for
        poll_interval (int): poll interval for job information in seconds. Defaults to 10.
    """
    print("Wait for jobs to finish, poll interval: {0} sec".format(poll_interval))
    while True:
        statuses = get_jobs_statuses_by_ids(endpoint, job_ids)

        errored_jobs = len([status for status in statuses if status == "error"])
        active_jobs = len([status for status in statuses if status == "active"])
        finished_jobs = len([status for status in statuses if status == "finished"])
        submitted_jobs = len([status for status in statuses if status == "submitted"])

        headers = ["TIME", "SUBMITTED-JOBS", "ACTIVE-JOBS", "FINISHED-JOBS", "ERRORED-JOBS"]
        now = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
        row = [now, submitted_jobs, active_jobs, finished_jobs, errored_jobs]
        print(tabulate([row], headers, tablefmt="grid", stralign="center"))

        if all([status not in ["pre-submission", "submitted", "active"] for status in statuses]):
            break
        time.sleep(poll_interval)


# WORKFLOW


def copy_bank_workflow_by_system_name(endpoint: BankWorkflowEndpoints, system_name: str, account_id: str) -> dict:
    """
    Copies a bank workflow with given ID into the account's workflows.

    Args:
        endpoint (endpoints.bank_workflows.BankWorkflowEndpoints): an instance of BankWorkflowEndpoints class
        system_name (str): workflow system name.
        account_id (str): ID of account to copy the bank workflow into.

    Returns:
        dict: new account's workflow
    """
    bank_workflow_id = endpoint.list({"systemName": system_name})[0]["_id"]
    return endpoint.copy(bank_workflow_id, account_id)["_id"]


# PROPERTY


def get_property_by_subworkflow_and_unit_indicies(
    endpoint: RawPropertiesEndpoints, property_name: str, job: dict, subworkflow_index: int, unit_index: int
) -> dict:
    """
    Returns the property extracted in the given unit of the job's subworkflow.

    Args:
        endpoint (endpoints.raw_properties.RawPropertiesEndpoints): an instance of RawPropertiesEndpoints class.
        property_name (str): name of property to extract.
        job (dict): job config to extract the property from.
        subworkflow_index (int): index of subworkflow to extract the property from.
        unit_index (int): index of unit to extract the property from.

    Returns:
        dict: extracted property
    """
    unit_flowchart_id = job["workflow"]["subworkflows"][subworkflow_index]["units"][unit_index]["flowchartId"]
    return endpoint.get_property(job["_id"], unit_flowchart_id, property_name)


# DISPLAY UTILITIES


def dataframe_to_html(df: DataFrame, text_align: str = "center") -> Styler:
    """
    Converts Pandas dataframe to HTML.
    See https://pandas.pydata.org/pandas-docs/stable/style.html for more information about styling.

    Args:
        df (pd.DataFrame): Pandas dataframe.
        text_align (str): text align. Defaults to center.
    """
    styles = [
        dict(selector="th", props=[("text-align", text_align)]),
        dict(selector="td", props=[("text-align", text_align)]),
    ]
    return df.style.set_table_styles(styles)


def display_JSON(
    obj: Union[dict, list], interactive_viewer: bool = settings.use_interactive_JSON_viewer, level: int = 2
) -> None:
    """
    Displays JSON, either interactively or via a text dump to Stdout.

    The interactive viewer is based on https://github.com/mljar/mercury/blob/main/mercury/widgets/json.py.

    Args:
        obj (dict): Object to display as nicely-formatted JSON
        interactive (bool): Whether to use the interactive viewer or not
    """
    if interactive_viewer:
        JSON_CSS = """.renderjson a              { text-decoration: none; }
        .renderjson .disclosure    { color: grey; font-size: 125%; }
        .renderjson .syntax        { color: grey; }
        .renderjson .string        { color: #fe46a5; }
        .renderjson .number        { color: #0f9b8e; }
        .renderjson .boolean       { color: black; }
        .renderjson .key           { color: #2684ff; }
        .renderjson .keyword       { color: gray; }
        .renderjson .object.syntax { color: gray; }
        .renderjson .array.syntax  { color: gray; }"""

        if isinstance(obj, (dict, list)):
            json_str = json.dumps(obj)
        else:
            json_str = obj
        id = str(uuid.uuid4())

        # Minimized version based on https://raw.githubusercontent.com/caldwell/renderjson/master/renderjson.js
        # Changes:
        # 1. Replace \n with \\n
        # 2. Remove module exports
        js = """var renderjson=function(){var t=function(){for(var t=[];arguments.length;)t.push(n(s(Array.prototype.shift.call(arguments)),o(Array.prototype.shift.call(arguments))));return t},n=function(){for(var t=Array.prototype.shift.call(arguments),e=0;e<arguments.length;e++)arguments[e].constructor==Array?n.apply(this,[t].concat(arguments[e])):t.appendChild(arguments[e]);return t},e=function(t,n){return t.insertBefore(n,t.firstChild),t},r=function(t,n){var e=n||Object.keys(t);for(var r in e)if(Object.hasOwnProperty.call(t,e[r]))return!1;return!0},o=function(t){return document.createTextNode(t)},s=function(t){var n=document.createElement("span");return t&&(n.className=t),n},l=function(t,n,e){var r=document.createElement("a");return n&&(r.className=n),r.appendChild(o(t)),r.href="#",r.onclick=function(t){return e(),t&&t.stopPropagation(),!1},r};function a(i,c,u,p,y){var _=u?"":c,f=function(r,a,i,c,u){var f,g=s(c),h=function(){f||n(g.parentNode,f=e(u(),l(y.hide,"disclosure",(function(){f.style.display="none",g.style.display="inline"})))),f.style.display="inline",g.style.display="none"};n(g,l(y.show,"disclosure",h),t(c+" syntax",r),l(a,null,h),t(c+" syntax",i));var d=n(s(),o(_.slice(0,-1)),g);return p>0&&"string"!=c&&h(),d};return null===i?t(null,_,"keyword","null"):void 0===i?t(null,_,"keyword","undefined"):"string"==typeof i&&i.length>y.max_string_length?f('"',i.substr(0,y.max_string_length)+" ...",'"',"string",(function(){return n(s("string"),t(null,_,"string",JSON.stringify(i)))})):"object"!=typeof i||[Number,String,Boolean,Date].indexOf(i.constructor)>=0?t(null,_,typeof i,JSON.stringify(i)):i.constructor==Array?0==i.length?t(null,_,"array syntax","[]"):f("[",y.collapse_msg(i.length),"]","array",(function(){for(var e=n(s("array"),t("array syntax","[",null,"\\n")),r=0;r<i.length;r++)n(e,a(y.replacer.call(i,r,i[r]),c+"    ",!1,p-1,y),r!=i.length-1?t("syntax",","):[],o("\\n"));return n(e,t(null,c,"array syntax","]")),e})):r(i,y.property_list)?t(null,_,"object syntax","{}"):f("{",y.collapse_msg(Object.keys(i).length),"}","object",(function(){var e=n(s("object"),t("object syntax","{",null,"\\n"));for(var r in i)var l=r;var u=y.property_list||Object.keys(i);for(var _ in y.sort_objects&&(u=u.sort()),u){(r=u[_])in i&&n(e,t(null,c+"    ","key",'"'+r+'"',"object syntax",": "),a(y.replacer.call(i,r,i[r]),c+"    ",!0,p-1,y),r!=l?t("syntax",","):[],o("\\n"))}return n(e,t(null,c,"object syntax","}")),e}))}var i=function t(e){var r=new Object(t.options);r.replacer="function"==typeof r.replacer?r.replacer:function(t,n){return n};var o=n(document.createElement("pre"),a(e,"",!1,r.show_to_level,r));return o.className="renderjson",o};return i.set_icons=function(t,n){return i.options.show=t,i.options.hide=n,i},i.set_show_to_level=function(t){return i.options.show_to_level="string"==typeof t&&"all"===t.toLowerCase()?Number.MAX_VALUE:t,i},i.set_max_string_length=function(t){return i.options.max_string_length="string"==typeof t&&"none"===t.toLowerCase()?Number.MAX_VALUE:t,i},i.set_sort_objects=function(t){return i.options.sort_objects=t,i},i.set_replacer=function(t){return i.options.replacer=t,i},i.set_collapse_msg=function(t){return i.options.collapse_msg=t,i},i.set_property_list=function(t){return i.options.property_list=t,i},i.set_show_by_default=function(t){return i.options.show_to_level=t?Number.MAX_VALUE:0,i},i.options={},i.set_icons("⊕","⊖"),i.set_show_by_default(!1),i.set_sort_objects(!1),i.set_max_string_length("none"),i.set_replacer(void 0),i.set_property_list(void 0),i.set_collapse_msg((function(t){return t+" item"+(1==t?"":"s")})),i}();"""  # noqa: E501

        display(HTML(f'<style>{JSON_CSS}</style><div id="{id}"></div>'))
        display(
            HTML(
                f'<script>{js} renderjson.set_show_to_level({str(level)}); document.getElementById("{id}").appendChild(renderjson({json_str}))</script>'  # noqa: E501
            )
        )
    else:
        print(json.dumps(obj, indent=4))
