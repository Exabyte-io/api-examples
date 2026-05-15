from mat3ra.notebooks_utils import display_JSON


def visualize_workflow(workflow, level: int = 2) -> None:
    """
    Visualize a workflow by displaying its JSON configuration.

    Args:
        workflow: Workflow object with a to_dict() method
        level: Expansion level for the JSON viewer (default: 2)

    Returns:
        None
    """
    workflow_config = workflow.to_dict()
    display_JSON(workflow_config, level=level)
