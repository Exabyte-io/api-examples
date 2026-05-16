from pandas import DataFrame
from pandas.io.formats.style import Styler


def dataframe_to_html(df: DataFrame, text_align: str = "center") -> Styler:
    """
    Converts Pandas dataframe to HTML.

    Args:
        df (pd.DataFrame): Pandas dataframe.
        text_align (str): text align. Defaults to center.
    """
    styles = [
        dict(selector="th", props=[("text-align", text_align)]),
        dict(selector="td", props=[("text-align", text_align)]),
    ]
    return df.style.set_table_styles(styles)
