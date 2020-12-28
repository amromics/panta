

import re


def valid_id(sid):
    """
    Check if string can be a valid id, that is contains only alphanumerical and underscore characters
    Parameters
    ----------
    sid

    Returns
    -------

    """
    return re.match(r'^[.\w]+$', sid)