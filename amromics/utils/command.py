import os
import logging
logger = logging.getLogger(__name__)
def run_command(cmd, timing_log=None):
    """
    Run a command line, return the returning code of the command
    :param cmd:
    :param timing_log:
    :return:
    """
    if timing_log is not None:
        cmd = '/usr/bin/time --append -v -o {} bash -c "{}"'.format(timing_log, cmd)
    logger.info('Running "{}'.format(cmd))
    ret = os.system(cmd)
    return ret
