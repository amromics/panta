import multiprocessing
from threading import Thread
from six.moves.queue import Queue, Empty
import math
import platform
import sys

# alternative: https://github.com/jreese/aiomultiprocess
def run_parallel(target, tasks, n_workers = None, use_processes = False):
    """
    From James: https://github.com/nesp-tsr3-1/tsx/blob/master/tsx/util.py
    Runs tasks in parallel
    `target` is a function
    `tasks` is a list of argument tuples passed to the function. If `target` only takes one argument, then it doesn't need to
    be wrapped in a tuple.
    A generator yielding the result of each task is returned, in the form of a (result, error) tuple which allows errors to be
    handled. The generator must be consumed in order to ensure all tasks are processed.
    Results are yielded in the order they are completed, which is generally not the same as the order in which they are supplied.
    Example::
        def do_hard_work(a, b):
            ...
        tasks = [(1,2), (5,2), (3,4) ....]
        for result, error in run_parallel(do_hard_work, tasks):
            print result
    A pool of worker threads (or processes if `use_processes = True`) is used to process the tasks.
    Threads may not always be able to achieve parallelism due to Python GIL.
    If using processes, be careful not to use shared global resources such as database connection pools in the target function.
    The number of workers defaults to the number of cpu cores as reported by `multiprocessing.cpu_count`, but can be set
    using the `n_workers` parameter.
    """
    if n_workers is None:
        n_workers = multiprocessing.cpu_count()

    # Multiprocessing has issues on Windows
    if platform.system() == 'Windows':
        use_processes = False

    Q = multiprocessing.Queue if use_processes else Queue

    # Setup queues

    work_q = Q()
    result_q = Q()
    # Helper to get next item from queue without constantly blocking
    def next(q):
        while True:
            try:
                return q.get(True, 1) # Get with timeout so thread isn't constantly blocked
            except Empty:
                pass
            except:
                # logger.exception("Exception getting item from queue")
                raise

    # Setup worker threads
    def worker(work_q, result_q):
        # faulthandler.enable()
        while True:
            task = next(work_q)
            if task is None:
                break
            try:
                result_q.put((target(*task), None))
            except:
                e = sys.exc_info()[0]
                # logger.exception("Exception in worker")
                result_q.put((None, e))

    for i in range(0, n_workers):
            if use_processes:
                p = multiprocessing.Process(target = worker, args = (work_q, result_q))
                p.daemon = True # Kill process if parent terminates early
                p.start()
            else:
                t = Thread(target = worker, args = (work_q, result_q))
                t.daemon = True
                t.start()

    # Feed in tasks and yield results
    i = 0
    for task in tasks:
        # Tasks must always be tuples
        if type(task) != tuple:
            task = (task,)
        work_q.put(task)
        i += 1
        # Start getting results once all threads have something to do
        if i > n_workers:
            yield next(result_q)
            i -= 1

    # Signal threads to stop
    for j in range(0, n_workers):
        work_q.put(None)

    # Finish collecting results
    while i > 0:
        yield next(result_q)
        i -= 1