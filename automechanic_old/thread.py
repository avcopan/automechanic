""" thread I/O-bound tasks through a queue with named workers
"""
from queue import Queue
from queue import Empty as EmptyQueue
from threading import Thread
from threading import Event


def tag_team_starmap(function, iterable, worker_ids):
    """ tag-team a mapping on named worker threads

    :param function: takes *args and one kwarg with key 'worker_id'
    :type function: callable
    """

    def _ifunction(idx, args, worker_id):
        ret = function(*args, worker_id=worker_id)
        return (idx, ret)

    args_queue = Queue()
    ret_queue = Queue()

    args_lst = list(iterable)
    nargs = len(args_lst)

    team = []
    for worker_id in worker_ids:
        worker = _Worker(worker_id=worker_id, function=_ifunction,
                         args_queue=args_queue, ret_queue=ret_queue)
        team.append(worker)

    for worker in team:
        worker.start()

    for idx, args in enumerate(args_lst):
        args_queue.put((idx, args))

    ret_lst = [None] * nargs
    for _ in range(nargs):
        idx, ret = ret_queue.get()
        ret_lst[idx] = ret

    return ret_lst


class _Worker(Thread):
    """ take arguments from a queue and pass them to a fucntion

    Based on Eli Bendersky's post at
    eli.thegreenplace.net/2011/12/27/python-threads-communication-and-stopping
    """

    def __init__(self, worker_id, function, args_queue, ret_queue):
        """ initialize the worker thread
        """
        self.worker_id = worker_id
        self.function = function
        self.args_queue = args_queue
        self.ret_queue = ret_queue
        self.term_signal = Event()
        super(_Worker, self).__init__(daemon=True)

    def run(self):
        """ runs when start() is called until join() sets the term_signal
        """
        while not self.term_signal.is_set():
            try:
                args = self.args_queue.get(block=True, timeout=0.05)
                ret = self.function(*args, worker_id=self.worker_id)
                self.ret_queue.put(ret)
            except EmptyQueue:
                continue

    def join(self, timeout=None):
        """ wind the worker down
        """
        self.term_signal.set()
        super(_Worker, self).join(timeout)
