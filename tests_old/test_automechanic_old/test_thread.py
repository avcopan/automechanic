""" test the automechanic_old.thread module
"""
import time
from automechanic_old import thread


def test__tag_team_starmap():
    """ test thread.tag_team_starmap
    """

    def func(wait_time, ret_val, worker_id):
        print("Waiting {:d} seconds to return {:s} from worker {:s}"
              .format(wait_time, ret_val, worker_id))
        time.sleep(wait_time)
        print("Returning {:s} from worker {:s} after waiting {:d} seconds"
              .format(ret_val, worker_id, wait_time))
        return ret_val

    wait_times = (15, 4, 16, 1, 9, 10, 8, 7)
    ret_vals = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')

    iterable = tuple(zip(wait_times, ret_vals))
    wait_times, ret_vals = tuple(zip(*iterable))

    start_time = time.time()
    worker_ids = ('b444', 'b445', 'b446', 'b447', 'b448')
    thread.tag_team_starmap(func, iterable, worker_ids)
    run_time = time.time() - start_time
    print(run_time)
    assert run_time < 17


if __name__ == '__main__':
    test__tag_team_starmap()
