import time
import sys


def timer(func):
    """
    Execution time profiling decorator
    :param func:
    :return:
    """
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        sys.stderr.write(f'Execution time of "{func.__name__}" function: '
                         f'{round(time.time() - start, 2)} sec\n')
        return result

    return wrapper
