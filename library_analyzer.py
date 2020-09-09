import time
import pandas as pd
import sys


def timer(func):

    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        sys.stderr.write(f'Execution time of "{func.__name__}" function: '
                         f'{round(time.time() - start, 2)} sec\n')
        return result

    return wrapper


@timer
def main():
    junctions = pd.read_table(sys.argv[1], header=None)
    junctions.columns = ['junction_id', 'offset', 'F1', 'R1', 'F2', 'R2']
    junctions = junctions.set_index(['junction_id', 'offset'])
    # column sums
    col_sums = junctions.sum(axis=0)
    print(f'Sum of counts for each read type:\n{col_sums.to_string()}')
    paired = all(col_sums > 0)
    print(f'Guessing the data is {("single", "pair")[paired]}-end')
    # computing sums of pairs
    junctions['F1+R2'] = junctions['F1'] + junctions['R2']
    junctions['F2+R1'] = junctions['F2'] + junctions['R1']
    # correlation
    print(f'\n{"-" * 60}\n')
    correlation_matrix = junctions.corr()
    print(f'Correlation matrix:\n{correlation_matrix.to_string()}')
    if paired:
        stranded = bool(correlation_matrix.loc['F1+R2', 'F2+R1'] <= 0.5)
    else:
        if all(col_sums[:2]):
            stranded = bool(correlation_matrix.loc['F1', 'R1'] <= 0.5)
        else:
            stranded = bool(correlation_matrix.loc['F2', 'R2'] <= 0.5)
    print(f'Guessing the data is {("un", "")[stranded]}stranded')
    print(f'\n{"-" * 60}\n')
    # offset distribution
    print('Offset distribution:')
    print(junctions.groupby(by='offset').sum().iloc[:, :4].to_string())


if __name__ == '__main__':
    main()
