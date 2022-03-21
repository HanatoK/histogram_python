#!/usr/bin/env python3
# from histogram import Axis
from histogram import HistogramScalar
from boltzmann_constant import boltzmann_constant_kcalmolk
import argparse


class GetTrajWeight:

    def __init__(self, position_cols, pmf_filename, kbt=300.0*boltzmann_constant_kcalmolk):
        import numpy as np
        import logging
        import copy
        self.logger = logging.getLogger(self.__class__.__name__)
        logging_handler = logging.StreamHandler()
        logging_formatter = logging.Formatter('[%(name)s %(levelname)s]: %(message)s')
        logging_handler.setFormatter(logging_formatter)
        self.logger.addHandler(logging_handler)
        self.logger.setLevel(logging.INFO)
        self.positionColumns = copy.deepcopy(position_cols)
        self.maxColumn = max(self.positionColumns)
        self.probability = None
        with open(pmf_filename, 'r') as f_pmf:
            pmf = HistogramScalar()
            pmf.read_from_stream(f_pmf)
            pmf.data = np.exp(-1.0 * pmf.data / kbt)
            self.probability = pmf
        if self.probability.get_dimension() != len(self.positionColumns):
            self.logger.warning('the number of columns selected does not match the dimension of the PMF!')

    def parse_traj(self, f_traj, f_output):
        total_lines = 0
        valid_lines = 0
        for line in f_traj:
            if line.strip().startswith('#'):
                continue
            total_lines = total_lines + 1
            tmp_fields = line.split()
            if len(tmp_fields) > self.maxColumn:
                tmp_position = [float(tmp_fields[i]) for i in self.positionColumns]
                # check if the position is in boundary
                if self.probability.is_in_grid(tmp_position):
                    valid_lines = valid_lines + 1
                    weight = self.probability[tmp_position]
                    f_output.write(' '.join(tmp_fields) + f' {weight:20.15f}\n')
                else:
                    self.logger.warning(f'position {tmp_position} is not in the boundary.')
            else:
                raise RuntimeError(f'Maximum column ({self.maxColumn}) is out of bound ({len(tmp_fields)})!')
        self.logger.info(f'Total data lines: {total_lines}')
        self.logger.info(f'Valid data lines: {valid_lines}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Print the weights of a Colvars trajectory')
    required_args = parser.add_argument_group('required named arguments')
    required_args.add_argument('--pmf', help='the PMF file', required=True)
    required_args.add_argument('--traj', nargs='+', help='the Colvars trajectory file', required=True)
    required_args.add_argument('--column', type=int, nargs='+', help='the columns in the trajectory', required=True)
    required_args.add_argument('--output', help='the output file with weights', required=True)
    parser.add_argument('--kbt', default=300.0*boltzmann_constant_kcalmolk, type=float, help='KbT')
    args = parser.parse_args()
    # all the arguments are mandatory
    get_weight_traj = GetTrajWeight(args.column, args.pmf, args.kbt)
    with open(args.output, 'w') as f_output:
        for traj_file in args.traj:
            with open(traj_file, 'r') as f_traj:
                get_weight_traj.parse_traj(f_traj, f_output)
