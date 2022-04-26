#!/usr/bin/env python3
# from histogram import Axis
from histogram import HistogramScalar
from boltzmann_constant import boltzmann_constant_kcalmolk
import argparse


class GetTrajWeightEGABF:

    def __init__(self, position_cols, pmf_filenames, kbt=300.0*boltzmann_constant_kcalmolk):
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
        self.pmfs = list()
        self.kbt = kbt
        self.logger.warning('The weight will not be normalized!')
        # for egABF simulations, all PMFs must be 1D, and the number of PMFs should match the number of position columns
        if len(position_cols) != len(pmf_filenames):
            raise RuntimeError('The number of PMFs mismatches the number of columns')
        for pmf_filename in pmf_filenames:
            with open(pmf_filename, 'r') as f_pmf:
                self.logger.info(f'Reading {pmf_filename}')
                pmf = HistogramScalar()
                pmf.read_from_stream(f_pmf)
                if pmf.get_dimension() != 1:
                    raise RuntimeError(f'PMF should be 1D in egABF reweighting (file: {pmf_filename})')
                self.pmfs.append(copy.deepcopy(pmf))

    def parse_traj(self, f_traj, f_output):
        import numpy as np
        total_lines = 0
        valid_lines = 0
        for line in f_traj:
            if line.strip().startswith('#'):
                continue
            total_lines = total_lines + 1
            tmp_fields = line.split()
            if len(tmp_fields) > self.maxColumn:
                sum_delta_G = 0
                position_in_grid = True
                tmp_positions = list()
                for pmf, col_index in zip(self.pmfs, self.positionColumns):
                    pos = float(tmp_fields[col_index])
                    tmp_positions.append(pos)
                    tmp_position = [pos]
                    if pmf.is_in_grid(tmp_position):
                        valid_lines = valid_lines + 1
                        sum_delta_G += pmf[tmp_position]
                    else:
                        sum_delta_G = 0
                        position_in_grid = False
                        continue
                weight = np.exp(-1.0 * sum_delta_G / self.kbt)
                if position_in_grid:
                    f_output.write(' '.join(tmp_fields) + f' {weight:20.15f}\n')
                else:
                    self.logger.warning(f'position {tmp_positions} is not in the boundary.')
            else:
                raise RuntimeError(f'Maximum column ({self.maxColumn}) is out of bound ({len(tmp_fields)})!')
        self.logger.info(f'Total data lines: {total_lines}')
        self.logger.info(f'Valid data lines: {valid_lines}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Print the weights of a Colvars trajectory from an egABF simulation')
    required_args = parser.add_argument_group('required named arguments')
    required_args.add_argument('--pmfs', nargs='+', help='egABF 1D PMF file(s)', required=True)
    required_args.add_argument('--traj', nargs='+', help='Colvars trajectory files', required=True)
    required_args.add_argument('--columns', type=int, nargs='+', help='the columns in the trajectory matching the CVs '
                                                                      'of the PMF', required=True)
    required_args.add_argument('--output', help='the output file with weights', required=True)
    parser.add_argument('--kbt', default=300.0*boltzmann_constant_kcalmolk, type=float, help='KbT')
    args = parser.parse_args()
    get_weight_traj = GetTrajWeightEGABF(args.columns, args.pmfs, args.kbt)
    with open(args.output, 'w') as f_output:
        for traj_file in args.traj:
            with open(traj_file, 'r') as f_traj:
                get_weight_traj.parse_traj(f_traj, f_output)
