#!/usr/bin/env python3

class BuildHistogramFromTraj:

    def __init__(self, json_file, position_cols=None):
        import copy
        import logging
        from histogram import HistogramScalar
        self.logger = logging.getLogger(self.__class__.__name__)
        logging_handler = logging.StreamHandler()
        logging_formatter = logging.Formatter('[%(name)s %(levelname)s]: %(message)s')
        logging_handler.setFormatter(logging_formatter)
        self.logger.addHandler(logging_handler)
        self.logger.setLevel(logging.INFO)
        self.histogram = HistogramScalar.from_json_file(json_file)
        if position_cols is None:
            position_cols = list(range(0, self.histogram.get_dimension()))
        self.positionColumns = copy.deepcopy(position_cols)
        self.maxColumn = max(self.positionColumns)
        self.logger.info(f'Using columns: {self.positionColumns}')

    def read_traj(self, f_traj):
        total_lines = 0
        valid_lines = 0
        for line in f_traj:
            if line.startswith('#'):
                continue
            tmp_fields = line.split()
            if len(tmp_fields) > self.maxColumn:
                total_lines += 1
                tmp_position = [float(tmp_fields[i]) for i in self.positionColumns]
                if self.histogram.is_in_grid(tmp_position):
                    self.histogram[tmp_position] += 1.0
                    valid_lines += 1
                else:
                    self.logger.warning(f'Position {tmp_position} is not in the boundary!')
        self.logger.info(f'Total data lines: {total_lines}')
        self.logger.info(f'Valid data lines: {valid_lines}')

    def get_histogram(self):
        return self.histogram


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Build histogram from colvars trajectories')
    required_args = parser.add_argument_group('required named arguments')
    required_args.add_argument('--axis', help='json file to setup axes')
    required_args.add_argument('--traj', nargs='+', help='the Colvars trajectory file', required=True)
    required_args.add_argument('--output', help='the output file with weights', required=True)
    parser.add_argument('--columns', type=int, nargs='+', help='the columns in the trajectory')
    args = parser.parse_args()
    build_histogram = BuildHistogramFromTraj(json_file=args.axis, position_cols=args.columns)
    for traj_file in args.traj:
        with open(traj_file, 'r') as f_traj:
            build_histogram.read_traj(f_traj=f_traj)
    with open(args.output, 'w') as f_output:
        build_histogram.get_histogram().write_to_stream(stream=f_output)
