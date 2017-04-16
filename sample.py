import sys
if sys.version_info[0] < 3:
    raise Exception("Please use python3, not {}!".format(sys.version))

import argparse

from sdf import plotting
from sdf import tables
from sdf import templates
from sdf import www
from sdf import utils
from sdf import config as cfg


# run from the command line
if __name__ == "__main__":
    """Generate HTML sample pages to browse database."""
    
    # inputs
    parser = argparse.ArgumentParser(description='Update web pages')
    parser.add_argument('--tables','-t',action='store_true',
                        help='Update sample tables')
    parser.add_argument('--plots','-p',action='store_true',
                        help='Update sample plots')
    parser.add_argument('--calibration','-l',action='store_true',
                        help='Update calibration plots')
    parser.add_argument('--cleanup','-c',action='store_true',
                        help='Remove unneccessary sample dirs')
    args = parser.parse_args()

    if args.tables:
        print("Updating sample tables")
        tables.sample_tables()

    if args.plots:
        print("Updating sample plots")
        plotting.flux_size_plots()
        plotting.sample_plots()

    if args.calibration:
        print("Updating calibration plots")
        for sample in cfg.www['cal_samples']:
            plotting.calibration(sample=sample)

    if args.cleanup:
        print("Cleaning up")
        www.cleanup_sample_dirs()
        www.cleanup_calibration_dirs()
