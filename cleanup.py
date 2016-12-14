import os
import glob
import argparse

from sdf import db
from sdf import config as cfg

def rm_all_mn(sample='zpo_cal_'):
    """Remove all multinest files."""

    # get the targets from a sample
    ids = db.sample_targets(sample)
    fs = []
    for id in ids:
        fs += glob.glob( cfg.file['sdb_root']+'masters/'\
                        +id+'/*/*mnest/*' )

    # to remove everything
#    fs = glob.glob(cfg.file['sdb_root']+'masters/*/*/*mnest/*')

    removed = [os.unlink(f) for f in fs]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Remove multinest files')
    parser.add_argument('--sample','-s',nargs='+',
                        help='Sample(s) to clean up')
    args = parser.parse_args()

    if args.sample is not None:
        for s in args.sample:
            rm_all_mn(s)
