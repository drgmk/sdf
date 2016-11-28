import os
import glob

from sdf import config as cfg

def rm_all_mn():
    """Remove all multinest files."""

    fs = glob.glob(cfg.file['sdb_root']+'masters/*/*/*mnest/*')
    removed = [os.unlink(f) for f in fs]
    
if __name__ == '__main__':

    rm_all_mn()
