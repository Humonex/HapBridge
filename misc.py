import sys
import gzip
import logging

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt = '%Y-%m-%d %H:%M:%S')

def create_logger(name):
    return logging.getLogger(name)

def open_file(fname):
    if fname == "-":
        return sys.stdin
    elif fname.endswith(".gz"):
        return gzip.open(fname, "rt") 
    else:
        return open(fname, "r")



