import sys
import argparse
from collections import defaultdict
import traceback
def fx_hamming_error(args):

    try:

        counts = defaultdict(lambda: [0,0])

        for i, line in enumerate(open(args.hapmer_count)):
            if i == 0: continue

            its = line.split()
            ss = (int(its[2]), int(its[3]))
            counts[its[0]][0] += max(ss)
            counts[its[0]][1] += min(ss)

        for k,v in counts.items():
            print(k, "(%):", 100*v[1]/(v[0]+v[1]))

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()
        exit(-1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--hapmer_count", help="xxx.hapmer.count in merqury result",  dest="hapmer_count", type=str)
    args = parser.parse_args()

    fx_hamming_error(args)