# Script to collect MAD-X beam and MAD-X lattice input.

def parse(**kwargs):
    import argparse
    parser = argparse.ArgumentParser(**kwargs)
    parser.add_argument('latticefile', type=str, help="MAD-X lattice filename")
    parser.add_argument('--beamfile', dest='beamfile', default='', type=str, help="MAD-X beam input filename")
    parser_namespace = parser.parse_args()