import os
abspath = os.path.abspath('.')

# read pyproject.toml as the single reference for the version number:
with open(f'{abspath}/pyproject.toml', 'r') as f:
    for line in f.readlines():
        if 'version' in line:
            __version__ = line.split('=')[1][:-1]
            break