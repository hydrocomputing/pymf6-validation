
from pathlib import Path


def rmtree(root):
    """Remove directory tree"""
    for p in root.iterdir():
        if p.is_dir():
            rmtree(p)
        else:
            p.unlink()
    root.rmdir()


def get_full_model_path(path_name):
    return Path(__file__).parent / 'models' / path_name
