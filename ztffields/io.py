import sys
if sys.version_info[:2] >= (3, 10):
    from importlib.resources import files  # Python 3.10+
else:
    from importlib_resources import files  # External backport

