# Module to launch python tutorials.
# It serves the purpose of enabling batch mode,
# since PyROOT does not parse anymore by default
# the command line arguments

import sys, os

if len(sys.argv) < 2:
    raise RuntimeError('Please specify tutorial file path as only argument of this script')

# Get path to the tutorial file
file_path = os.path.expanduser(sys.argv[1])
sys.argv.remove(sys.argv[1])

if os.path.exists(file_path):
    # Module needs to run as main.
    # Some tutorials have "if __name__ == '__main__'"
    module_name = "__main__"

    # Ensure batch mode (some tutorials use graphics)
    import ROOT
    ROOT.gROOT.SetBatch(True)

    # Prevent import from generating .pyc files in source directory
    sys.dont_write_bytecode = True

    # Execute test
    if sys.version_info >= (3,5):
        import importlib.util
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
    else:
        import imp
        imp.load_module(module_name, open(file_path, 'r'), file_path, ('.py','r',1))
else:
    raise RuntimeError('Cannot execute test, {} is not a valid tutorial file path'.format(file_path))
