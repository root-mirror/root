.. _testing:


Test suite
==========

The cppyy tests live in the top-level cppyy package, can be run for
both CPython and PyPy, and exercises the full setup, including the backend.
Most tests are standalone and can be run independently, with a few exceptions
in the template tests (see file ``test_templates.py``).

To run the tests, first install cppyy by any usual means, then clone the
cppyy repo, and enter the ``test`` directory::

 $ git clone https://bitbucket.org/wlav/cppyy.git
 $ cd cppyy/test

Next, build the dictionaries, the manner of which depends on your platform.
On Linux or MacOS-X, run ``make``::

 $ make all

On Windows, run the dictionary building script::

 $ python make_dict_win32.py all

Next, make sure you have `pytest`_ installed, for example with ``pip``::

 $ python -m pip install pytest

and finally run the tests::

 $ python -m pytest -sv

On Linux and MacOS-X, all tests should succeed.
On MS Windows 32bit there are 4 failing tests, on 64bit there are 5 still
failing.


.. _`pytest`: https://docs.pytest.org/en/latest/
