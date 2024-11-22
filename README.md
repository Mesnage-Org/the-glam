# glam
GLycan Analysis Module

# Documentation

Docstrings are used to automatically generate documentation for `glam`.
This project will follow the docstring standard embraced and defined by [SciPy/NumPy](https://numpydoc.readthedocs.io/en/latest/format.html).
These docstrings look something like:

```
"""Gets and prints the spreadsheet's header columns

Parameters
----------
file_loc : str
    The file location of the spreadsheet
print_cols : bool, optional
    A flag used to print the columns to the console (default is False)

Returns
-------
list
    a list of strings representing the header columns
"""
```

# Related Work
https://github.com/SugarPy/SugarPy
https://github.com/ursgal/ursgal
https://github.com/pyQms/pyqms

# Findings
The MeHex mass was wrong, so structures with more MeHex were more incorrect in their masses.
Also, something horrific happened with the core structure — it's off by 12+ Da
