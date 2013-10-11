# PyMOL-RR

PyMOL-RR is a PyMOL plugin that can be used for visualizing residue-residue contact predictions in the format required by the CASP experiment in PyMOL. You can filter predictions using a GUI and map them onto a provided PDB structure. 


## Usage
The plugin takes a contact prediction in CASP RR format and maps the highest-scoring predictions onto a protein structure.

## Requirements
For easier installation, I recommend PyMOL >= 1.5.0.6 . Older versions will not have the plugin manager described in the installation procedure.

## Installation
Checkout the project somewhere on your hard disk. In PyMOL, go to Plugin -> Plugin Manager -> Settings and ensure that the root of the PyMOL-RR git repository (NOT THE `casprr` SUBDIRECTORY!) is in the plugin search path.

I recommend also installing the [ZeroResidues](http://www.pymolwiki.org/index.php/Zero_residues) script from the PyMOLWiki to renumber residues so they are starting from index 0.

## License
PyMOL-RR is licensed under the MIT License.

Copyright (c) 2012 Stefan Seemayer

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
