# pylint: disable=too-few-public-methods,no-member
"""File parser for Chemstation files

.. note:: This file parser went through a large re-write on ??? which
   changed the data structures of the resulting objects. This means
   that upon upgrading it *will* be necessary to update code. The
   re-write was done to fix some serious errors from the first
   version, like relying on the Report.TXT file for injections
   summaries. These are now fetched from the more ordered CSV files.

"""

from __future__ import print_function, unicode_literals, division

from collections import defaultdict
import codecs
import os
from itertools import islice
from io import BytesIO
import time
import struct
from struct import unpack

# The standard library csv module is no good for encoded CSV, which is
# kind of annoying
import unicodecsv as csv
import numpy

class CHFile(object):
    """Class that implementats the Agilent .ch file format version 179

    .. warning:: Not all aspects of the file header is understood, so there may and probably
       is information that is not parsed. See the method :meth:`._parse_header_status` for
       an overview of which parts of the header is understood.

    .. note:: Although the fundamental storage of the actual data has change, lots of
       inspiration for the parsing of the header has been drawn from the parser in the
       `ImportAgilent.m file <https://github.com/chemplexity/chromatography/blob/dev/
       Methods/Import/ImportAgilent.m>`_ in the `chemplexity/chromatography project
       <https://github.com/chemplexity/chromatography>`_ project. All credit for the parts
       of the header parsing that could be reused goes to the author of that project.

    Attributes:
        values (numpy.array): The internsity values (y-value) or the spectrum. The unit
            for the values is given in `metadata['units']`
        metadata (dict): The extracted metadata
        filepath (str): The filepath this object was loaded from

    """

    # Fields is a table of name, offset and type. Types 'x-time' and 'utf16' are specially
    # handled, the rest are format arguments for struct unpack
    fields = (
        ('sequence_line_or_injection', 252, UINT16),
        ('injection_or_sequence_line', 256, UINT16),
        ('start_time', 282, 'x-time'),
        ('end_time', 286, 'x-time'),
        ('version_string', 326, 'utf16'),
        ('description', 347, 'utf16'),
        ('sample', 858, 'utf16'),
        ('operator', 1880, 'utf16'),
        ('date', 2391, 'utf16'),
        ('inlet', 2492, 'utf16'),
        ('instrument', 2533, 'utf16'),
        ('method', 2574, 'utf16'),
        ('software version', 3601, 'utf16'),
        ('software name', 3089, 'utf16'),
        ('software revision', 3802, 'utf16'),
        ('units', 4172, 'utf16'),
        ('detector', 4213, 'utf16'),
        ('yscaling', 4732, ENDIAN + 'd')
    )
    # The start position of the data
    data_start = 6144
    # The versions of the file format supported by this implementation
    supported_versions = {179}

    def __init__(self, filepath):
        """Instantiate object

        Args:
            filepath (str): The path of the data file
        """
        self.filepath = filepath
        self.metadata = {}
        with open(self.filepath, 'rb') as file_:
            self._parse_header(file_)
            self.values = self._parse_data(file_)


    def _parse_header(self, file_):
        """Parse the header"""
        # Parse and check version
        length = unpack(UINT8, file_.read(1))[0]
        parsed = unpack(STRING.format(length), file_.read(length))
        version = int(parsed[0])
        if version not in self.supported_versions:
            raise ValueError('Unsupported file version {}'.format(version))
        self.metadata['magic_number_version'] = version

        # Parse all metadata fields
        for name, offset, type_ in self.fields:
            file_.seek(offset)
            if type_ == 'utf16':
                self.metadata[name] = parse_utf16_string(file_)
            elif type_ == 'x-time':
                self.metadata[name] = unpack(ENDIAN + 'f', file_.read(4))[0] / 60000
            else:
                self.metadata[name] = unpack(type_, file_.read(struct.calcsize(type_)))[0]

        # Convert date
        self.metadata['datetime'] = time.strptime(self.metadata['date'], '%d-%b-%y, %H:%M:%S')

    def _parse_header_status(self):
        """Print known and unknown parts of the header"""
        file_ = open(self.filepath, 'rb')
        # Map positions to fields for all the known fields
        knowns = {item[1]: item for item in self.fields}
        # A couple of places has a \x01 byte before a string, these we simply skip
        skips = {325, 3600}
        # Jump to after the magic number version
        file_.seek(4)

        # Initialize variables for unknown bytes
        unknown_start = None
        unknown_bytes = b''
        # While we have not yet reached the data
        while file_.tell() < self.data_start:
            current_position = file_.tell()
            # Just continue on skip bytes
            if current_position in skips:
                file_.read(1)
                continue

            # If we know about a data field that starts at this point
            if current_position in knowns:
                # If we have collected unknown bytes, print them out and reset
                if unknown_bytes != b'':
                    print('Unknown at', unknown_start, repr(unknown_bytes.rstrip(b'\x00')))
                    unknown_bytes = b''
                    unknown_start = None

                # Print out the position, type, name and value of the known value
                print('Known field at {: >4},'.format(current_position), end=' ')
                name, _, type_ = knowns[current_position]
                if type_ == 'x-time':
                    print('x-time, "{: <19}'.format(name + '"'),
                          unpack(ENDIAN + 'f', file_.read(4))[0] / 60000)
                elif type_ == 'utf16':
                    print(' utf16, "{: <19}'.format(name + '"'),
                          parse_utf16_string(file_))
                else:
                    size = struct.calcsize(type_)
                    print('{: >6}, "{: <19}'.format(type_, name + '"'),
                          unpack(type_, file_.read(size))[0])
            else:  # We do not know about a data field at this position If we have already
                # collected 4 zero bytes, assume that we are done with this unkonw field,
                # print and reset
                if unknown_bytes[-4:] == b'\x00\x00\x00\x00':
                    print('Unknown at', unknown_start, repr(unknown_bytes.rstrip(b'\x00')))
                    unknown_bytes = b''
                    unknown_start = None

                # Read one byte and save it
                one_byte = file_.read(1)
                if unknown_bytes == b'':
                    # Only start a new collection of unknown bytes, if this byte is not a
                    # zero byte
                    if one_byte != b'\x00':
                        unknown_bytes = one_byte
                        unknown_start = file_.tell() - 1
                else:
                    unknown_bytes += one_byte

        file_.close()

    def _parse_data(self, file_):
        """Parse the data"""
        # Go to the end of the file and calculate how many points 8 byte floats there are
        file_.seek(0, 2)
        n_points = (file_.tell() - self.data_start) // 8

        # Read the data into a numpy array
        file_.seek(self.data_start)
        return numpy.fromfile(file_, dtype='<d', count=n_points) * self.metadata['yscaling']

    @cached_property
    def times(self):
        """The time values (x-value) for the data set in minutes"""
        return numpy.linspace(self.metadata['start_time'], self.metadata['end_time'],
                              len(self.values))
