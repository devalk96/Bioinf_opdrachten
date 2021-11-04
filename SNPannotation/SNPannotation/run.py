#!/usr/bin/python3

"""
Main class to run the SNPAnnotation, for more info visit the README.md
"""

__author__ = "S.J.Bouwman"
__version__ = "1.0"

from arg_parser import ArgParser
from msa_parser import MSAParser


class Run:
    """SNPAnnotation run class"""
    parsed_args: dict
    parser: ArgParser = ArgParser()
    msa = MSAParser()

    def __init__(self):
        self.parsed_args = self.parser.parsed
        self.msa.generate_msa_standalone(self.parser.ref_file,
                                         self.parser.target_file,
                                         self.parser.output)


if __name__ == '__main__':
    start = Run()
