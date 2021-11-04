"""
Argparse class used by run to parse arguments
"""

__author__ = "S.J.Bouwman"
__version__ = "1.0"

import argparse
import os


def _filecheck(path) -> None:
    """
    Checks if file exists, if not raises Error
    :param path: path to file
    :return: None
    """
    exists = os.path.exists(path)
    if not exists:
        raise FileExistsError(f"File: {path} does not exist")


class ArgParser:
    """
    ArgParse class that parses arguments and saves them as variables so they can easily be accessed
    later in the process.
    """
    ref_file: str
    target_file: str
    output: str

    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-r', '--ref', type=str, required=True,
                            help="Filepath to file containg "
                                 "reference proteins in fasta format")

        parser.add_argument('-t', '--target', type=str, required=True,
                            help="Filepath to file containing target "
                                 "proteins in fasta format")

        parser.add_argument('-o', '--output', type=str, required=False,
                            help="If provided output will be written to file instead")

        self.parsed = parser.parse_args()

        _filecheck(self.parsed.ref)
        self.ref_file = self.parsed.ref

        _filecheck(self.parsed.target)
        self.target_file = self.parsed.target

        self.output = self.parsed.output
