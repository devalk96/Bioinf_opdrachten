"""
Main computing class, more info can be found in README.md
"""

__author__ = "S.J.Bouwman"
__version__ = "1.0"

import sys
import os
from collections import Counter

from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline


def _perform_msa(filename: str, outfile: str) -> None:
    """
    Uses ClustalO commandline to perform MSA
    :param filename: filepath to file containing reference proteins (in fasta format)
    :param outfile: filepath where created fasta will be exported to
    :return: None
    """
    clustalomega_cline = ClustalOmegaCommandline("clustalo", infile=filename,
                                                 outfile=outfile,
                                                 force=True, outfmt="fasta")
    clustalomega_cline()


def _count_reference(reference) -> list:
    """
    Counts the unique amino acids in each column
    :param reference: list containing tuples of amino acids [('M', '-'), ('M', 'K'), ....]
    :return: list[[letter, count],....]
    """
    counted_reference = []
    for column in reference:
        letters = list(Counter(column).keys())
        occurrence = list(Counter(column).values())

        # Remove gaps
        if "-" in letters:
            ind = letters.index("-")
            letters.pop(ind)
            occurrence.pop(ind)

        counted_reference.append(list(zip(letters, occurrence)))
    return counted_reference


class MSAParser:
    """MSAParser parses protein files and creates MSA from them, wherafter it creates a score for
    SNPS"""
    _tmp_folder_path: str = "data/tmp/"
    _protein_count: int = 0
    _protein_headers: list = []
    _output_provided = None

    def generate_msa_standalone(self, refprotein_path: str,
                                target_gene_protein: str,
                                output=None) -> None:
        """
        Creates a MSA from ref proteins and target proteins
        and returns score of these target proteins. More info in README.md
        :param refprotein_path: filepath to file containing reference proteins (in fasta format)
        :param target_gene_protein: filepath to file containing target proteins (in fasta format)
        :return: None
        """
        # Setup for run
        tmp_fasta_path: str = self._tmp_folder_path + "msa_multi.fasta"
        tmp_msa_path: str = self._tmp_folder_path + "mas_complete.fasta"
        self._output_provided = output
        if os.path.exists(tmp_fasta_path):
            os.remove(tmp_fasta_path)

        # Create combined file so msa might be performed
        with open(tmp_fasta_path, "w") as combined_f:
            for file in [target_gene_protein, refprotein_path]:
                with open(file, "r") as read_f:
                    for line in read_f:
                        combined_f.write(line)

        _perform_msa(tmp_fasta_path, tmp_msa_path)

        # Get amount of proteins and their names
        with open(target_gene_protein, "r") as stream:
            for line in stream:
                if line.startswith(">"):
                    self._protein_count += 1
                    self._protein_headers.append(line.replace("\n", ""))

        self._calculate_conservation(tmp_msa_path)

    def _calculate_conservation(self, file) -> None:
        """
        Calculates the conservation of individual amino acids in the target sequences.
        More info about score calculation can be found in the README.md
        :param file: MSA file containing both target proteins and reference proteins
        :return: None
        """
        try:
            pars = AlignIO.read(file, format="fasta")

        except ValueError as exception:
            print("Error: File doesn't seem to contain alligned data", file=sys.stderr)
            raise exception

        # Get number of unique letters and their count
        counted_reference = _count_reference(list(zip(*pars[self._protein_count:])))

        # Construct full protein sequences from zip data
        sequences = self._create_seq(list(zip(*pars[:self._protein_count])))

        results = {}
        for seq_nr, test_seq in enumerate(sequences):
            summary = {}
            deleterious = []
            possible_deleterious = []
            scores = {}

            for column, _ in enumerate(test_seq):
                current_protein = test_seq[column]

                # Gaps should be skipped
                if current_protein == "-":
                    continue

                counted = counted_reference[column]

                letters = [x[0] for x in counted]
                letter_count = [x[1] for x in counted]

                # Set default score to maximum
                score = 9

                if current_protein in letters:
                    count = letter_count[letters.index(current_protein)]
                    score = 9 - (count / sum(letter_count) * 9)

                scores[column] = score

                if score >= 8.1:
                    deleterious.append(column)
                    continue

                if score >= 6.75:
                    possible_deleterious.append(column)
                    continue

            summary["deleterious"] = deleterious
            summary["possible_deleterious"] = possible_deleterious
            summary["scores"] = scores
            results[seq_nr] = summary
        self._clear_tmp()
        self.report(results)

    def _create_seq(self, target_protein) -> list:
        """
        Transforms column data to full protein sequences
        :param target_protein: list of sets containing amino acids
        :return: list of protein sequences
        """
        sequences = []
        for protein_nr in range(self._protein_count):
            protein = ""
            for col in target_protein:
                protein += col[protein_nr]
            sequences.append(protein)
        return sequences

    def _clear_tmp(self) -> None:
        """Empties temporary folder"""
        for file in os.listdir(self._tmp_folder_path):
            os.remove(self._tmp_folder_path + file)

    def report(self, results: dict[int: dict]) -> None:
        """
        Provides user feedback
        :param results: conservation calculation results
        :return: None
        """

        if self._output_provided:
            outfile = open(self._output_provided, "w")
        else:
            outfile = sys.stdout

        for seq_id in results:
            deleterious_c = results[seq_id]["deleterious"]
            pos_deleterious_c = results[seq_id]["possible_deleterious"]

            print(f"Protein sequence ({seq_id}):\n{self._protein_headers[seq_id]}: ", file=outfile)
            if deleterious_c:
                print(f"Contains: {len(deleterious_c)} deleterious (score > 8.1) "
                      f"snps at index: {' '.join([str(x) for x in deleterious_c])}", file=outfile)
            else:
                print("Contains no deleterious snps", file=outfile)

            if pos_deleterious_c:
                print(f"Contains: {len(pos_deleterious_c)} possible deleterious (score > 6.75) "
                      f"snps at index: {' '.join([str(x) for x in pos_deleterious_c])}",
                      file=outfile)

            tot_amino = len(results[seq_id]['scores'])
            tot_dele = len(deleterious_c) + len(pos_deleterious_c)
            percentage = (tot_dele / tot_amino) * 100

            print(f"Percentage of amino acids marked as "
                  f"(possible) deleterious: "
                  f"{len(deleterious_c) + len(pos_deleterious_c)} "
                  f"of {len(results[seq_id]['scores'])} amino acids ({round(percentage, 3)}%)\n",
                  file=outfile)

        if self._output_provided:
            outfile.close()
            path = self._output_provided
            self._output_provided = None
            self.report(results)
            print(f"Output written to {path}")
