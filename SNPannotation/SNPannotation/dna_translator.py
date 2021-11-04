from Bio.Seq import Seq
import os
import sys


class DNATranslator:
    _dna_seq: str = None
    _reading_frames: list[str] = []

    def __init__(self):
        raise NotImplemented("Not implemented yet")

    def parse_from_string(self, dna_string) -> str:
        self._set_dna_seq(dna_string)
        return self._dna_seq

    def parse_from_file(self, path):
        if not os.path.exists(path):
            raise FileExistsError(f"File to DNA sequence does not exists!\nPath: {path}")

        with open(path, "r") as stream:
            tempseq = ""
            for line in stream:
                if not line.startswith(">"):
                    tempseq += line

        tempseq = tempseq.replace("\n", "")
        return self.parse_from_string(tempseq)

    def get_reading_frame_protein(self) -> list[Seq]:
        self._empty_dnaseq_check()
        return [Seq(dna).translate() for dna in self._reading_frames]

    def get_reading_frame_dna(self) -> list[Seq]:
        self._empty_dnaseq_check()
        return [Seq(dna) for dna in self._reading_frames]

    def _empty_dnaseq_check(self):
        if not self._dna_seq:
            print(f"No DNA string provided, will return empty Seq() object",
                  file=sys.stderr)

    def _create_reading_frame(self):
        for frame in range(0, 3):
            length = 3 * ((len(self._dna_seq) - frame) // 3)
            self._reading_frames.append(self._dna_seq[frame:frame + length])

        for x in self._reading_frames:
            print(f"Length: {len(x)} --> {len(x) % 3}\t{x}")

    def get_dna_seq(self):
        return self._dna_seq

    def _set_dna_seq(self, seq: str):
        self._dna_seq = seq
        self._create_reading_frame()


if __name__ == '__main__':
    dn = DNATranslator()
    dn.parse_from_file("data/examplegene.txt")
    proteins = dn.get_reading_frame_protein()

