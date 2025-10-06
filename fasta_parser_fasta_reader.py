"""
Простая реализация класса FastaReader
"""

from fasta_parser.seq import Seq

class FastaReader:
    """Читает FASTA-файл и возвращает Seq через итератор."""

    def __init__(self, filepath):
        self.filepath = filepath

    def __iter__(self):
        header = None
        seq_lines = []
        with open(self.filepath, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        yield Seq("".join(seq_lines), header)
                    header = line[1:]
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if header:
                yield Seq("".join(seq_lines), header)
