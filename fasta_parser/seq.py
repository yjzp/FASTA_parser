"""
Реализация класса Seq
"""

class Seq:
    """Класс для работы с биологических последовательностями."""

    def __init__(self, seq, header=""):
        if not seq.strip():
            raise ValueError("Последовательность не может быть пустой")
        # Убираем пробелы и переносы строк
        self.sequence = seq.replace("\n", "").replace(" ", "").upper()
        self.header = header

    def __str__(self):
        return f">{self.header}\n{self.sequence}"

    def __len__(self):
        return len(self.sequence)

    def alphabet_type(self):
        """Определяет тип последовательности: DNA, RNA или PROTEIN."""
        s = set(self.sequence)
        if s.issubset(set("ATCG")):
            return "DNA"
        if s.issubset(set("AUCG")):
            return "RNA"
        if s.issubset(set("ACDEFGHIKLMNPQRSTVWY")):
            return "PROTEIN"
        return "UNKNOWN"

    def gc_content(self):
        """Вычисляет процент GC для DNA/RNA."""
        atype = self.alphabet_type()
        if atype not in ("DNA", "RNA"):
            raise ValueError("GC-состав доступен только для DNA или RNA")
        g = self.sequence.count("G")
        c = self.sequence.count("C")
        return round(100 * (g + c) / len(self.sequence), 2)
