"""
Тесты для класса Seq

Модульные тесты для проверки функциональности класса Seq
и его методов работы с биологическими последовательностями.
"""

import unittest
import sys
import os

# Добавляем родительскую директорию в путь
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fasta_parser.seq import Seq
from fasta_parser.exceptions import InvalidSequenceError


class TestSeq(unittest.TestCase):
    """Тесты для класса Seq."""
    
    def setUp(self):
        """Подготовка тестовых данных."""
        self.dna_seq = Seq("ATGCGTAG", "Test DNA")
        self.rna_seq = Seq("AUGCGUAG", "Test RNA")
        self.protein_seq = Seq("MKFG", "Test protein")
        
    def test_initialization(self):
        """Тест инициализации объекта Seq."""
        seq = Seq("ATGC", "Test sequence")
        self.assertEqual(seq.sequence, "ATGC")
        self.assertEqual(seq.header, "Test sequence")
        
        # Тест с пустым заголовком
        seq_no_header = Seq("ATGC")
        self.assertEqual(seq_no_header.sequence, "ATGC")
        self.assertEqual(seq_no_header.header, "")
        
    def test_invalid_sequence(self):
        """Тест обработки некорректных последовательностей."""
        with self.assertRaises(InvalidSequenceError):
            Seq("")  # Пустая последовательность
            
        with self.assertRaises(InvalidSequenceError):
            Seq("   ")  # Только пробелы
            
        with self.assertRaises(InvalidSequenceError):
            Seq("ATGC123")  # Недопустимые символы
    
    def test_sequence_cleaning(self):
        """Тест очистки последовательности от пробелов."""
        seq = Seq("  ATG CGT AG  \n", "Test")
        self.assertEqual(seq.sequence, "ATGCGTAG")
        
    def test_len(self):
        """Тест определения длины последовательности."""
        self.assertEqual(len(self.dna_seq), 8)
        self.assertEqual(len(Seq("A")), 1)
        
    def test_str_repr(self):
        """Тест строковых представлений."""
        str_repr = str(self.dna_seq)
        self.assertIn(">Test DNA", str_repr)
        self.assertIn("ATGCGTAG", str_repr)
        
        repr_str = repr(self.dna_seq)
        self.assertIn("Seq(", repr_str)
        self.assertIn("ATGCGTAG", repr_str)
        
    def test_equality(self):
        """Тест сравнения последовательностей."""
        seq1 = Seq("ATGC", "Header1")
        seq2 = Seq("ATGC", "Header2")
        seq3 = Seq("CGTA", "Header1")
        
        self.assertEqual(seq1, seq2)  # Одинаковые последовательности
        self.assertNotEqual(seq1, seq3)  # Разные последовательности
        self.assertNotEqual(seq1, "ATGC")  # Сравнение с строкой
        
    def test_get_alphabet_type(self):
        """Тест определения типа алфавита."""
        self.assertEqual(self.dna_seq.get_alphabet_type(), "DNA")
        self.assertEqual(self.rna_seq.get_alphabet_type(), "RNA")
        self.assertEqual(self.protein_seq.get_alphabet_type(), "PROTEIN")
        
        # Тест с расширенными алфавитами
        extended_dna = Seq("ATGCRYSWN")
        self.assertEqual(extended_dna.get_alphabet_type(), "DNA")
        
        # Неизвестный тип
        unknown_seq = Seq("XYZ", validate=False)
        # Это создаст ошибку, так как валидация включена по умолчанию
        
    def test_get_composition(self):
        """Тест анализа состава последовательности."""
        composition = self.dna_seq.get_composition()
        expected = {'A': 2, 'T': 2, 'G': 3, 'C': 1}
        self.assertEqual(composition, expected)
        
    def test_gc_content(self):
        """Тест вычисления GC-состава."""
        # ДНК: ATGCGTAG - 4 GC из 8 = 50%
        self.assertEqual(self.dna_seq.get_gc_content(), 50.0)
        
        # РНК
        self.assertEqual(self.rna_seq.get_gc_content(), 50.0)
        
        # Ошибка для белковой последовательности
        with self.assertRaises(InvalidSequenceError):
            self.protein_seq.get_gc_content()
            
        # Тест с последовательностью без GC
        at_seq = Seq("ATATAT")
        self.assertEqual(at_seq.get_gc_content(), 0.0)
        
        # Тест с только GC
        gc_seq = Seq("GCGCGC")
        self.assertEqual(gc_seq.get_gc_content(), 100.0)
        
    def test_reverse_complement(self):
        """Тест создания обратной комплементарной последовательности."""
        # ATGCGTAG -> CTACGCAT
        rev_comp = self.dna_seq.reverse_complement()
        self.assertEqual(rev_comp.sequence, "CTACGCAT")
        self.assertIn("reverse complement", rev_comp.header)
        
        # Ошибка для не-ДНК последовательности
        with self.assertRaises(InvalidSequenceError):
            self.rna_seq.reverse_complement()
            
        with self.assertRaises(InvalidSequenceError):
            self.protein_seq.reverse_complement()
            
    def test_translate(self):
        """Тест трансляции ДНК в белок."""
        # ATG AAA TTT GGA TAA (15 нуклеотидов, 5 кодонов)
        coding_seq = Seq("ATGAAATTTGGATAA", "Coding DNA")
        protein = coding_seq.translate()
        
        # ATG=M, AAA=K, TTT=F, GGA=G, TAA=*
        expected_protein = "MKFG*"
        self.assertEqual(protein.sequence, expected_protein)
        self.assertIn("translated", protein.header)
        
        # Ошибка для последовательности неправильной длины
        with self.assertRaises(InvalidSequenceError):
            Seq("ATGAA").translate()  # 5 нуклеотидов, не кратно 3
            
        # Ошибка для не-ДНК последовательности
        with self.assertRaises(InvalidSequenceError):
            self.protein_seq.translate()
            
    def test_find_orfs(self):
        """Тест поиска открытых рамок считывания."""
        # Последовательность с ORF: ATGAAATTTGGATAA
        # ATG-AAA-TTT-GGA-TAA (старт ATG, стоп TAA)
        test_seq = Seq("ATGAAATTTGGATAA", "Test ORF")
        orfs = test_seq.find_orfs(min_length=12)  # Минимум 12 нуклеотидов
        
        self.assertEqual(len(orfs), 1)
        start, end, frame, orf_seq = orfs[0]
        self.assertEqual(start, 0)
        self.assertEqual(end, 15)
        self.assertEqual(frame, 1)
        self.assertEqual(orf_seq, "ATGAAATTTGGATAA")
        
        # Тест без ORF
        no_orf_seq = Seq("AAATTTCCCGGG")
        orfs_empty = no_orf_seq.find_orfs()
        self.assertEqual(len(orfs_empty), 0)
        
        # Ошибка для не-ДНК последовательности
        with self.assertRaises(InvalidSequenceError):
            self.protein_seq.find_orfs()


class TestSeqEdgeCases(unittest.TestCase):
    """Тесты граничных случаев для класса Seq."""
    
    def test_case_insensitive(self):
        """Тест нечувствительности к регистру."""
        seq_upper = Seq("ATGC")
        seq_lower = Seq("atgc")
        seq_mixed = Seq("AtGc")
        
        self.assertEqual(seq_upper.sequence, "ATGC")
        self.assertEqual(seq_lower.sequence, "ATGC")
        self.assertEqual(seq_mixed.sequence, "ATGC")
        
    def test_whitespace_handling(self):
        """Тест обработки пробелов."""
        seq_with_spaces = Seq("ATG CGT AG")
        self.assertEqual(seq_with_spaces.sequence, "ATGCGTAG")
        
        seq_with_newlines = Seq("ATG\nCGT\nAG")
        self.assertEqual(seq_with_newlines.sequence, "ATGCGTAG")
        
    def test_single_character(self):
        """Тест с однобуквенными последовательностями."""
        single_a = Seq("A")
        self.assertEqual(len(single_a), 1)
        self.assertEqual(single_a.get_alphabet_type(), "DNA")
        self.assertEqual(single_a.get_gc_content(), 0.0)
        
    def test_ambiguous_nucleotides(self):
        """Тест с неоднозначными нуклеотидами."""
        amb_seq = Seq("ATGCRYWSMK")
        self.assertEqual(amb_seq.get_alphabet_type(), "DNA")
        
        # GC content должен учитывать только G и C
        self.assertEqual(amb_seq.get_gc_content(), 20.0)  # 2 из 10


if __name__ == '__main__':
    unittest.main()