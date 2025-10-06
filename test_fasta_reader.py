"""
Тесты для класса FastaReader

Модульные тесты для проверки функциональности класса FastaReader
и его методов чтения FASTA файлов.
"""

import unittest
import os
import tempfile
import sys

# Добавляем родительскую директорию в путь
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fasta_parser.fasta_reader import FastaReader
from fasta_parser.seq import Seq
from fasta_parser.exceptions import FastaFormatError


class TestFastaReader(unittest.TestCase):
    """Тесты для класса FastaReader."""
    
    def setUp(self):
        """Подготовка тестовых данных."""
        # Создаем временные FASTA файлы для тестирования
        self.temp_dir = tempfile.mkdtemp()
        
        # Корректный FASTA файл
        self.valid_fasta_content = """>seq1 First DNA sequence
ATGCGTACGTAGCTA
ACGTACGTACGTACG
>seq2 Second protein sequence
MKFGSTOP
>seq3 Third RNA sequence
AUGCGUACGUAGCUA
"""
        
        self.valid_fasta_file = os.path.join(self.temp_dir, "valid.fasta")
        with open(self.valid_fasta_file, 'w') as f:
            f.write(self.valid_fasta_content)
            
        # Некорректный файл (не FASTA)
        self.invalid_file = os.path.join(self.temp_dir, "invalid.txt")
        with open(self.invalid_file, 'w') as f:
            f.write("This is not a FASTA file")
            
        # Пустой FASTA файл
        self.empty_fasta_file = os.path.join(self.temp_dir, "empty.fasta")
        with open(self.empty_fasta_file, 'w') as f:
            f.write(">seq1 Empty sequence\n")
            
    def tearDown(self):
        """Очистка временных файлов."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_file_not_found(self):
        """Тест обработки несуществующего файла."""
        with self.assertRaises(FileNotFoundError):
            FastaReader("nonexistent.fasta")
    
    def test_invalid_format(self):
        """Тест обработки файла неверного формата."""
        with self.assertRaises(FastaFormatError):
            FastaReader(self.invalid_file)
    
    def test_valid_format_detection(self):
        """Тест определения корректного FASTA формата."""
        reader = FastaReader(self.valid_fasta_file)
        self.assertTrue(reader.is_fasta_format())
        
    def test_read_sequences_generator(self):
        """Тест чтения последовательностей через генератор."""
        reader = FastaReader(self.valid_fasta_file)
        sequences = list(reader.read_sequences())
        
        self.assertEqual(len(sequences), 3)
        
        # Проверяем первую последовательность
        seq1 = sequences[0]
        self.assertEqual(seq1.header, "seq1 First DNA sequence")
        self.assertEqual(seq1.sequence, "ATGCGTACGTAGCTAACGTACGTACGTACG")
        self.assertEqual(seq1.get_alphabet_type(), "DNA")
        
        # Проверяем вторую последовательность
        seq2 = sequences[1]
        self.assertEqual(seq2.header, "seq2 Second protein sequence")
        self.assertEqual(seq2.sequence, "MKFGSTOP")
        
        # Проверяем третью последовательность
        seq3 = sequences[2]
        self.assertEqual(seq3.header, "seq3 Third RNA sequence")
        self.assertEqual(seq3.get_alphabet_type(), "RNA")
    
    def test_iterator_interface(self):
        """Тест интерфейса итератора."""
        reader = FastaReader(self.valid_fasta_file)
        sequences = []
        
        for seq in reader:
            sequences.append(seq)
            
        self.assertEqual(len(sequences), 3)
        self.assertIsInstance(sequences[0], Seq)
    
    def test_get_sequence_count(self):
        """Тест подсчета количества последовательностей."""
        reader = FastaReader(self.valid_fasta_file)
        count = reader.get_sequence_count()
        self.assertEqual(count, 3)
    
    def test_get_file_stats(self):
        """Тест сбора статистики файла."""
        reader = FastaReader(self.valid_fasta_file)
        stats = reader.get_file_stats()
        
        # Проверяем основные поля статистики
        self.assertEqual(stats['sequence_count'], 3)
        self.assertGreater(stats['total_length'], 0)
        self.assertGreater(stats['max_length'], stats['min_length'])
        self.assertIn('DNA', stats['alphabet_types'])
        self.assertIn('RNA', stats['alphabet_types'])
        self.assertIn('PROTEIN', stats['alphabet_types'])
        
    def test_get_sequence_by_id(self):
        """Тест поиска последовательности по идентификатору."""
        reader = FastaReader(self.valid_fasta_file)
        
        # Поиск существующей последовательности
        seq = reader.get_sequence_by_id("seq2")
        self.assertIsNotNone(seq)
        self.assertEqual(seq.header, "seq2 Second protein sequence")
        
        # Поиск по части заголовка
        seq_by_desc = reader.get_sequence_by_id("protein")
        self.assertIsNotNone(seq_by_desc)
        
        # Поиск несуществующей последовательности
        not_found = reader.get_sequence_by_id("nonexistent")
        self.assertIsNone(not_found)
    
    def test_filter_sequences(self):
        """Тест фильтрации последовательностей."""
        reader = FastaReader(self.valid_fasta_file)
        
        # Фильтр по длине
        def long_sequences(seq):
            return len(seq) > 10
            
        long_seqs = list(reader.filter_sequences(long_sequences))
        self.assertGreater(len(long_seqs), 0)
        
        # Все отфильтрованные последовательности должны быть длинными
        for seq in long_seqs:
            self.assertGreater(len(seq), 10)
    
    def test_get_sequences_by_type(self):
        """Тест группировки по типу последовательностей."""
        reader = FastaReader(self.valid_fasta_file)
        
        # Получаем ДНК последовательности
        dna_sequences = list(reader.get_sequences_by_type("DNA"))
        self.assertGreater(len(dna_sequences), 0)
        
        for seq in dna_sequences:
            self.assertEqual(seq.get_alphabet_type(), "DNA")
        
        # Получаем белковые последовательности
        protein_sequences = list(reader.get_sequences_by_type("PROTEIN"))
        self.assertGreater(len(protein_sequences), 0)
        
        for seq in protein_sequences:
            self.assertEqual(seq.get_alphabet_type(), "PROTEIN")
    
    def test_write_filtered_fasta(self):
        """Тест записи отфильтрованных последовательностей."""
        reader = FastaReader(self.valid_fasta_file)
        output_file = os.path.join(self.temp_dir, "filtered.fasta")
        
        # Фильтр для ДНК последовательностей
        def dna_filter(seq):
            return seq.get_alphabet_type() == "DNA"
        
        count = reader.write_filtered_fasta(output_file, dna_filter)
        
        # Проверяем что файл создан и содержит правильное количество записей
        self.assertTrue(os.path.exists(output_file))
        self.assertGreater(count, 0)
        
        # Читаем созданный файл
        filtered_reader = FastaReader(output_file)
        filtered_sequences = list(filtered_reader.read_sequences())
        
        self.assertEqual(len(filtered_sequences), count)
        for seq in filtered_sequences:
            self.assertEqual(seq.get_alphabet_type(), "DNA")
    
    def test_extract_subsequences(self):
        """Тест извлечения подпоследовательностей."""
        reader = FastaReader(self.valid_fasta_file)
        
        # Извлекаем первые 10 символов
        subseqs = list(reader.extract_subsequences(0, 10))
        
        for subseq in subseqs:
            self.assertEqual(len(subseq), 10)
            self.assertIn("[0:10]", subseq.header)


class TestFastaReaderEdgeCases(unittest.TestCase):
    """Тесты граничных случаев для FastaReader."""
    
    def setUp(self):
        """Подготовка тестовых данных."""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Очистка временных файлов."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_empty_lines_handling(self):
        """Тест обработки пустых строк."""
        fasta_with_empty_lines = """>seq1 Test sequence

ATGCGTAG

>seq2 Another sequence


MKFG

"""
        test_file = os.path.join(self.temp_dir, "empty_lines.fasta")
        with open(test_file, 'w') as f:
            f.write(fasta_with_empty_lines)
        
        reader = FastaReader(test_file)
        sequences = list(reader.read_sequences())
        
        self.assertEqual(len(sequences), 2)
        self.assertEqual(sequences[0].sequence, "ATGCGTAG")
        self.assertEqual(sequences[1].sequence, "MKFG")
    
    def test_single_sequence(self):
        """Тест файла с одной последовательностью."""
        single_seq_fasta = """>single_seq Test
ATGCGTAG"""
        
        test_file = os.path.join(self.temp_dir, "single.fasta")
        with open(test_file, 'w') as f:
            f.write(single_seq_fasta)
        
        reader = FastaReader(test_file)
        sequences = list(reader.read_sequences())
        
        self.assertEqual(len(sequences), 1)
        self.assertEqual(sequences[0].sequence, "ATGCGTAG")
    
    def test_very_long_sequence(self):
        """Тест с очень длинной последовательностью."""
        long_sequence = "A" * 10000
        long_seq_fasta = f">long_seq Very long sequence\n{long_sequence}"
        
        test_file = os.path.join(self.temp_dir, "long.fasta")
        with open(test_file, 'w') as f:
            f.write(long_seq_fasta)
        
        reader = FastaReader(test_file)
        sequences = list(reader.read_sequences())
        
        self.assertEqual(len(sequences), 1)
        self.assertEqual(len(sequences[0]), 10000)
    
    def test_multiline_sequence(self):
        """Тест последовательности на нескольких строках."""
        multiline_fasta = """>multiline Test sequence
ATGCGTAG
CGTACGTA
CGTAGCTA"""
        
        test_file = os.path.join(self.temp_dir, "multiline.fasta")
        with open(test_file, 'w') as f:
            f.write(multiline_fasta)
        
        reader = FastaReader(test_file)
        sequences = list(reader.read_sequences())
        
        self.assertEqual(len(sequences), 1)
        expected_sequence = "ATGCGTAGCGTACGTACGTAGCTA"
        self.assertEqual(sequences[0].sequence, expected_sequence)


if __name__ == '__main__':
    unittest.main()