"""
Демонстрационная программа для FASTA Parser

Этот скрипт демонстрирует основные возможности библиотеки fasta_parser,
включая работу с классами Seq и FastaReader.
"""

import os
import sys

# Добавляем родительскую директорию в путь для импорта
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fasta_parser import Seq, FastaReader, FastaFormatError, InvalidSequenceError


def demonstrate_seq_class():
    """Демонстрация возможностей класса Seq."""
    print("=" * 60)
    print("ДЕМОНСТРАЦИЯ КЛАССА Seq")
    print("=" * 60)
    
    # Создание ДНК последовательности
    print("\n1. Работа с ДНК последовательностью:")
    dna_seq = Seq("ATGCGTACGTAGCTA", "Example DNA sequence")
    print(f"Последовательность: {dna_seq.sequence}")
    print(f"Заголовок: {dna_seq.header}")
    print(f"Длина: {len(dna_seq)}")
    print(f"Тип алфавита: {dna_seq.get_alphabet_type()}")
    print(f"GC-состав: {dna_seq.get_gc_content()}%")
    print(f"Состав: {dna_seq.get_composition()}")
    
    # Обратная комплементарная последовательность
    print(f"\nОбратная комплементарная: {dna_seq.reverse_complement().sequence}")
    
    # Трансляция
    try:
        # Создаем ДНК с длиной кратной 3 для трансляции
        coding_dna = Seq("ATGAAATTTGGATAA", "Coding sequence")
        protein = coding_dna.translate()
        print(f"Трансляция: {protein.sequence}")
    except InvalidSequenceError as e:
        print(f"Ошибка трансляции: {e}")
    
    # Работа с белковой последовательностью
    print("\n2. Работа с белковой последовательностью:")
    protein_seq = Seq("MKFGSTOP", "Example protein")
    print(f"Последовательность: {protein_seq.sequence}")
    print(f"Тип алфавита: {protein_seq.get_alphabet_type()}")
    print(f"Состав: {protein_seq.get_composition()}")
    
    # Работа с РНК последовательностью
    print("\n3. Работа с РНК последовательностью:")
    rna_seq = Seq("AUGCGUACGUAGCUA", "Example RNA")
    print(f"Последовательность: {rna_seq.sequence}")
    print(f"Тип алфавита: {rna_seq.get_alphabet_type()}")
    print(f"GC-состав: {rna_seq.get_gc_content()}%")
    
    # Строковые представления
    print("\n4. Строковые представления:")
    print("str(dna_seq):")
    print(str(dna_seq))
    print(f"repr(dna_seq): {repr(dna_seq)}")


def create_sample_fasta():
    """Создает образец FASTA файла для демонстрации."""
    sample_content = """>seq1 Homo sapiens DNA sequence
ATGCGTACGTAGCTAACGTACGTACGTACGTACGTACGTACGTACGTACG
TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG

>seq2 Mus musculus protein sequence  
MKFGSTOPQWERTYUIOPASDFGHJKLZXCVBNMQWERTYUIOPASDFGH
JKLZXCVBNMQWERTYUIOPASDFGHJKLZXCVBNM

>seq3 E. coli RNA sequence
AUGCGUACGUAGCUAACGUACGUACGUACGUACGUACGUACGUACGUACG
UACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUA

>seq4 Short DNA fragment
ATGCGTACGTAGCTA

>seq5 Another protein example
ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKL
MNPQRSTVWYACDEFGHIKLMNPQRSTVWY
"""
    
    with open("sample_sequences.fasta", "w", encoding="utf-8") as f:
        f.write(sample_content)
    
    return "sample_sequences.fasta"


def demonstrate_fasta_reader():
    """Демонстрация возможностей класса FastaReader."""
    print("\n" + "=" * 60)
    print("ДЕМОНСТРАЦИЯ КЛАССА FastaReader")
    print("=" * 60)
    
    # Создаем образец файла
    fasta_file = create_sample_fasta()
    print(f"\nСоздан образец FASTA файла: {fasta_file}")
    
    try:
        # Создание объекта для чтения
        reader = FastaReader(fasta_file)
        print(f"Успешно открыт файл: {reader.file_path}")
        
        # Проверка формата
        print(f"Файл соответствует формату FASTA: {reader.is_fasta_format()}")
        
        # Чтение всех последовательностей
        print("\n1. Все последовательности в файле:")
        for i, seq in enumerate(reader, 1):
            print(f"\nПоследовательность {i}:")
            print(f"  Заголовок: {seq.header}")
            print(f"  Тип: {seq.get_alphabet_type()}")
            print(f"  Длина: {len(seq)}")
            print(f"  Первые 50 символов: {seq.sequence[:50]}...")
        
        # Статистика файла
        print("\n2. Статистика файла:")
        stats = reader.get_file_stats()
        for key, value in stats.items():
            print(f"  {key}: {value}")
        
        # Поиск по идентификатору
        print("\n3. Поиск последовательности по ID:")
        found_seq = reader.get_sequence_by_id("seq3")
        if found_seq:
            print(f"  Найдена: {found_seq.header}")
            print(f"  Тип: {found_seq.get_alphabet_type()}")
        
        # Фильтрация по длине
        print("\n4. Фильтрация длинных последовательностей (>50 нуклеотидов):")
        def long_sequences(seq):
            return len(seq) > 50
        
        long_count = 0
        for seq in reader.filter_sequences(long_sequences):
            print(f"  {seq.header} (длина: {len(seq)})")
            long_count += 1
        print(f"  Найдено длинных последовательностей: {long_count}")
        
        # Группировка по типам
        print("\n5. Группировка по типам последовательностей:")
        for seq_type in ["DNA", "RNA", "PROTEIN"]:
            print(f"\n  {seq_type} последовательности:")
            type_count = 0
            for seq in reader.get_sequences_by_type(seq_type):
                print(f"    {seq.header}")
                type_count += 1
            print(f"    Всего: {type_count}")
        
        # Создание отфильтрованного файла
        print("\n6. Создание файла с ДНК последовательностями:")
        def dna_filter(seq):
            return seq.get_alphabet_type() == "DNA"
        
        dna_count = reader.write_filtered_fasta("dna_only.fasta", dna_filter)
        print(f"  Сохранено ДНК последовательностей: {dna_count}")
        
    except (FastaFormatError, FileNotFoundError) as e:
        print(f"Ошибка при работе с FASTA файлом: {e}")
    
    finally:
        # Очистка временных файлов
        for temp_file in ["sample_sequences.fasta", "dna_only.fasta"]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                print(f"Удален временный файл: {temp_file}")


def demonstrate_error_handling():
    """Демонстрация обработки ошибок."""
    print("\n" + "=" * 60)
    print("ДЕМОНСТРАЦИЯ ОБРАБОТКИ ОШИБОК")
    print("=" * 60)
    
    # Ошибки в классе Seq
    print("\n1. Ошибки класса Seq:")
    
    try:
        invalid_seq = Seq("")
    except InvalidSequenceError as e:
        print(f"  Пустая последовательность: {e}")
    
    try:
        invalid_chars = Seq("ATGC123XYZ")
    except InvalidSequenceError as e:
        print(f"  Недопустимые символы: {e}")
    
    try:
        protein_seq = Seq("MKFG", "protein")
        protein_seq.get_gc_content()  # Ошибка: GC только для нуклеотидов
    except InvalidSequenceError as e:
        print(f"  GC-состав для белка: {e}")
    
    # Ошибки в классе FastaReader
    print("\n2. Ошибки класса FastaReader:")
    
    try:
        reader = FastaReader("nonexistent.fasta")
    except FileNotFoundError as e:
        print(f"  Файл не найден: {e}")
    
    # Создаем некорректный FASTA файл
    with open("bad_format.fasta", "w") as f:
        f.write("This is not a FASTA file")
    
    try:
        reader = FastaReader("bad_format.fasta")
    except FastaFormatError as e:
        print(f"  Неверный формат: {e}")
    
    # Очистка
    if os.path.exists("bad_format.fasta"):
        os.remove("bad_format.fasta")


def main():
    """Главная функция демонстрации."""
    print("ДЕМОНСТРАЦИЯ БИБЛИОТЕКИ FASTA PARSER")
    print("Версия: 1.0.0")
    print("=" * 60)
    
    try:
        demonstrate_seq_class()
        demonstrate_fasta_reader()
        demonstrate_error_handling()
        
        print("\n" + "=" * 60)
        print("ДЕМОНСТРАЦИЯ ЗАВЕРШЕНА УСПЕШНО!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\nНепредвиденная ошибка: {e}")
        print("Демонстрация прервана.")


if __name__ == "__main__":
    main()