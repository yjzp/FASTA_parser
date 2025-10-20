# FASTA Parser

Python библиотека для работы с биологическими последовательностями в формате FASTA.

## Основные возможности

- **Класс Seq**: Работа с биологическими последовательностями
- **Класс FastaReader**: Эффективное чтение и анализ FASTA файлов
- **Валидация данных**: Проверка формата и корректности последовательностей

## Установка

```bash
git clone https://github.com/yjzp/FASTA_parser.git
cd fasta_parser
```

## Быстрый старт

```python
from fasta_parser import Seq, FastaReader

# Работа с последовательностями
seq = Seq("ATGCGTAG", "Example DNA")
print(seq.get_alphabet_type())
print(seq.get_gc_content())
print(len(seq))

# Чтение FASTA файла
reader = FastaReader("sequences.fasta")
for sequence in reader:
    print(f">{sequence.header}")
    print(f"Type: {sequence.get_alphabet_type()}")
    print(f"Length: {len(sequence)}")
```

### Класс Seq

Основной класс для работы с биологическими последовательностями.

**Основные методы:**
- `get_alphabet_type()` - определение типа последовательности (DNA/RNA/PROTEIN)
- `get_composition()` - анализ состава последовательности
- `get_gc_content()` - вычисление GC-состава для нуклеотидов
- `reverse_complement()` - обратная комплементарность для ДНК
- `translate()` - трансляция ДНК в белок
- `find_orfs()` - поиск открытых рамок считывания


### Класс FastaReader

Класс для эффективного чтения FASTA файлов.

**Основные методы:**
- `read_sequences()` - генератор чтения последовательностей
- `get_file_stats()` - статистика файла
- `get_sequence_by_id()` - поиск по идентификатору
- `filter_sequences()` - фильтрация последовательностей
- `write_filtered_fasta()` - запись отфильтрованных данных


## Лицензия

MIT License - см. файл LICENSE для подробностей.
