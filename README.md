# FASTA Parser

Python библиотека для работы с биологическими последовательностями в формате FASTA.

## Основные возможности

- **Класс Seq**: Работа с биологическими последовательностями (ДНК, РНК, белки)
- **Класс FastaReader**: Эффективное чтение и анализ FASTA файлов
- **Валидация данных**: Проверка формата и корректности последовательностей

## Установка

### Из репозитория

```bash
git clone https://github.com/yjzp/FASTA_parser.git
cd fasta-parser
```

### Для разработки

```bash
git clone https://github.com/username/FASTA_parser.git
cd fasta-parser
```

## Быстрый старт

```python
from fasta_parser import Seq, FastaReader

# Работа с последовательностями
seq = Seq("ATGCGTAG", "Example DNA")
print(seq.get_alphabet_type())  # "DNA"
print(seq.get_gc_content())     # 50.0
print(len(seq))                 # 8

# Чтение FASTA файла
reader = FastaReader("sequences.fasta")
for sequence in reader:
    print(f">{sequence.header}")
    print(f"Type: {sequence.get_alphabet_type()}")
    print(f"Length: {len(sequence)}")
```

## Документация

### Сборка HTML документации

```bash
cd docs
sphinx-build -b html . _build/html
```

Документация будет доступна в `docs/_build/html/index.html`

### Быстрая сборка

```bash
cd docs
make html  # На Unix/Linux/macOS
```

## Тестирование

### Запуск всех тестов

```bash
pytest
```

### Запуск с покрытием кода

```bash
pytest --cov=fasta_parser --cov-report=html
```

### Запуск конкретных тестов

```bash
pytest tests/test_seq.py
pytest tests/test_fasta_reader.py
```

## Демонстрация

Запустите демонстрационную программу:

```bash
python examples/demo.py
```

Или через установленную команду:

```bash
fasta-demo
```

## Структура проекта

```
fasta-parser/
├── fasta_parser/           # Основной пакет
│   ├── __init__.py        # Инициализация пакета
│   ├── seq.py             # Класс Seq
│   ├── fasta_reader.py    # Класс FastaReader
│   └── exceptions.py      # Пользовательские исключения
├── tests/                 # Модульные тесты
│   ├── test_seq.py        # Тесты для класса Seq
│   └── test_fasta_reader.py # Тесты для FastaReader
├── examples/              # Примеры использования
│   └── demo.py           # Демонстрационная программа
├── docs/                  # HTML документация
│   ├── conf.py           # Конфигурация Sphinx
│   └── index.rst         # Главная страница документации
├── README.md             # Этот файл
├── requirements.txt      # Зависимости
├── setup.py             # Скрипт установки
├── pytest.ini          # Конфигурация pytest
└── LICENSE              # Лицензия
```

## API Reference

### Класс Seq

Основной класс для работы с биологическими последовательностями.

**Основные методы:**
- `get_alphabet_type()` - определение типа последовательности (DNA/RNA/PROTEIN)
- `get_composition()` - анализ состава последовательности
- `get_gc_content()` - вычисление GC-состава для нуклеотидов
- `reverse_complement()` - обратная комплементарность для ДНК
- `translate()` - трансляция ДНК в белок
- `find_orfs()` - поиск открытых рамок считывания

**Пример:**
```python
# Создание ДНК последовательности
dna = Seq("ATGAAATTTGGATAA", "Coding sequence")

# Анализ последовательности
print(f"Тип: {dna.get_alphabet_type()}")
print(f"GC-состав: {dna.get_gc_content()}%")
print(f"Состав: {dna.get_composition()}")

# Обратная комплементарная последовательность
rev_comp = dna.reverse_complement()
print(f"Обратная комплементарная: {rev_comp.sequence}")

# Трансляция в белок
protein = dna.translate()
print(f"Белковая последовательность: {protein.sequence}")
```

### Класс FastaReader

Класс для эффективного чтения FASTA файлов.

**Основные методы:**
- `read_sequences()` - генератор чтения последовательностей
- `get_file_stats()` - статистика файла
- `get_sequence_by_id()` - поиск по идентификатору
- `filter_sequences()` - фильтрация последовательностей
- `write_filtered_fasta()` - запись отфильтрованных данных

**Пример:**
```python
# Чтение файла
reader = FastaReader("sequences.fasta")

# Итерация по всем последовательностям
for seq in reader:
    print(f"Заголовок: {seq.header}")
    print(f"Длина: {len(seq)}")
    print(f"Тип: {seq.get_alphabet_type()}")

# Получение статистики
stats = reader.get_file_stats()
print(f"Количество последовательностей: {stats['sequence_count']}")
print(f"Общая длина: {stats['total_length']}")

# Фильтрация длинных последовательностей
def long_sequences(seq, min_length=100):
    return len(seq) >= min_length

count = reader.write_filtered_fasta(
    "long_sequences.fasta", 
    long_sequences,
    min_length=500
)
print(f"Сохранено {count} длинных последовательностей")
```

## Примеры использования

### Анализ GC-состава

```python
from fasta_parser import FastaReader

reader = FastaReader("sequences.fasta")
gc_contents = []

for seq in reader.get_sequences_by_type("DNA"):
    try:
        gc_content = seq.get_gc_content()
        gc_contents.append(gc_content)
        print(f"{seq.header}: {gc_content}%")
    except Exception as e:
        print(f"Ошибка для {seq.header}: {e}")

if gc_contents:
    avg_gc = sum(gc_contents) / len(gc_contents)
    print(f"Средний GC-состав: {avg_gc:.2f}%")
```

### Поиск и трансляция ORF

```python
from fasta_parser import Seq

# ДНК последовательность с открытой рамкой считывания
dna_seq = Seq("ATGAAATTTGGATAAGGCATGCGTACG", "Example with ORFs")

# Поиск ORF
orfs = dna_seq.find_orfs(min_length=15)

for i, (start, end, frame, orf_sequence) in enumerate(orfs, 1):
    print(f"ORF {i}:")
    print(f"  Позиция: {start}-{end}")
    print(f"  Рамка: {frame}")
    print(f"  Последовательность: {orf_sequence}")
    
    # Трансляция ORF
    orf = Seq(orf_sequence, f"ORF{i}")
    protein = orf.translate()
    print(f"  Белок: {protein.sequence}")
```

## Обработка ошибок

Библиотека определяет специальные исключения:

- `FastaFormatError` - ошибки формата FASTA файла
- `InvalidSequenceError` - ошибки в последовательностях

```python
from fasta_parser import FastaReader, Seq
from fasta_parser import FastaFormatError, InvalidSequenceError

try:
    # Попытка прочитать некорректный файл
    reader = FastaReader("bad_file.txt")
except FastaFormatError as e:
    print(f"Ошибка формата: {e}")
    
try:
    # Создание некорректной последовательности
    seq = Seq("ATGC123XYZ")
except InvalidSequenceError as e:
    print(f"Ошибка последовательности: {e}")
```

## Производительность

Библиотека оптимизирована для работы с большими файлами:

- **Генераторы**: Чтение по одной записи для экономии памяти
- **Валидация**: Проверка формата перед обработкой
- **Эффективность**: Минимальные зависимости и быстрая обработка

## Разработка

### Запуск тестов

```bash
# Все тесты
pytest

# С подробным выводом
pytest -v

# Конкретный тест
pytest tests/test_seq.py::TestSeq::test_gc_content

# С покрытием кода
pytest --cov=fasta_parser
```

### Форматирование кода

```bash
# Форматирование с black
black fasta_parser/ tests/ examples/

# Проверка стиля с flake8
flake8 fasta_parser/ tests/ examples/
```

### Сборка документации

```bash
cd docs
sphinx-build -b html . _build/html
```

## Лицензия

MIT License - см. файл LICENSE для подробностей.
