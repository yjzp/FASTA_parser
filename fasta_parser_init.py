"""
FASTA Parser - Python библиотека для работы с биологическими последовательностями

Модули:
    seq: Класс Seq для работы с биологическими последовательностями
    fasta_reader: Класс FastaReader для чтения FASTA файлов
    exceptions: Пользовательские исключения

Автор: Биологический факультет
Версия: 1.0.0
"""

from .seq import Seq
from .fasta_reader import FastaReader
from .exceptions import FastaFormatError, InvalidSequenceError

__version__ = "1.0.0"
__author__ = "Biology Department"
__all__ = ["Seq", "FastaReader", "FastaFormatError", "InvalidSequenceError"]