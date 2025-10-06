"""
Пользовательские исключения для FASTA Parser

Этот модуль определяет специализированные исключения для обработки ошибок
при работе с FASTA файлами и биологическими последовательностями.
"""


class FastaFormatError(Exception):
    """
    Исключение для ошибок формата FASTA файла.
    
    Вызывается при обнаружении некорректного формата FASTA файла,
    например, отсутствии заголовков или неправильной структуре записей.
    
    Args:
        message (str): Описание ошибки
        line_number (int, optional): Номер строки с ошибкой
    """
    
    def __init__(self, message, line_number=None):
        self.message = message
        self.line_number = line_number
        
        if line_number:
            super().__init__(f"Строка {line_number}: {message}")
        else:
            super().__init__(message)


class InvalidSequenceError(Exception):
    """
    Исключение для некорректных биологических последовательностей.
    
    Вызывается при обнаружении недопустимых символов в последовательности
    или других проблем с валидностью биологических данных.
    
    Args:
        message (str): Описание ошибки
        sequence (str, optional): Проблемная последовательность
        invalid_chars (set, optional): Набор некорректных символов
    """
    
    def __init__(self, message, sequence=None, invalid_chars=None):
        self.message = message
        self.sequence = sequence
        self.invalid_chars = invalid_chars
        
        full_message = message
        if invalid_chars:
            full_message += f" Некорректные символы: {sorted(invalid_chars)}"
        if sequence and len(sequence) <= 50:
            full_message += f" Последовательность: {sequence}"
            
        super().__init__(full_message)