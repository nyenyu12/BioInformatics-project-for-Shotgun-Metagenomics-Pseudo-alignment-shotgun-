"""
@file datafile.py
@brief Provides file I/O utilities for loading and dumping record data.
@details Defines custom exceptions and classes for parsing FASTA/FASTQ files into record containers.
"""

import gzip
import pickle
from typing import Set, Optional
from records import (
    FASTARecordContainer,
    FASTAQRecordContainer,
    NoRecordsInData,
    RecordContainer,
)

## ===========================================================
## Custom Exceptions
## ===========================================================

class InvalidExtensionError(Exception):
    """
    @brief Exception raised when the file extension is not valid.
    """
    def __init__(self, message: str = "") -> None:
        super().__init__(message)

class NoRecordsInDataFile(Exception):
    """
    @brief Exception raised when a file does not contain any valid records.
    """
    def __init__(self, message: str = "") -> None:
        super().__init__(message)

## ===========================================================
## Class: DataFile
## ===========================================================

class DataFile:
    """
    @brief Abstract base class for data files.
    """
    EXTENSIONS: Optional[Set[str]] = None

    def __init__(self, file_path: str) -> None:
        """
        @brief Initializes the DataFile by checking its extension and parsing its content.
        @param file_path Path to the data file.
        @exception InvalidExtensionError if the file extension is not acceptable.
        @exception NotImplementedError if EXTENSIONS is not defined.
        """
        if not self.__class__.EXTENSIONS:
            raise NotImplementedError("EXTENSIONS must be defined.")

        if not any(file_path.endswith(ext) for ext in self.__class__.EXTENSIONS):
            raise InvalidExtensionError(
                f"Invalid file extension. Expected one of {self.__class__.EXTENSIONS}, got {file_path}"
            )

        self.container: RecordContainer = self.get_container_type()
        self.parse_file(file_path)

    def get_container_type(self) -> RecordContainer:
        """
        @brief Returns the record container type.
        @return An instance of RecordContainer.
        @exception NotImplementedError must be implemented in subclasses.
        """
        raise NotImplementedError("This method must be implemented in subclasses.")

    def parse_file(self, file_path: str) -> None:
        """
        @brief Parses the file and populates the container with records.
        @param file_path Path to the data file.
        @exception NoRecordsInDataFile if no valid records are found.
        """
        try:
            data: str = self.load_file(file_path)
            self.container.parse_records(data)
        except NoRecordsInData:
            raise NoRecordsInDataFile(f"No valid records found in file: {file_path}")

    def load_file(self, file_path: str) -> str:
        """
        @brief Loads the file content.
        @param file_path Path to the data file.
        @return The file content as a string.
        @exception NotImplementedError must be implemented by subclasses.
        """
        raise NotImplementedError("load_file must be implemented by subclasses.")

    def dump(self, output_file: str) -> None:
        """
        @brief Dumps the parsed container into a pickle file.
        @param output_file Path to the output file.
        """
        with open(output_file, "wb") as f:
            pickle.dump(self.container, f)

## ===========================================================
## Class: FASTAFile
## ===========================================================

class FASTAFile(DataFile):
    """
    @brief Data file handler for FASTA files.
    """
    EXTENSIONS = {".fa", ".fa.gz"}

    def get_container_type(self) -> FASTARecordContainer:
        """
        @brief Returns a new FASTARecordContainer.
        @return A FASTARecordContainer instance.
        """
        return FASTARecordContainer()

    def load_file(self, file_path: str) -> str:
        """
        @brief Loads FASTA file content, decompressing if necessary.
        @param file_path Path to the FASTA file.
        @return File content as a string.
        """
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt", encoding="utf-8") as file:
                return file.read()
        else:
            with open(file_path, "r", encoding="utf-8") as file:
                return file.read()

## ===========================================================
## Class: FASTAQFile
## ===========================================================

class FASTAQFile(DataFile):
    """
    @brief Data file handler for FASTQ files.
    """
    EXTENSIONS = {".fq", ".fq.gz"}

    def get_container_type(self) -> FASTAQRecordContainer:
        """
        @brief Returns a new FASTAQRecordContainer.
        @return A FASTAQRecordContainer instance.
        """
        return FASTAQRecordContainer()

    def load_file(self, file_path: str) -> str:
        """
        @brief Loads FASTQ file content, decompressing if necessary.
        @param file_path Path to the FASTQ file.
        @return File content as a string.
        """
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt", encoding="utf-8") as file:
                return file.read()
        else:
            with open(file_path, "r", encoding="utf-8") as file:
                return file.read()
