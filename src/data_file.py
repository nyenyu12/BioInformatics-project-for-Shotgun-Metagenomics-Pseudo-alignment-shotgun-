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

## @brief Exception raised when the file has an invalid extension.
class InvalidExtensionError(Exception):
    """
    @brief Exception raised when the file extension is not valid.
    """
    def __init__(self, message: str = "") -> None:
        super().__init__(message)

## @brief Exception raised when no records are found in the data file.
class NoRecordsInDataFile(Exception):
    """
    @brief Exception raised when a file does not contain any valid records.
    """
    def __init__(self, message: str = "") -> None:
        super().__init__(message)

## ===========================================================
## Class: DataFile
## ===========================================================

## @brief Abstract base class for data files.
class DataFile:
    EXTENSIONS: Optional[Set[str]] = None

    ## @brief Initializes the DataFile by checking its extension and parsing its content.
    #  @param file_path Path to the data file.
    #  @exception InvalidExtensionError if the file extension is not acceptable.
    #  @exception NotImplementedError if EXTENSIONS is not defined.
    def __init__(self, file_path: str) -> None:
        if not self.__class__.EXTENSIONS:
            raise NotImplementedError("EXTENSIONS must be defined.")

        if not any(file_path.endswith(ext) for ext in self.__class__.EXTENSIONS):
            raise InvalidExtensionError(
                f"Invalid file extension. Expected one of {self.__class__.EXTENSIONS}, got {file_path}"
            )

        self.container: RecordContainer = self.get_container_type()
        self.parse_file(file_path)

    ## @brief Returns the record container type.
    #  @return An instance of RecordContainer.
    #  @exception NotImplementedError must be implemented in subclasses.
    def get_container_type(self) -> RecordContainer:
        raise NotImplementedError("This method must be implemented in subclasses.")

    ## @brief Parses the file and populates the container with records.
    #  @param file_path Path to the data file.
    #  @exception NoRecordsInDataFile if no valid records are found.
    def parse_file(self, file_path: str) -> None:
        try:
            data: str = self.load_file(file_path)
            self.container.parse_records(data)
        except NoRecordsInData:
            raise NoRecordsInDataFile(f"No valid records found in file: {file_path}")

    ## @brief Loads the file content.
    #  @param file_path Path to the data file.
    #  @return The file content as a string.
    #  @exception NotImplementedError must be implemented by subclasses.
    def load_file(self, file_path: str) -> str:
        raise NotImplementedError("load_file must be implemented by subclasses.")

    ## @brief Dumps the parsed container into a pickle file.
    #  @param output_file Path to the output file.
    def dump(self, output_file: str) -> None:
        with open(output_file, "wb") as f:
            pickle.dump(self.container, f)

## ===========================================================
## Class: FASTAFile
## ===========================================================

## @brief Data file handler for FASTA files.
class FASTAFile(DataFile):
    EXTENSIONS = {".fa", ".fa.gz"}

    ## @brief Returns a new FASTARecordContainer.
    #  @return A FASTARecordContainer instance.
    def get_container_type(self) -> FASTARecordContainer:
        return FASTARecordContainer()

    ## @brief Loads FASTA file content, decompressing if necessary.
    #  @param file_path Path to the FASTA file.
    #  @return File content as a string.
    def load_file(self, file_path: str) -> str:
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt", encoding="utf-8") as file:
                return file.read()
        else:
            with open(file_path, "r", encoding="utf-8") as file:
                return file.read()

## ===========================================================
## Class: FASTAQFile
## ===========================================================

## @brief Data file handler for FASTQ files.
class FASTAQFile(DataFile):
    EXTENSIONS = {".fq", ".fq.gz"}

    ## @brief Returns a new FASTAQRecordContainer.
    #  @return A FASTAQRecordContainer instance.
    def get_container_type(self) -> FASTAQRecordContainer:
        return FASTAQRecordContainer()

    ## @brief Loads FASTQ file content, decompressing if necessary.
    #  @param file_path Path to the FASTQ file.
    #  @return File content as a string.
    def load_file(self, file_path: str) -> str:
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt", encoding="utf-8") as file:
                return file.read()
        else:
            with open(file_path, "r", encoding="utf-8") as file:
                return file.read()
