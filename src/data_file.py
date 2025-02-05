import gzip
import pickle
from typing import Set, Union
from records import (
    FASTARecordContainer,
    FASTAQRecordContainer,
    NoRecordsInData,
    RecordContainer,
)

class InvalidExtensionError(Exception):
    def __init__(self, message=""):
        super().__init__(message)
        
class NoRecordsInDataFile(Exception):
    def __init__(self, message=""):
        super().__init__(message)


class DataFile:
    EXTENSIONS: Union[Set[str], None] = None

    def __init__(self, file_path: str):
        if not self.__class__.EXTENSIONS:
            raise NotImplementedError("EXTENSIONS must be defined.")

        if not any(file_path.endswith(ext) for ext in self.__class__.EXTENSIONS):
            raise InvalidExtensionError(
                f"Invalid file extension. Expected one of {self.__class__.EXTENSIONS}, got {file_path}"
            )

        self.container: RecordContainer = self.get_container_type()
        self.parse_file(file_path)

    def get_container_type(self):
        raise NotImplementedError("This method must be implemented in subclasses.")

    def parse_file(self, file_path: str):
        try:
            data = self.load_file(file_path)
            self.container.parse_records(data)
        except NoRecordsInData:
            raise NoRecordsInDataFile(f"No valid records found in file: {file_path}")

    def load_file(self, file_path):
        """Hook method to be implemented by child classes for handling different file extensions."""
        raise NotImplementedError("load_file must be implemented by subclasses.")

    def dump(self, output_file):
        """Dumps the parsed data into a JSON file."""
        with open(output_file, "wb") as f:
            pickle.dump(self.container, f)


class FASTAFile(DataFile):
    EXTENSIONS = {".fa", ".fa.gz"}

    def get_container_type(self) -> FASTARecordContainer:
        return FASTARecordContainer()

    def load_file(self, file_path: str) -> str:
        """Handles loading for FASTA files, including decompression if needed."""
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt", encoding="utf-8") as file:
                return file.read()
        else:
            with open(file_path, "r", encoding="utf-8") as file:
                return file.read()


class FASTAQFile(DataFile):
    EXTENSIONS = {".fq", ".fq.gz"}

    def get_container_type(self) -> FASTAQRecordContainer:
        return FASTAQRecordContainer()

    def load_file(self, file_path: str) -> str:
        """Handles loading for FASTQ files, including decompression if needed."""
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt", encoding="utf-8") as file:
                return file.read()
        else:
            with open(file_path, "r", encoding="utf-8") as file:
                return file.read()
