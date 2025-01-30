import gzip
import pickle
from records import FASTARecordContainer, FASTAQRecordContainer, NoRecordsInData


class NoRecordsInDataFile(Exception):
    def __init__(self, message=""):
        super().__init__(message)


class DataFile:
    EXTENSIONS = None

    def __init__(self, file_path):
        if not self.__class__.EXTENSIONS:
            raise NotImplementedError("EXTENSIONS must be defined.")

        if not any(file_path.endswith(ext) for ext in self.__class__.EXTENSIONS):
            raise ValueError(
                f"Invalid file extension. Expected one of {self.__class__.EXTENSIONS}, got {file_path}"
            )

        self.container = self.get_container_type()
        self.parse_file(file_path)

    def get_container_type(self):
        raise NotImplementedError("This method must be implemented in subclasses.")

    def parse_file(self, file_path):
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

    def get_container_type(self):
        return FASTARecordContainer()

    def load_file(self, file_path):
        """Handles loading for FASTA files, including decompression if needed."""
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt", encoding="utf-8") as file:
                return file.read()
        else:
            with open(file_path, "r", encoding="utf-8") as file:
                return file.read()

    def load(self, file_path):
        """Loads a FASTA file and parses it."""
        self.parse_file(file_path)


class FASTAQFile(DataFile):
    EXTENSIONS = {".fq", ".fq.gz"}

    def get_container_type(self):
        return FASTAQRecordContainer()

    def load_file(self, file_path):
        """Handles loading for FASTQ files, including decompression if needed."""
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt", encoding="utf-8") as file:
                return file.read()
        else:
            with open(file_path, "r", encoding="utf-8") as file:
                return file.read()

    def load(self, file_path):
        """Loads a FASTQ file and parses it."""
        self.parse_file(file_path)
