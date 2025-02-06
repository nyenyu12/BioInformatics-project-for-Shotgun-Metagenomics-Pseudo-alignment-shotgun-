"""
@file data_parser.py
@brief Provides record parsing utilities for FASTA/FASTQ files.
@details Contains exceptions, record structures, and record container classes that parse
         input text into structured records.
"""

import re
import constants
from collections import namedtuple
from typing import Iterator, List, Tuple, Union, Dict, Set, Any, Optional

## Regular expression snippet used to delimit sections until the next header or EOF.
UNTIL_NEXT_HEADER_OR_EOF: str = r"(?=(?=\r?\n{section_header})|(?=(?:\r?\n)?\Z))"
## Number of characters used to show a snippet of unparsed data in error messages.
UNPARSED_SNIPPET_LEN: int = 20


## ===========================================================
## Custom Exceptions
## ===========================================================

## @brief Exception raised when no valid records are found in the data.
class NoRecordsInData(Exception):
    """
    @brief Exception for missing valid records in the input data.
    """
    def __init__(self, message: str = "No valid records found in the data.") -> None:
        super().__init__(message)


## @brief Exception raised when record data is invalid.
class InvalidRecordData(Exception):
    """
    @brief Exception for invalid record data.
    """
    def __init__(self, message: str = "") -> None:
        super().__init__(message)


## @brief Exception raised when duplicate records are found for a unique index.
class DuplicateRecordError(Exception):
    """
    @brief Exception for duplicate records (based on unique index).
    """
    def __init__(self, message: str = "Duplicate records found for the unique index.") -> None:
        super().__init__(message)


## @brief Exception raised when unparsed data remains in the input.
class UnparsedDataError(Exception):
    """
    @brief Exception for unparsed portions of the input data.
    """
    def __init__(self, message: str = "Unparsed data found in the input.") -> None:
        super().__init__(message)


## ===========================================================
## Named Tuples for Record Structure
## ===========================================================

## @brief A section of a record.
Section = namedtuple("Section", ["name", "data"])

## @brief Specification for a record section.
SectionSpecification = namedtuple(
    "SectionSpecification",
    [
        "section_name",
        "section_header",
        "must_have_data",
        "section_legal_chars",
        "chars_to_remove",
        "is_unique_index",
    ],
)


## ===========================================================
## Class: Record
## ===========================================================

## @brief Represents a parsed record.
class Record(object):
    ## @brief Initializes a Record from a list of sections.
    #  @param sections List of Section namedtuples.
    #  @exception InvalidRecordData if no sections are provided or a duplicate section name is found.
    def __init__(self, sections: List[Section]) -> None:
        if len(sections) == 0:
            raise InvalidRecordData("The data given to construct record has no sections.")
        self.identifier: str = sections[0].data
        self.__sections: Dict[str, str] = {}
        for section in sections:
            if section.name in self.__sections:
                raise InvalidRecordData(f"Section header: {section.name} has appeared twice in the given data.")
            self.__sections[section.name] = section.data

    ## @brief Retrieves the data associated with a given section key.
    #  @param key The section name.
    #  @return The corresponding section data.
    def __getitem__(self, key: str) -> str:
        return self.__sections[key]

    ## @brief Returns a string representation of the record.
    #  @return A multi-line string with section names and their data.
    def __str__(self) -> str:
        self_str: List[str] = []
        for section in self.__sections.keys():
            self_str.append(f"{section}: {self.__sections[section]}")
        return "\n".join(self_str)

    ## @brief Returns the representation of the record.
    #  @return String representation of the record.
    def __repr__(self) -> str:
        return self.__str__()


## ===========================================================
## Class: RecordContainer
## ===========================================================

## @brief Abstract base class for record containers.
class RecordContainer(object):
    ## @brief Should be defined in subclasses with the section specifications.
    SECTION_SPECIFICATIONS: Optional[Tuple[SectionSpecification, ...]] = None

    ## @brief Initializes the RecordContainer.
    #  @exception NotImplementedError if SECTION_SPECIFICATIONS is not defined.
    def __init__(self) -> None:
        if self.__class__.SECTION_SPECIFICATIONS is None:
            raise NotImplementedError("SECTION_SPECIFICATIONS must be defined.")
        self.__re_pattern: Optional[str] = None
        self.create_record_re_string()
        self._unique_index_values: Set[str] = set()
        self._records: List[Record] = []

    ## @brief Creates the regular expression string for record parsing.
    def create_record_re_string(self) -> None:
        re_parts: List[str] = []
        is_first_header: bool = True
        for (_, section_header, must_have_data, section_legal_chars, chars_to_remove, _) in self.__class__.SECTION_SPECIFICATIONS:
            if is_first_header:
                re_parts.append(f"^{re.escape(section_header)}")
                is_first_header = False
            else:
                re_parts.append(rf"\r?\n{re.escape(section_header)}")
            if must_have_data:
                re_parts.append(f"((?:[{section_legal_chars}{chars_to_remove}])+?)")
            else:
                re_parts.append(f"((?:[{section_legal_chars}{chars_to_remove}])*?)")
        first_section_header: str = re.escape(self.__class__.SECTION_SPECIFICATIONS[0].section_header)
        re_parts.append(UNTIL_NEXT_HEADER_OR_EOF.format(section_header=first_section_header))
        self.__re_pattern = "".join(re_parts)

    ## @brief Parses records from the input data string.
    #  @param data The input string to parse.
    #  @exception NoRecordsInData if no valid records are found.
    #  @exception UnparsedDataError if any unparsed data remains.
    def parse_records(self, data: str) -> None:
        parsed_indices: Set[int] = set()
        match_spans: List[Tuple[int, int]] = []
        for match in re.finditer(self.__re_pattern, data, flags=re.MULTILINE):
            if any(match.groups()):
                parsed_indices.update(range(match.start(), match.end()))
                match_spans.append((match.start(), match.end()))
                self.create_record(match.groups())
        if len(self._records) == 0:
            raise NoRecordsInData
        unparsed_data = [i for i in range(len(data)) if i not in parsed_indices and data[i].strip()]
        if unparsed_data:
            unparsed_snippet = data[min(unparsed_data): min(unparsed_data) + UNPARSED_SNIPPET_LEN]
            raise UnparsedDataError(f"Unparsed data found at index {min(unparsed_data)}: {unparsed_snippet}...")

    ## @brief Creates a Record from matched groups.
    #  @param record_match_groups Tuple of matched strings.
    def create_record(self, record_match_groups: Union[Tuple[str, Any], Tuple[str, str, str, str]]) -> None:
        sections: List[Section] = []
        for i, spec in enumerate(self.__class__.SECTION_SPECIFICATIONS):
            raw_data: str = record_match_groups[i] or ""
            cleaned_data: str = re.sub(spec.chars_to_remove, "", raw_data).strip()
            sections.append(Section(name=spec.section_name, data=cleaned_data))
            if spec.is_unique_index:
                if cleaned_data in self._unique_index_values:
                    raise DuplicateRecordError(f"Duplicate record found with unique index: {cleaned_data}")
                self._unique_index_values.add(cleaned_data)
        self._records.append(Record(sections))

    ## @brief Returns an iterator over the parsed records.
    #  @return Iterator of Record objects.
    def __iter__(self) -> Iterator[Record]:
        return iter(self._records)


## ===========================================================
## Class: FASTARecordContainer
## ===========================================================

## @brief Container for FASTA records.
class FASTARecordContainer(RecordContainer):
    SECTION_SPECIFICATIONS = (
        SectionSpecification(
            section_name="description",
            section_header=">",
            must_have_data=True,
            section_legal_chars=r"\S\t ",
            chars_to_remove="",
            is_unique_index=False,
        ),
        SectionSpecification(
            section_name="genome",
            section_header="",
            must_have_data=True,
            section_legal_chars=constants.NUCLEOTIDES_CHARS,
            chars_to_remove=r"\s",
            is_unique_index=False,
        ),
    )

    ## @brief Initializes a FASTARecordContainer.
    def __init__(self) -> None:
        super().__init__()


## ===========================================================
## Class: FASTAQRecordContainer
## ===========================================================

## @brief Container for FASTQ records.
class FASTAQRecordContainer(RecordContainer):
    SECTION_SPECIFICATIONS = (
        SectionSpecification(
            section_name="identifier",
            section_header="@",
            must_have_data=True,
            section_legal_chars=r"\S\t ",
            chars_to_remove="",
            is_unique_index=True,
        ),
        SectionSpecification(
            section_name="sequence",
            section_header="",
            must_have_data=True,
            section_legal_chars=re.escape("".join(constants.REAL_NUCLEOTIDES_CHARS)),
            chars_to_remove="",
            is_unique_index=False,
        ),
        SectionSpecification(
            section_name="space",
            section_header="+",
            must_have_data=False,
            section_legal_chars=".",
            chars_to_remove="",
            is_unique_index=False,
        ),
        SectionSpecification(
            section_name="quality_sequence",
            section_header="",
            must_have_data=True,
            section_legal_chars=re.escape("".join(constants.PHRED33_SCORES.keys())),
            chars_to_remove="",
            is_unique_index=False,
        ),
    )

    ## @brief Initializes a FASTAQRecordContainer.
    def __init__(self) -> None:
        super().__init__()

    ## @brief Parses FASTQ records and validates sequence and quality lengths.
    #  @param data Input FASTQ data as a string.
    #  @exception InvalidRecordData if sequence and quality lengths mismatch.
    def parse_records(self, data: str) -> None:
        super().parse_records(data)
        for i, record in enumerate(self):
            if len(record["sequence"]) != len(record["quality_sequence"]):
                raise InvalidRecordData(
                    f"Mismatch in record {i + 1} between nucleotide length: {len(record['sequence'])} "
                    f"and PHRED section lengths: {len(record['quality_sequence'])}"
                )
