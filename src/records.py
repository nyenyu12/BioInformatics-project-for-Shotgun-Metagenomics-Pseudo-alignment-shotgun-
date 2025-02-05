import re
import constants
from collections import namedtuple
from typing import Iterator, List, Tuple, Union, Dict, Set, Any

UNTIL_NEXT_HEADER_OR_EOF = r"(?=(?=\r?\n{section_header})|(?=(?:\r?\n)?\Z))"
UNPARSED_SNIPPET_LEN = 20


class NoRecordsInData(Exception):
    def __init__(self, message="No valid records found in the data."):
        super().__init__(message)
        

class InvalidRecordData(Exception):
    def __init__(self, message=""):
        super().__init__(message)


class DuplicateRecordError(Exception):
    def __init__(self, message="Duplicate records found for the unique index."):
        super().__init__(message)


class UnparsedDataError(Exception):
    def __init__(self, message="Unparsed data found in the input."):
        super().__init__(message)


Section = namedtuple("Section", ["name", "data"])
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


class Record(object):

    def __init__(self, sections: List[Section]):
        if len(sections) == 0:
            raise InvalidRecordData(
                "The data given to construct record has no sections."
            )

        self.identifier: str = sections[0].data
        self.__sections: Dict[str, str] = {}

        for section in sections:
            if section.name in self.__sections:
                raise InvalidRecordData(
                    f"Section header: {section.name} has appeared twice in the given data."
                )

            self.__sections[section.name] = section.data

    def __getitem__(self, key: str) -> str:
        return self.__sections[key]

    def __str__(self):
        self_str = []
        for section in self.__sections.keys():
            self_str.append(f"{section}: {self.__sections[section]}")

        return "\n".join(self_str)

    def __repr__(self):
        return self.__str__()


class RecordContainer(object):

    SECTION_SPECIFICATIONS = None

    def __init__(self):
        if self.__class__.SECTION_SPECIFICATIONS == None:
            raise NotImplementedError("SECTION_SPECIFICATIONS must be defined.")

        self.__re_pattern: str = None
        self.create_record_re_string()
        self._unique_index_values: Set[str] = set()
        self._records: List[Record] = []

    def create_record_re_string(self):
        self.__re_pattern = []
        is_first_header = True

        for (
            _,
            section_header,
            must_have_data,
            section_legal_chars,
            chars_to_remove,
            _,
        ) in self.__class__.SECTION_SPECIFICATIONS:
            if is_first_header:
                self.__re_pattern.append(f"^{re.escape(section_header)}")
                is_first_header = False
            else:
                self.__re_pattern.append(rf"\r?\n{re.escape(section_header)}")
            if must_have_data:
                self.__re_pattern.append(
                    f"((?:[{section_legal_chars}{chars_to_remove}])+?)"
                )
            else:
                self.__re_pattern.append(
                    f"((?:{section_legal_chars}{chars_to_remove}])*?)"
                )

        first_section_header = re.escape(
            self.__class__.SECTION_SPECIFICATIONS[0].section_header
        )
        self.__re_pattern.append(
            UNTIL_NEXT_HEADER_OR_EOF.format(section_header=first_section_header)
        )

        self.__re_pattern = "".join(self.__re_pattern)
        # print(self.__re_pattern)

    def parse_records(self, data: str):
        parsed_indices: Set[int] = set()
        match_spans = []

        for match in re.finditer(self.__re_pattern, data, flags=re.MULTILINE):
            if any(match.groups()):
                parsed_indices.update(range(match.start(), match.end()))
                match_spans.append((match.start(), match.end()))
                self.create_record(match.groups())

        if len(self._records) == 0:
            raise NoRecordsInData

        unparsed_data = [
            i for i in range(len(data)) if i not in parsed_indices and data[i].strip()
        ]
        if unparsed_data:
            unparsed_snippet = data[
                min(unparsed_data) : min(unparsed_data) + UNPARSED_SNIPPET_LEN
            ]
            raise UnparsedDataError(
                f"Unparsed data found at index {min(unparsed_data)}: {unparsed_snippet}..."
            )

    def create_record(self, record_match_groups: Union[Tuple[str, Any], Tuple[str, str, str, str]]):
        sections = []

        for i, spec in enumerate(self.__class__.SECTION_SPECIFICATIONS):
            raw_data = record_match_groups[i] or ""
            cleaned_data = re.sub(spec.chars_to_remove, "", raw_data).strip()

            sections.append(
                Section(
                    name=spec.section_name,
                    data=cleaned_data,
                )
            )

            if spec.is_unique_index:
                if cleaned_data in self._unique_index_values:
                    raise DuplicateRecordError(
                        f"Duplicate record found with unique index: {cleaned_data}"
                    )
                self._unique_index_values.add(cleaned_data)

        self._records.append(Record(sections))

    def __iter__(self) -> Iterator[Record]:
        for record in self._records:
            yield record


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

    def __init__(self):
        super().__init__()


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
            section_legal_chars=f"{re.escape("".join(constants.REAL_NUCLEOTIDES_CHARS))}",
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
            section_legal_chars=f"{re.escape("".join(constants.PHRED33_SCORES.keys()))}",
            chars_to_remove="",
            is_unique_index=False,
        ),
    )

    def __init__(self):
        super().__init__()

    def parse_records(self, data: str):
        super().parse_records(data)
        for i, record in enumerate(self):
            if len(record["sequence"]) != len(record["quality_sequence"]):
                raise InvalidRecordData(
                    f"Mismatch in record {i + 1} between nucleotide length: {len(record["sequence"])} and PHRED section lengths: {len(record["quality_sequence"])}"
                )

    # TODO: maybe need to seperate identifier from description in @
