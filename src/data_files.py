import re
import constants
from collections import namedtuple

UNTIL_NEXT_HEADER_OR_EOF = r"(?=(?=\r?\n{section_header})|(?=(?:\r?\n)?\Z))"


class NoRecordsInData(Exception):
    def __init__(self, message="No valid records found in the data."):
        super().__init__(message)


class NoRecordsInDataFile(Exception):
    def __init__(self, message=""):
        super().__init__(message)


class InvalidRecordData(Exception):
    def __init__(self, message=""):
        super().__init__(message)


Section = namedtuple("Section", ["name", "data"])
SectionSpecification = namedtuple(
    "SectionSpecifications",
    [
        "section_name",
        "section_header",
        "must_have_data",
        "section_legal_chars",
        "chars_to_remove",
    ],
)


class Record(object):

    def __init__(self, sections):
        if len(sections) == 0:
            raise InvalidRecordData(
                "The data given to construct record has no sections."
            )

        self.identifier = sections[0].data
        self.__sections = {}

        for section in sections:
            if section.name in self.__sections:
                raise InvalidRecordData(
                    f"Section header: {section.name} has appeared twice in the given data."
                )

            self.__sections[section.name] = section.data

    def __getitem__(self, key):
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

        self.__re_pattern = None
        self.create_record_re_string()
        self._records = []

    def create_record_re_string(self):
        self.__re_pattern = []
        is_first_header = True

        for (
            _,
            section_header,
            must_have_data,
            section_legal_chars,
            chars_to_remove,
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
        print(self.__re_pattern)

    def parse_records(self, data):
        for record_match in re.finditer(self.__re_pattern, data, flags=re.MULTILINE):
            if any(record_match.groups()):
                self.create_record(record_match.groups())

        if len(self._records) == 0:
            raise NoRecordsInData

    def create_record(self, record_match_groups):
        sections = []

        for i, spec in enumerate(self.__class__.SECTION_SPECIFICATIONS):
            # print (spec.section_name, i, record_match_groups[i])
            raw_data = record_match_groups[i] or ""
            cleaned_data = re.sub(spec.chars_to_remove, "", raw_data)

            sections.append(
                Section(
                    name=spec.section_name,
                    data=cleaned_data.strip(),
                )
            )
        self._records.append(Record(sections))

    def __iter__(self):
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
        ),
        SectionSpecification(
            section_name="genome",
            section_header="",
            must_have_data=True,
            section_legal_chars=constants.NUCLEOTIDES_CHARS,
            chars_to_remove=r"\s",
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
        ),
        SectionSpecification(
            section_name="sequence",
            section_header="",
            must_have_data=True,
            section_legal_chars=f"{re.escape("".join(constants.REAL_NUCLEOTIDES_CHARS))}",
            chars_to_remove="",
        ),
        SectionSpecification(
            section_name="space",
            section_header="+",
            must_have_data=False,
            section_legal_chars=".",
            chars_to_remove="",
        ),
        SectionSpecification(
            section_name="quality_sequence",
            section_header="",
            must_have_data=True,
            section_legal_chars=f"{re.escape("".join(constants.PHRED33_SCORES.keys()))}",
            chars_to_remove="",
        ),
    )

    def __init__(self):
        super().__init__()

    def parse_records(self, data):
        super().parse_records(data)
        for i, record in enumerate(self):
            if len(record["sequence"]) != len(record["quality_sequence"]):
                raise InvalidRecordData(
                    f"Mismatch in record {i + 1} between nucleotide length: {len(record["sequence"])} and PHRED section lengths: {len(record["quality_sequence"])}"
                )

    # TODO: maybe need to seperate identifier from description in @


class DataFile:
    EXTENSIONS = None

    def __init__(self, file_path):
        if not self.__class__.EXTENSIONS:
            raise NotImplementedError("EXTENSIONS must be defined.")

        self.container = self.get_container()
        self.parse_file(file_path)

    def get_container(self):
        raise NotImplementedError("This method must be implemented in subclasses.")

    def parse_file(self, file_path):
        try:
            with open(file_path, "r", encoding="utf-8") as file:
                data = file.read()
                self.container.parse_records(data)
        except NoRecordsInData:
            raise NoRecordsInDataFile(f"No valid records found in file: {file_path}")

    def __iter__(self):
        return iter(self.container)


class FASTAFile(DataFile):
    EXTENSIONS = {".fa", ".fa.gz"}

    def get_container(self):
        return FASTARecordContainer()


class FASTAQFile(DataFile):
    EXTENSIONS = {".fq", ".fq.gz"}

    def get_container(self):
        return FASTAQRecordContainer()
