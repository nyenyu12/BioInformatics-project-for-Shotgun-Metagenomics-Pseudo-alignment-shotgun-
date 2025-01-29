import re
import constants
from collections import namedtuple


class NoRecordsInData(Exception):
    def __init__(self, message=""):
        super().__init__(message)


class NoRecordsInDataFile(Exception):
    def __init__(self, message=""):
        super().__init__(message)


class InvalidRecordData(Exception):
    def __init__(self, message=""):
        super().__init__(message)


# TODO: NEED TO ADD SECTION NAME to both and implement in code so we can access fields easily
Section = namedtuple("Section", ["name", "data"])
SectionSpecification = namedtuple(
    "SectionSpecifications",
    ["section_name", "section_header", "must_have_data", "section_legal_chars", "chars_to_remove"],
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
            
        return '\n'.join(self_str)
    
    def __repr__(self):
        return self.__str__()


class RecordContainer(object):

    SECTION_SPECIFICATIONS = None

    def __init__(self):
        if self.__class__.SECTION_SPECIFICATIONS == None:
            raise NotImplementedError("SECTION_SPECIFICATIONS must be defined.")

        self.__re_pattern = None
        self.create_record_re_string()
        self.__records = []

    def create_record_re_string(self):
        self.__re_pattern = []

        for (
            _,
            section_header,
            must_have_data,
            section_legal_chars,
            chars_to_remove,
        ) in self.__class__.SECTION_SPECIFICATIONS:
            self.__re_pattern.append(f"(^{re.escape(section_header)})")
            if must_have_data:
                self.__re_pattern.append(
                    f"((?:(?:{section_legal_chars})|(?:{chars_to_remove}))+?)"
                )
            else:
                self.__re_pattern.append(
                    f"((?:(?:{section_legal_chars})|(?:{chars_to_remove}))*?)"
                )

        first_section_header = self.__class__.SECTION_SPECIFICATIONS[0].section_header
        self.__re_pattern.append(f"(?:(?:^{re.escape(first_section_header)})|\\Z)")

        self.__re_pattern = "".join(self.__re_pattern)
        print(self.__re_pattern)

    def parse_records(self, data):
        for record_match in re.finditer(self.__re_pattern, data, flags=re.MULTILINE):
            self.create_record(record_match)

        if len(self.__records) == 0:
            raise NoRecordsInData

    def create_record(self, record_match):
        record_sections = []
        for specification, i in zip(
            self.__class__.SECTION_SPECIFICATIONS,
            range(1, len(record_match.groups()), 3),
        ):
            record_sections.append(
                Section(
                    name=specification.section_name,
                    data=re.sub(
                        specification.chars_to_remove,
                        "",
                        record_match.group(i + 1) or "",
                    ),
                )
            )

        self.__records.append(Record(record_sections))

    def __iter__(self):
        for record in self.__records:
            yield record


class FASTARecordContainer(RecordContainer):

    SECTION_SPECIFICATIONS = (
        SectionSpecification(
            section_name="description",
            section_header=">",
            must_have_data=True,
            section_legal_chars=".",
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
            section_legal_chars=".",
            chars_to_remove="",
        ),
        SectionSpecification(
            section_name="sequence",
            section_header="",
            must_have_data=True,
            section_legal_chars=f"[{re.escape("".join(constants.REAL_NUCLEOTIDES_CHARS))}]",
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
            section_legal_chars=f"[{re.escape("".join(constants.PHRED33_SCORES.keys()))}]",
            chars_to_remove="",
        ),
    )

    def __init__(self):
        super().__init__()
        for record in self.__records:
            if len(record["sequence"]) != len(record["quality_sequence"]):
                raise InvalidRecordData("Mismatch between nucleotide and PHRED section lengths.")

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
