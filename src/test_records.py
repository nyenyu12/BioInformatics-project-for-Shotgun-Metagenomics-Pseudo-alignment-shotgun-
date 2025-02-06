"""
@file test_records.py
@brief Unit tests for the record parsing functionality.
@details Tests for Section, SectionSpecification, Record, RecordContainer, FASTARecordContainer,
         and FASTAQRecordContainer classes.
"""

from records import (
    Section,
    SectionSpecification,
    Record,
    RecordContainer,
    InvalidRecordData,
    DuplicateRecordError,
    NoRecordsInData,
    UnparsedDataError,
    FASTARecordContainer,
    FASTAQRecordContainer,
)
import pytest
from typing import List, Any

## ===========================================================
## Section and SectionSpecification Tests
## ===========================================================

def test_section_creation() -> None:
    """
    @brief Test creation of a Section.
    @details Verifies that a Section stores its name and data correctly.
    """
    section: Section = Section(name="header", data="AGCT")
    assert section.name == "header"
    assert section.data == "AGCT"


def test_section_specification_creation() -> None:
    """
    @brief Test creation of a SectionSpecification.
    @details Verifies that SectionSpecification is created with the correct attributes.
    """
    spec: SectionSpecification = SectionSpecification(
        section_name="header",
        section_header=">",
        must_have_data=True,
        section_legal_chars="[AGCT]",
        chars_to_remove="\\s",
        is_unique_index=False,
    )
    assert spec.section_name == "header"
    assert spec.section_header == ">"
    assert spec.must_have_data is True
    assert spec.section_legal_chars == "[AGCT]"
    assert spec.chars_to_remove == "\\s"
    assert spec.is_unique_index is False

## ===========================================================
## Record Class Tests
## ===========================================================

def test_record_creation() -> None:
    """
    @brief Test creation of a Record.
    @details Verifies that the Record stores its identifier and section data correctly.
    """
    sections = [
        Section(name="header", data=">Record1"),
        Section(name="sequence", data="AGCTAGCT"),
    ]
    record: Record = Record(sections)
    assert record.identifier == ">Record1"
    assert record["header"] == ">Record1"
    assert record["sequence"] == "AGCTAGCT"


def test_record_creation_empty_sections() -> None:
    """
    @brief Test that creating a Record with no sections raises an exception.
    """
    with pytest.raises(
        InvalidRecordData, match="The data given to construct record has no sections."
    ):
        Record([])


def test_record_duplicate_sections() -> None:
    """
    @brief Test that creating a Record with duplicate section names raises an exception.
    """
    sections = [
        Section(name="header", data=">Record1"),
        Section(name="header", data=">DuplicateHeader"),
    ]
    with pytest.raises(
        InvalidRecordData,
        match="Section header: header has appeared twice in the given data.",
    ):
        Record(sections)

## ===========================================================
## RecordContainer and MockRecordContainer Tests
## ===========================================================

class MockRecordContainer(RecordContainer):
    SECTION_SPECIFICATIONS = (
        SectionSpecification(
            section_name="header",
            section_header=">",
            must_have_data=True,
            section_legal_chars="AGCT",
            chars_to_remove="",
            is_unique_index=True,
        ),
        SectionSpecification(
            section_name="sequence",
            section_header="",
            must_have_data=True,
            section_legal_chars="AGCT",
            chars_to_remove="\\s",
            is_unique_index=False,
        ),
    )


def test_record_container_re_pattern() -> None:
    """
    @brief Test that the RecordContainer builds the correct regular expression pattern.
    """
    container: MockRecordContainer = MockRecordContainer()
    expected_pattern: str = (
        r"^>((?:[AGCT])+?)" r"\r?\n((?:[AGCT\s])+?)" r"(?=(?=\r?\n>)|(?=(?:\r?\n)?\Z))"
    )
    # Access the private regex pattern using name mangling.
    assert container._RecordContainer__re_pattern == expected_pattern


def test_record_container_parse_records() -> None:
    """
    @brief Test that RecordContainer correctly parses multiple records.
    """
    container: MockRecordContainer = MockRecordContainer()
    data: str = (
        ">AGCTAGCT\nAGCTAGCT\n"
        ">AGCTTGCT\nGCGCGCGC\n"
        ">AAAAAAG\nAG A\t\tAGCT\n"
        ">AGCTGGCT\nGCGCGCGC"
    )
    container.parse_records(data)
    records: List[Record] = list(container)
    assert len(records) == 4
    assert records[0]["header"] == "AGCTAGCT"
    assert records[0]["sequence"] == "AGCTAGCT"
    assert records[1]["header"] == "AGCTTGCT"
    assert records[1]["sequence"] == "GCGCGCGC"
    # Whitespace removal test: "AG A\t\tAGCT" becomes "AGAAGCT"
    assert records[2]["header"] == "AAAAAAG"
    assert records[2]["sequence"] == "AGAAGCT"
    assert records[3]["header"] == "AGCTGGCT"
    assert records[3]["sequence"] == "GCGCGCGC"


def test_record_container_no_records() -> None:
    """
    @brief Test that parsing an empty string raises a NoRecordsInData exception.
    """
    container: MockRecordContainer = MockRecordContainer()
    with pytest.raises(NoRecordsInData, match="No valid records found in the data."):
        container.parse_records("")


def test_record_container_create_record_invalid_data() -> None:
    """
    @brief Test that invalid record data causes a parsing error.
    """
    container: MockRecordContainer = MockRecordContainer()
    invalid_data: str = ">Header1\nAGCTAGCT\nInvalidHeader"
    with pytest.raises(NoRecordsInData):
        container.parse_records(invalid_data)


def test_mock_record_container_duplicate() -> None:
    """
    @brief Test that duplicate records based on unique index raise an error.
    """
    container: MockRecordContainer = MockRecordContainer()
    data: str = ">AGCTAGCT\nAGCTAGCT\n" ">AGCTAGCT\nGCGCGCGC\n"
    with pytest.raises(DuplicateRecordError):
        container.parse_records(data)


def test_mock_record_container_valid() -> None:
    """
    @brief Test that valid data is parsed correctly by the RecordContainer.
    """
    container: MockRecordContainer = MockRecordContainer()
    data: str = ">AGCTAGCT\nAGCTAGCT\n" ">AGCTTGCT\nGCGCGCGC\n"
    container.parse_records(data)
    records: List[Record] = list(container)
    assert len(records) == 2
    assert records[0]["header"] == "AGCTAGCT"
    assert records[1]["header"] == "AGCTTGCT"


def test_mock_record_container_valid_endings() -> None:
    """
    @brief Test that valid records are parsed correctly with various line endings.
    """
    data_cases: List[str] = [
        ">AGCTAGCT\nAGCTAGCT\n",
        ">AGCTAGCT\nAGCTAGCT\r\n",
        ">AGCTAGCT\nAGCTAGCT",
    ]
    for data in data_cases:
        container: MockRecordContainer = MockRecordContainer()
        container.parse_records(data)
        records: List[Record] = list(container)
        assert len(records) == 1
        assert records[0]["header"] == "AGCTAGCT"
        assert records[0]["sequence"] == "AGCTAGCT"

## ===========================================================
## FASTARecordContainer Tests
## ===========================================================

def test_FASTARecordContainer_parse_records_valid_data() -> None:
    """
    @brief Test that FASTARecordContainer correctly parses valid FASTA data.
    """
    container: FASTARecordContainer = FASTARecordContainer()
    data: str = (
        ">Header1\nAGCTAGCT\n"
        ">Header2\nGCGCGCGC\n"
        ">Header3\nAAAAAAANN\n"
        ">Header4\nAG N \r\nC\t \nTAGCT\n"
        ">Hea\tde r5\r\nAGCTAGCT\n"
    )
    container.parse_records(data)
    records: List[Record] = list(container)
    assert len(records) == 5
    assert records[0]["description"] == "Header1"
    assert records[0]["genome"] == "AGCTAGCT"
    assert records[1]["description"] == "Header2"
    assert records[1]["genome"] == "GCGCGCGC"
    assert records[2]["description"] == "Header3"
    assert records[2]["genome"] == "AAAAAAANN"
    assert records[3]["description"] == "Header4"
    assert records[3]["genome"] == "AGNCTAGCT"
    assert records[4]["description"] == "Hea\tde r5"
    assert records[4]["genome"] == "AGCTAGCT"


def test_FASTARecordContainer_parse_records_invalid_data() -> None:
    """
    @brief Test that FASTARecordContainer raises UnparsedDataError on invalid data.
    """
    container: FASTARecordContainer = FASTARecordContainer()
    data: str = (
        ">Header1\nAGCTAGCT\n"
        ">Header2\nGCGCGCGC\n"
        ">Header3\nAAAAAAANN\n"
        ">Header4\nAG N \r\nC\t \nTAGCT\n"
        ">FalseHeader1\nAG a \r\nC\t \nTAGCT\n"
        ">Hea\tde r5\r\nAGCTAGCT\n"
    )
    with pytest.raises(UnparsedDataError):
        container.parse_records(data)

## ===========================================================
## FASTAQRecordContainer Tests
## ===========================================================

def test_FASTAQRecordContainer_valid() -> None:
    """
    @brief Test valid FASTQ parsing with various records.
    """
    container: FASTAQRecordContainer = FASTAQRecordContainer()
    data: str = (
        "@Read1\nGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read2\nATCGATCGATCGATCGATCGAA\n+\nIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read3\nGCGCGCGCGCGCGCGCGCGCGG\n+\nIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read4\nAGCTAGCTAGCTAGCTAGCTTT\n+\nIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read5\nTTTTTTTTTTTTTTTTTTTTAA\n+\nIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read6\nAGGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read7\nTTTTTTTTTTTTTTTTTGCTGCAGATCGTGGGTTTATGGATGATGTAGTGTAGAGTGAGTAGTAGTGATGGATTATGGATTGATTGAGTCAGCCG\n+\n"
        r"`1234567890-=qwertyuiop[]\asdfghjkl;'zxcvbnm,./~!@#$%^&*()_+QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>?"
        "\n"
        "@Read8\nTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAACCAGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTT\n+\n"
        r"`1234567890-=qwertyuiop[]\asdfghjkl;'zxcvbnm,./~!@#$%^&*()_+QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>?IIIIIIIIIIIIIIIIIIIIII"
        "\n"
        "@Read9\nTTTTTTTTTTTTTTTTTTTTAA\n+\nIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read10\nTTTTTTTTTTTTTTTTTTTTAA\n+\nIIIIIIIIIIIIIIIIIIIIII"
    )
    container.parse_records(data)
    records: List[Any] = list(container)
    print(records)
    assert len(records) == 10
    assert records[0]["identifier"] == "Read1"
    assert records[0]["sequence"] == "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAA"
    assert records[0]["quality_sequence"] == "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    assert records[1]["identifier"] == "Read2"
    assert records[1]["sequence"] == "ATCGATCGATCGATCGATCGAA"
    assert records[1]["quality_sequence"] == "IIIIIIIIIIIIIIIIIIIIII"
    assert records[2]["identifier"] == "Read3"
    assert records[2]["sequence"] == "GCGCGCGCGCGCGCGCGCGCGG"
    assert records[2]["quality_sequence"] == "IIIIIIIIIIIIIIIIIIIIII"
    assert records[3]["identifier"] == "Read4"
    assert records[3]["sequence"] == "AGCTAGCTAGCTAGCTAGCTTT"
    assert records[3]["quality_sequence"] == "IIIIIIIIIIIIIIIIIIIIII"
    assert records[4]["identifier"] == "Read5"
    assert records[4]["sequence"] == "TTTTTTTTTTTTTTTTTTTTAA"
    assert records[4]["quality_sequence"] == "IIIIIIIIIIIIIIIIIIIIII"
    assert records[5]["identifier"] == "Read6"
    assert records[5]["sequence"] == "AGGGGGGGGGGGGGGGGGGGGG"
    assert records[5]["quality_sequence"] == "IIIIIIIIIIIIIIIIIIIIII"
    assert records[6]["identifier"] == "Read7"
    assert (
        records[6]["sequence"]
        == "TTTTTTTTTTTTTTTTTGCTGCAGATCGTGGGTTTATGGATGATGTAGTGTAGAGTGAGTAGTAGTGATGGATTATGGATTGATTGAGTCAGCCG"
    )
    assert (
        records[6]["quality_sequence"]
        == r"`1234567890-=qwertyuiop[]\asdfghjkl;'zxcvbnm,./~!@#$%^&*()_+QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>?"
    )
    assert records[7]["identifier"] == "Read8"
    assert (
        records[7]["sequence"]
        == "TTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAACCAGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTT"
    )
    assert (
        records[7]["quality_sequence"]
        == r"`1234567890-=qwertyuiop[]\asdfghjkl;'zxcvbnm,./~!@#$%^&*()_+QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>?IIIIIIIIIIIIIIIIIIIIII"
    )
    assert records[8]["identifier"] == "Read9"
    assert records[8]["sequence"] == "TTTTTTTTTTTTTTTTTTTTAA"
    assert records[8]["quality_sequence"] == "IIIIIIIIIIIIIIIIIIIIII"
    assert records[9]["identifier"] == "Read10"
    assert records[9]["sequence"] == "TTTTTTTTTTTTTTTTTTTTAA"
    assert records[9]["quality_sequence"] == "IIIIIIIIIIIIIIIIIIIIII"


def test_FASTAQRecordContainer_mismatched_sequence_quality() -> None:
    """
    @brief Test that FASTAQRecordContainer raises an error when sequence and quality lengths mismatch.
    """
    container: FASTAQRecordContainer = FASTAQRecordContainer()
    data: str = (
        "@Read1\nGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC\n+\nIIII#V&*TDGGUIF*(IIIIIIIIIIIIIIII\n"
    )
    with pytest.raises(InvalidRecordData):
        container.parse_records(data)


def test_FASTAQRecordContainer_empty_file() -> None:
    """
    @brief Test that parsing an empty FASTQ file raises a NoRecordsInData exception.
    """
    container: FASTAQRecordContainer = FASTAQRecordContainer()
    data: str = ""
    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_invalid_nucleotides() -> None:
    """
    @brief Test that FASTAQRecordContainer raises an error when sequence contains invalid characters.
    """
    container: FASTAQRecordContainer = FASTAQRecordContainer()
    data: str = (
        "@Read1\nGGGTGATXGCCGCTGCCGATGGCGTCAAATCCCACC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )
    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_invalid_quality_chars() -> None:
    """
    @brief Test that FASTAQRecordContainer raises an error when quality string contains invalid characters.
    """
    container: FASTAQRecordContainer = FASTAQRecordContainer()
    data: str = (
        "@Read1\nGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  \n"
    )
    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_missing_space_line() -> None:
    """
    @brief Test that FASTAQRecordContainer raises an error when the '+' separator is missing.
    """
    container: FASTAQRecordContainer = FASTAQRecordContainer()
    data: str = (
        "@Read1\nGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )
    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_whitespace_handling() -> None:
    """
    @brief Test that FASTAQRecordContainer handles unexpected whitespace in the sequence.
    """
    container: FASTAQRecordContainer = FASTAQRecordContainer()
    data: str = (
        "@Read1\nGGGGG TGA TGG CCG CTG CCG ATG GCG TCA AAT CCC ACC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )
    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_parse_records_check_unparsed() -> None:
    """
    @brief Test that unparsed data in FASTQ input raises an UnparsedDataError.
    """
    container: FASTAQRecordContainer = FASTAQRecordContainer()
    data: str = (
        "@Read1\nGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read2\nATCGATCGATCGATCGATCANA\n+\nIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read3\nGCGCGCGCGCGCGCGCGCGCGG\n+\nIIIIIIIIIIIIIIIIIIIIII\n"
    )
    with pytest.raises(UnparsedDataError):
        container.parse_records(data)
