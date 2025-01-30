from records import *
import pytest
from collections import namedtuple

# Assuming the classes Section, SectionSpecification, Record, RecordContainer are already defined


# Test Section and SectionSpecification
def test_section_creation():
    section = Section(name="header", data="AGCT")
    assert section.name == "header"
    assert section.data == "AGCT"


def test_section_specification_creation():
    spec = SectionSpecification(
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
    assert spec.is_unique_index == False


# Test Record class
def test_record_creation():
    sections = [
        Section(name="header", data=">Record1"),
        Section(name="sequence", data="AGCTAGCT"),
    ]
    record = Record(sections)
    assert record.identifier == ">Record1"
    assert record["header"] == ">Record1"
    assert record["sequence"] == "AGCTAGCT"


def test_record_creation_empty_sections():
    with pytest.raises(
        InvalidRecordData, match="The data given to construct record has no sections."
    ):
        Record([])


def test_record_duplicate_sections():
    sections = [
        Section(name="header", data=">Record1"),
        Section(name="header", data=">DuplicateHeader"),
    ]
    with pytest.raises(
        InvalidRecordData,
        match="Section header: header has appeared twice in the given data.",
    ):
        Record(sections)


# Test RecordContainer class
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


def test_record_container_re_pattern():
    container = MockRecordContainer()
    expected_pattern = (
        r"^>((?:[AGCT])+?)" r"\r?\n((?:[AGCT\s])+?)" r"(?=(?=\r?\n>)|(?=(?:\r?\n)?\Z))"
    )
    assert container._RecordContainer__re_pattern == expected_pattern


def test_record_container_parse_records():
    container = MockRecordContainer()
    data = (
        ">AGCTAGCT\nAGCTAGCT\n"
        ">AGCTTGCT\nGCGCGCGC\n"
        ">AAAAAAG\nAG A\t\tAGCT\n"
        ">AGCTGGCT\nGCGCGCGC"
    )
    container.parse_records(data)
    records = list(container)
    assert len(records) == 4
    assert records[0]["header"] == "AGCTAGCT"
    assert records[0]["sequence"] == "AGCTAGCT"
    assert records[1]["header"] == "AGCTTGCT"
    assert records[1]["sequence"] == "GCGCGCGC"
    assert records[2]["header"] == "AAAAAAG"
    assert records[2]["sequence"] == "AGAAGCT"
    assert records[3]["header"] == "AGCTGGCT"
    assert records[3]["sequence"] == "GCGCGCGC"


def test_record_container_no_records():
    container = MockRecordContainer()
    with pytest.raises(NoRecordsInData, match="No valid records found in the data."):
        container.parse_records("")


def test_record_container_create_record_invalid_data():
    container = MockRecordContainer()
    invalid_data = ">Header1\nAGCTAGCT\nInvalidHeader"
    with pytest.raises(NoRecordsInData):
        container.parse_records(invalid_data)


def test_mock_record_container_duplicate():
    container = MockRecordContainer()
    data = ">AGCTAGCT\nAGCTAGCT\n" ">AGCTAGCT\nGCGCGCGC\n"
    with pytest.raises(DuplicateRecordError):
        container.parse_records(data)


def test_mock_record_container_valid():
    container = MockRecordContainer()
    data = ">AGCTAGCT\nAGCTAGCT\n" ">AGCTTGCT\nGCGCGCGC\n"
    container.parse_records(data)
    records = list(container)
    assert len(records) == 2
    assert records[0]["header"] == "AGCTAGCT"
    assert records[1]["header"] == "AGCTTGCT"


def test_mock_record_container_duplicate():
    container = MockRecordContainer()
    data = ">AGCTAGCT\nAGCTAGCT\n" ">AGCTAGCT\nGCGCGCGC\n"
    with pytest.raises(DuplicateRecordError):
        container.parse_records(data)


def test_mock_record_container_valid():
    container = MockRecordContainer()
    data = ">AGCTAGCT\nAGCTAGCT\n" ">AGCTTGCT\nGCGCGCGC\n"
    container.parse_records(data)
    records = list(container)
    assert len(records) == 2
    assert records[0]["header"] == "AGCTAGCT"
    assert records[1]["header"] == "AGCTTGCT"


def test_mock_record_container_valid_endings():
    data_cases = [
        ">AGCTAGCT\nAGCTAGCT\n",
        ">AGCTAGCT\nAGCTAGCT\r\n",
        ">AGCTAGCT\nAGCTAGCT",
    ]

    for data in data_cases:
        container = MockRecordContainer()
        container.parse_records(data)
        records = list(container)
        assert len(records) == 1
        assert records[0]["header"] == "AGCTAGCT"
        assert records[0]["sequence"] == "AGCTAGCT"


def test_FASTARecordContainer_parse_records_valid_data():
    container = FASTARecordContainer()
    data = (
        ">Header1\nAGCTAGCT\n"
        ">Header2\nGCGCGCGC\n"
        ">Header3\nAAAAAAANN\n"
        ">Header4\nAG N \r\nC\t \nTAGCT\n"
        ">Hea\tde r5\r\nAGCTAGCT\n"
    )
    container.parse_records(data)

    records = list(container)
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


def test_FASTARecordContainer_parse_records_invalid_data():
    container = FASTARecordContainer()
    data = (
        ">Header1\nAGCTAGCT\n"
        ">Header2\nGCGCGCGC\n"
        ">Header3\nAAAAAAANN\n"
        ">Header4\nAG N \r\nC\t \nTAGCT\n"
        ">FalseHeader1\nAG a \r\nC\t \nTAGCT\n"
        ">Hea\tde r5\r\nAGCTAGCT\n"
    )
    with pytest.raises(UnparsedDataError):
        container.parse_records(data)


def test_FASTAQRecordContainer_valid():
    """Test valid FASTQ parsing with different cases."""
    container = FASTAQRecordContainer()
    data = (
        "@Read1\n"
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAA\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read2\n"
        "ATCGATCGATCGATCGATCGAA\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIII\n"
        "@Read3\n"
        "GCGCGCGCGCGCGCGCGCGCGG\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIII\n"
        "@Read4\n"
        "AGCTAGCTAGCTAGCTAGCTTT\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIII\n"
        "@Read5\n"
        "TTTTTTTTTTTTTTTTTTTTAA\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIII\n"
        "@Read6\n"
        "AGGGGGGGGGGGGGGGGGGGGG\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIII\n"
        "@Read7\n"
        "TTTTTTTTTTTTTTTTTGCTGCAGATCGTGGGTTTATGGATGATGTAGTGTAGAGTGAGTAGTAGTGATGGATTATGGATTGATTGAGTCAGCCG\n"
        "+\n"
        r"`1234567890-=qwertyuiop[]\asdfghjkl;'zxcvbnm,./~!@#$%^&*()_+QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>?"
        "\n"
        "@Read8\n"
        "TTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAACCAGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTT\n"
        "+\n"
        r"`1234567890-=qwertyuiop[]\asdfghjkl;'zxcvbnm,./~!@#$%^&*()_+QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>?IIIIIIIIIIIIIIIIIIIIII"
        "\n"
        "@Read9\n"
        "TTTTTTTTTTTTTTTTTTTTAA\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIII\n"
        "@Read10\n"
        "TTTTTTTTTTTTTTTTTTTTAA\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIII"  # Should work with \r\n,\n ending and without
    )

    container.parse_records(data)
    records = list(container)

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


def test_FASTAQRecordContainer_mismatched_sequence_quality():
    """Test case where sequence and quality lengths mismatch."""
    container = FASTAQRecordContainer()
    data = (
        "@Read1\n"
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC\n"
        "+\n"
        r"IIII#V&*TDGGUIF*(IIIIIIIIIIIIIIII\n"  # One char short
    )

    with pytest.raises(InvalidRecordData):
        container.parse_records(data)


def test_FASTAQRecordContainer_empty_file():
    """Test case for empty FASTQ file."""
    container = FASTAQRecordContainer()
    data = ""

    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_invalid_nucleotides():
    """Test case where sequence contains invalid characters."""
    container = FASTAQRecordContainer()
    data = (
        "@Read1\n"
        "GGGTGATXGCCGCTGCCGATGGCGTCAAATCCCACC\n"  # 'X' is not a valid nucleotide
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )

    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_invalid_quality_chars():
    """Test case where quality string contains invalid characters."""
    container = FASTAQRecordContainer()
    data = (
        "@Read1\n"
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  \n"  # ' ' is not a valid Phred33 char
    )

    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_missing_space_line():
    """Test case where the '+' space separator line is missing."""
    container = FASTAQRecordContainer()
    data = (
        "@Read1\n"
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"  # Missing '+'
    )

    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_whitespace_handling():
    """Test case where sequence and quality strings contain unexpected whitespace."""
    container = FASTAQRecordContainer()
    data = (
        "@Read1\n"
        "GGGGG TGA TGG CCG CTG CCG ATG GCG TCA AAT CCC ACC\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )

    with pytest.raises(NoRecordsInData):
        container.parse_records(data)


def test_FASTAQRecordContainer_parse_records_check_unparsed():
    container = FASTAQRecordContainer()
    data = (
        "@Read1\n"
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAA\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        "@Read2\n"
        "ATCGATCGATCGATCGATCANA\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIII\n"
        "@Read3\n"
        "GCGCGCGCGCGCGCGCGCGCGG\n"
        "+\n"
        "IIIIIIIIIIIIIIIIIIIIII\n"
    )
    with pytest.raises(UnparsedDataError):
        container.parse_records(data)


test_FASTAQRecordContainer_valid()
