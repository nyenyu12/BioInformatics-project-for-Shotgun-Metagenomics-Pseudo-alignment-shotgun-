from data_files import *

import pytest
import re
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
        chars_to_remove="\\s"
    )
    assert spec.section_name == "header"
    assert spec.section_header == ">"
    assert spec.must_have_data is True
    assert spec.section_legal_chars == "[AGCT]"
    assert spec.chars_to_remove == "\\s"

# Test Record class
def test_record_creation():
    sections = [
        Section(name="header", data=">Record1"),
        Section(name="sequence", data="AGCTAGCT")
    ]
    record = Record(sections)
    assert record.identifier == ">Record1"
    assert record["header"] == ">Record1"
    assert record["sequence"] == "AGCTAGCT"

def test_record_creation_empty_sections():
    with pytest.raises(InvalidRecordData, match="The data given to construct record has no sections."):
        Record([])

def test_record_duplicate_sections():
    sections = [
        Section(name="header", data=">Record1"),
        Section(name="header", data=">DuplicateHeader")
    ]
    with pytest.raises(InvalidRecordData, match="Section header: header has appeared twice in the given data."):
        Record(sections)

# Test RecordContainer class
class MockRecordContainer(RecordContainer):
    SECTION_SPECIFICATIONS = (
        SectionSpecification(
            section_name="header",
            section_header=">",
            must_have_data=True,
            section_legal_chars="[AGCT]",
            chars_to_remove="\\s",
        ),
        SectionSpecification(
            section_name="sequence",
            section_header="",
            must_have_data=True,
            section_legal_chars="[AGCT]",
            chars_to_remove="\\s",
        ),
    )

def test_record_container_re_pattern():
    container = MockRecordContainer()
    expected_pattern = (
        "(^>)((?:(?:[AGCT])|(?:\\s))+?)"
        "(^)((?:(?:[AGCT])|(?:\\s))+?)"
        "(?:(?:^>)|\\Z)"
    )
    assert container._RecordContainer__re_pattern == expected_pattern

def test_record_container_parse_records():
    container = MockRecordContainer()
    data = ">AGCTAGCT\nAGCTAGCT\n>AGCTAGCT\nGCGCGCGC\n>AAAAAAG\nAG A\t\tAGCT\n>AGCTAGCT\nGCGCGCGC"
    container.parse_records(data)
    print (list(container))
    records = list(container)
    assert len(records) == 2
    assert records[0]["header"] == ">AGCTAGCT"
    assert records[0]["sequence"] == "AGCTAGCT"
    assert records[1]["header"] == ">AGCTAGCT"
    assert records[1]["sequence"] == "GCGCGCGC"
"""
def test_record_container_parse_records():
    container = MockRecordContainer()
    data = ">Header1\nAGCTAGCT\n>Header2\nGCGCGCGC"
    container.parse_records(data)

    records = list(container)
    assert len(records) == 2
    assert records[0]["header"] == ">Header1"
    assert records[0]["sequence"] == "AGCTAGCT"
    assert records[1]["header"] == ">Header2"
    assert records[1]["sequence"] == "GCGCGCGC"
"""
def test_record_container_no_records():
    container = MockRecordContainer()
    with pytest.raises(NoRecordsInData, match="No valid records found in the data."):
        container.parse_records("")

def test_record_container_create_record_invalid_data():
    container = MockRecordContainer()
    invalid_data = ">Header1\nAGCTAGCT\nInvalidHeader"
    with pytest.raises(NoRecordsInData):
        container.parse_records(invalid_data)

test_record_container_parse_records()