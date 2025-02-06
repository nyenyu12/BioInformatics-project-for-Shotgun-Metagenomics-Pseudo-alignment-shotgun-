# BioSequence Project / [Elad Roichman]

## Extensions:

* **EXTSIM (Extension 4.4) – Filter Highly Similar Genomes**  
  This extension enhances the construction of the k-mer reference database by detecting and filtering out highly similar genomes. During reference building, the system computes per–genome statistics (such as unique and total k-mer counts) and applies a greedy filtering algorithm to remove nearly identical genomes. This filtering minimizes ambiguous mappings during downstream pseudo-alignment.

* **EXTQUALITY (Extension 4.1) – Pseudo-alignment with Quality Filtering**  
  This extension adds quality-based filtering to the pseudo-alignment process. In addition to the standard mapping steps, the algorithm applies three quality filters:
  1. **Read-Level Quality Filter:**  
     The mean quality of the entire read is computed (by averaging per-base quality scores). Reads whose mean quality falls below a specified threshold (`--min-read-quality` or MRQ) are completely ignored.
  2. **K-mer-Level Quality Filter:**  
     For each k-mer extracted from a read, the mean quality is calculated. K-mers with a mean quality below a given threshold (`--min-kmer-quality` or MKQ) are discarded from mapping.
  3. **Highly Redundant K-mer Filter:**  
     K-mers that map to more than a user-defined number of genomes (`--max-genomes` or MG) are excluded from the mapping process, reducing noise from overly common k-mers.
  
  When using the `dumpalign` command, the output summary is extended to include the following statistics:
  - `filtered_quality_reads` – the number of reads skipped due to low overall quality.
  - `filtered_quality_kmers` – the number of k-mers discarded due to low quality.
  - `filtered_hr_kmers` – the number of k-mers ignored because they mapped to too many genomes.

## Design:

My bioinformatics project is structured around a modular, object-oriented design that enhances maintainability and extensibility. Key design highlights include:

* **Records Module - records.py:**  
  I use classes such as `Record`, `RecordContainer`, and specialized subclasses for FASTA and FASTQ formats to encapsulate the parsing and storage of sequence data. This module enforces data integrity (e.g., unique record identifiers) and provides a uniform interface for downstream processing.

* **K-mer Module - kmer.py:**  
  The core functionality for constructing a k-mer reference database and performing pseudo-alignment is encapsulated in the `KmerReference`, `Read` and `PseudoAlignment` classes.  
  - The `KmerReference` class is responsible for extracting k-mers from genome sequences and organizing their positions using dictionaries and sets for efficient lookup.  
  - The `PseudoAlignment` class handles the mapping of sequencing reads to reference genomes by leveraging the k-mer reference.
  - The `Read` class represents a single read and handles most of the algorithmic burde, including aligning the read to a `KmerReference`  
  - **EXTQUALITY Integration:** The EXTQUALITY extension was incorporated directly into the pseudo-alignment workflow. New parameters (`--min-read-quality`, `--min-kmer-quality`, and `--max-genomes`) enable quality filtering at both the read level and the individual k-mer level. These filters are implemented as additional methods and conditions within the existing OOP structure, ensuring that the public API remains unchanged while extending functionality.

* **File I/O Module - data_file.py:**  
  The project defines an abstract `DataFile` class with concrete subclasses (`FASTAFile` and `FASTAQFile`) that manage file input/output. These classes support both plain text and gzip-compressed files. This design provides a consistent interface for reading and writing data, making it easier to add support for additional file formats in the future.

## Expanding on K-mer Module
Since the K-mer module is the heart of the project, Lets dive into it.
Since Records for FASTA records where so succint in describing genomes, I didn't define a new object for them. However for Reads I did define one since they are by far more complex objects. Here is an ascii design diagram for the K-mer module:

      +----------------------------------------+
      |    FASTARecordContainer                |
      |  (Iterable of Record objects)          |
      +------------------+---------------------+
                         |
                         v
                +-----------------+
                |     Record      |   <-- A Record object describing a FASTA record. Contains:
                |-----------------|        - identifier (e.g., "Genome1")
                | { "genome": ... }        - sections (e.g., description, genome)
                +-----------------+
                         |
                         v
                +--------------------------+
                |     KmerReference        |
                |--------------------------|
                | kmer_len: int            |
                | genomes: List[Record]    |    <-- Provided by FASTARecordContainer, the Set is a set of kmer positions in a genome
                | kmers: Dict[str,        |         e.g., {"AGCT": {Record1: {0, 4}, Record2: {10}}}
                |         Dict[Record, Set[int]]]
                | similarity_info (EXTSIM) |    <-- Optional filtering info
                +--------------------------+
                         |
                         v
                +--------------------------+
                |    Helper Functions      |
                |--------------------------|
                | extract_kmers_from_genome| <--- Yields (position, k-mer)
                +--------------------------+
                         |
                         v
                +--------------------------+
                |          Read            |
                |--------------------------|
                | identifier: str          |    <-- From FASTQ Record
                | __raw_read: str          |    <-- Raw sequence
                | kmers: Dict[str,        |         e.g., {"AGCT": ReadKmer(...)}
                |         ReadKmer]        |         where ReadKmer = (specifity, references)
                | mapping: ReadMapping     |         ReadMapping = (type, genomes_mapped_to)
                +--------------------------+
                         |
                         v
                +--------------------------+
                |      PseudoAlignment     |
                |--------------------------|
                | kmer_reference:          |
                |   KmerReference          |
                | reads: Dict[str,        |    <-- Mapping results aggregated per read
                |         { "mapping_type": ReadMappingType, "genomes_mapped_to": List[str] } |
                | Filtering counters:      |
                |  - filtered_quality_reads|
                |  - filtered_quality_kmers|
                |  - filtered_hr_kmers     |
                +--------------------------+

## Enumerations, Named Tuples, and Data Flow Summary

This is to expand on less important types in the file, those that are namedtuples or Enums.

### Enumerations

### `ReadMappingType`
Defines possible mapping statuses for a read:
- `UNMAPPED`
- `UNIQUELY_MAPPED`
- `AMBIGUOUSLY_MAPPED`

### `KmerSpecifity`
Indicates whether a k-mer is specific (appears in only one genome) or unspecific (appears in multiple genomes):
- `SPECIFIC`
- `UNSPECIFIC`

---

### Named Tuples

### `ReadKmer`
Represents information about a k-mer's occurrence in genome records.

**Fields:**
- `specifity` (`KmerSpecifity`): Indicates whether a k-mer is specific or unspecific.
- `references` (`Dict[Record, Set[int]]`): Maps each genome record to the positions where the k-mer occurs.

### `ReadMapping`
Represents the final mapping decision of a read.

**Fields:**
- `type` (`ReadMappingType`): The overall mapping status for the read.
- `genomes_mapped_to` (`List[str]`): The list of genome identifiers to which the read maps.

---

## Data Flow Summary

### **1. Input Data**
- A FASTA file is parsed by `FASTARecordContainer` to generate a list of `Record` objects.

### **2. K-mer Reference Construction**
- The `KmerReference` object is initialized with these records and a specified k-mer length.
- It extracts k-mers using the helper function and builds a nested dictionary (`self.kmers`).
- **Optional:** If **filtering is enabled (EXTSIM)**, similar genomes are filtered.

### **3. Read Processing**
- A FASTQ file is parsed to create `Read` objects.
- Each read extracts k-mers from its sequence and queries the `KmerReference` for matching records.
- The read constructs a dictionary of `ReadKmer` tuples that indicate the specificity of each k-mer.
- The read then decides its overall mapping (`ReadMapping`) based on the collected k-mer data.

### **4. Aggregation**
- The `PseudoAlignment` object collects these reads, aggregates their mapping results, and produces a summary.
- **If EXTQUALITY is applied**, additional filtering statistics are included.


## Expanding on Records Module
Since the parsing code is very generic and may be hard to comprehend, then I'll expand on it as well.

## **1. Named Tuples for Record Structure**
This is very important to understand BEFORE going into the flow unlike in the K-mer module. Records are designed 
to handle a generic record with the only expectation being that each section's data can be parsed as a regular expression (i.e. is regular).

### **`Section`**
Represents a section of a record.

- **Fields:**
  - `name`: Name of the section.
  - `data`: Content of the section.

### **`SectionSpecification`**
Defines the structure of a record section.

- **Fields:**
  - `section_name`: Name of the section.
  - `section_header`: The delimiter marking the start of the section.
  - `must_have_data`: Whether the section is required.
  - `section_legal_chars`: Allowed characters in this section.
  - `chars_to_remove`: Characters to remove from this section.
  - `is_unique_index`: Whether this section's data must be unique.

## **2. Core Record Classes**
Here I highlight key classes in the module and thier important attributes and methods .

### **`Record`**
Represents a structured data entry parsed from FASTA/FASTQ files. A `Record` consists of multiple named sections.

- **Attributes:**
  - `identifier`: The unique identifier for the record.
  - `__sections`: A dictionary mapping section names to their corresponding data.

### **`RecordContainer` (Abstract Base Class)**
Handles parsing input data into `Record` objects.

- **Attributes:**
  - `SECTION_SPECIFICATIONS`: Defines expected record sections.
  - `_records`: Stores parsed `Record` objects.

- **Methods:**
  - `parse_records(data)`: Extracts records from input text.
  - `create_record(match_groups)`: Constructs `Record` instances from parsed data.
  - `__iter__()`: Iterates over stored `Record` objects.

`SECTION_SPECIFICATIONS` is the most important connection etween `RecordContainer`, `Record`, and Child Classes of `RecordContainer`. It defines the structure of records by specifying the expected sections and their properties. This specification is crucial because it ensures that:

 **1. `RecordContainer` Can Parse Input Data into Sections**
- `SECTION_SPECIFICATIONS` provides the necessary details to extract meaningful sections from the input text.
- It enables `RecordContainer` to construct a **regular expression** (`self.__re_pattern`) that is used to split input data into structured components.

 **2. Consistency Between `RecordContainer` and `Record` Instances**
- Each parsed record is converted into a `Record` object, where sections are stored in a dictionary (`self.__sections`).
- Since `SECTION_SPECIFICATIONS` dictates which sections must exist, it ensures that every `Record` created from a `RecordContainer` follows the same structure.

 **3. Child Classes of `RecordContainer` Can Define Custom Parsing Rules**
- `FASTARecordContainer` and `FASTAQRecordContainer` inherit from `RecordContainer` but define their own `SECTION_SPECIFICATIONS`, allowing them to handle different file formats while keeping a **uniform interface for parsing**.


## **3. Data Flow Summary**

### **ASCII Diagram of Object Encapsulation and Data Flow**

```lua
      +--------------------------------------+
      |        RecordContainer               |
      |  (Abstract base class)               |
      +------------------+-------------------+
                         |
                         v
      +--------------------------------------+
      |           FASTARecordContainer       |
      |    (Parses FASTA format records)     |
      +--------------------------------------+
                         |
                         v
      +--------------------------------------+
      |           FASTAQRecordContainer      |
      |    (Parses FASTQ format records)     |
      +--------------------------------------+
                         |
                         v
      +--------------------------------------+
      |             Record                   |
      |  (Holds parsed sections from input)  |
      +--------------------------------------+
                         |
                         v
      +--------------------------------------+
      |        Section & SectionSpecification|
      |  (Define record structure details)   |
      +--------------------------------------+
```


**Command-Line Interface (CLI):**
---
  The `main.py` file serves as the command-line entry point, where arguments are parsed and dispatched to the appropriate modules. This design decouples the core bioinformatics logic from the user interface, allowing for flexible integration and testing.

